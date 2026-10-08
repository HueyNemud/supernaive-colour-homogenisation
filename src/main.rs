use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use rayon::prelude::*;

use homog::background::Lighting;
use homog::calibrate::{self, Proposal};
use homog::config::Config;
use homog::hints::Hints;
use homog::io::{Format, OutputOptions};
use homog::layout::{Atlas, Seams};
use homog::{process_file, Options};

/// Not-So-Naive Map Colour Homogenisation: makes the paper of scanned map sheets white,
/// evenly across each sheet and across adjacent sheets, while keeping the inks and washes.
/// Settings come from an optional `homog.toml` (-c), each overridden by the corresponding
/// option. Georeferencing is kept.
#[derive(Parser)]
#[command(version, args_conflicts_with_subcommands = true, subcommand_negates_reqs = true)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    process: ProcessArgs,
}

#[derive(Subcommand)]
enum Command {
    /// Analyse sheets without writing corrected images: contact sheet of what was taken
    /// for paper, report with warnings, and a proposed homog.toml
    Calibrate(CalibrateArgs),
}

#[derive(Args)]
struct Common {
    /// Input images (RGB/RGBA, 8 or 16 bits: JPEG, TIFF/GeoTIFF, PNG...)
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Configuration file (homog.toml)
    #[arg(short, long)]
    config: Option<PathBuf>,

    /// Target background colour, as r,g,b [default: 255,255,255]
    #[arg(short, long, value_parser = parse_rgb)]
    target: Option<[u8; 3]>,

    /// Unevenness of the lighting, i.e. flexibility of the paper surface [default: normal]
    #[arg(short, long, value_enum)]
    lighting: Option<Lighting>,

    /// Directory of imagettes of bare paper
    #[arg(long)]
    paper: Option<PathBuf>,

    /// Directory of imagettes of colours that must not become white
    #[arg(long)]
    keep: Option<PathBuf>,

    /// Directory of QGIS georeferencer GCP files (`<image file name>.points`) giving the
    /// layout of the sheets: seams with adjacent sheets are corrected (images are not resampled)
    #[arg(short, long)]
    gcp_dir: Option<PathBuf>,

    /// Ignore embedded ICC profiles (assume sRGB)
    #[arg(long)]
    ignore_icc: bool,

    /// Number of threads (default: all cores)
    #[arg(short, long)]
    jobs: Option<usize>,

    /// Print per-image details
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Args)]
struct ProcessArgs {
    #[command(flatten)]
    common: Common,

    /// Output directory [default: output]
    #[arg(short, long)]
    out_dir: Option<PathBuf>,

    /// Output format [default: auto, same family as the input]
    #[arg(short, long, value_enum)]
    format: Option<Format>,

    /// JPEG quality [default: 95]
    #[arg(long, value_parser = clap::value_parser!(u8).range(1..=100))]
    jpeg_quality: Option<u8>,

    /// TIFF compression: NONE, LZW, DEFLATE, ZSTD... [default: DEFLATE]
    #[arg(long)]
    compress: Option<String>,

    /// Write debug images (paper surface, gains, seam measurements...) to <out-dir>/debug/<image>/
    #[arg(long)]
    debug: bool,
}

#[derive(Args)]
struct CalibrateArgs {
    #[command(flatten)]
    common: Common,

    /// Output directory of the calibration report
    #[arg(short, long, default_value = "calib")]
    out_dir: PathBuf,
}

fn parse_rgb(s: &str) -> Result<[u8; 3], String> {
    let v: Vec<u8> = s
        .split(',')
        .map(|c| c.trim().parse::<u8>().map_err(|_| format!("invalid component '{c}' (expected 0-255)")))
        .collect::<Result<_, _>>()?;
    v.try_into().map_err(|_| format!("expected r,g,b, got '{s}'"))
}

/// Settings shared by both commands, after merging the configuration file and the options.
struct Settings {
    config: Config,
    target: [u8; 3],
    lighting: Lighting,
    icc: bool,
    paper: Option<PathBuf>,
    keep: Option<PathBuf>,
    gcp_dir: Option<PathBuf>,
}

impl Settings {
    fn new(c: &Common) -> Result<Self> {
        let config = c.config.as_deref().map(Config::load).transpose()?.unwrap_or_default();
        Ok(Self {
            target: c.target.or(config.target).unwrap_or([255, 255, 255]),
            lighting: c.lighting.or(config.lighting).unwrap_or_default(),
            icc: !c.ignore_icc && config.icc.unwrap_or(true),
            paper: c.paper.clone().or(config.hints.paper.clone()),
            keep: c.keep.clone().or(config.hints.keep.clone()),
            gcp_dir: c.gcp_dir.clone().or(config.seams.gcp_dir.clone()),
            config,
        })
    }

    fn options(&self, output: OutputOptions, debug: Option<PathBuf>) -> Result<Options> {
        Ok(Options {
            target: self.target,
            lighting: self.lighting,
            use_icc: self.icc,
            hints: Hints::load(self.paper.as_deref(), self.keep.as_deref())?,
            output,
            debug,
        })
    }

    /// Layout and shared edges of each input, from the GCP directory.
    fn seams(&self, files: &[&PathBuf]) -> Result<Vec<Option<Seams>>> {
        let Some(dir) = &self.gcp_dir else { return Ok(vec![None; files.len()]) };
        let atlas = Atlas::load(dir)?;
        if atlas.is_empty() {
            bail!("no .points file in {}", dir.display());
        }
        Ok(files
            .iter()
            .map(|f| {
                let name = f.file_name().and_then(|n| n.to_str()).unwrap_or_default();
                let s = atlas.seams(name);
                if s.is_none() {
                    eprintln!("warning: no GCPs for {name}, its seams are not corrected");
                }
                s
            })
            .collect())
    }
}

fn setup(c: &Common) -> Result<()> {
    if let Some(j) = c.jobs {
        rayon::ThreadPoolBuilder::new().num_threads(j).build_global()?;
    }
    Ok(())
}

fn process(args: ProcessArgs) -> Result<bool> {
    let c = &args.common;
    setup(c)?;
    let s = Settings::new(c)?;
    let out = &s.config.output;
    let out_dir = args.out_dir.clone().or(out.dir.clone()).unwrap_or_else(|| PathBuf::from("output"));
    let format = args.format.or(out.format).unwrap_or(Format::Auto);
    let output = OutputOptions {
        jpeg_quality: args.jpeg_quality.or(out.jpeg_quality).unwrap_or(95),
        tiff_compress: args.compress.clone().or(out.compress.clone()).unwrap_or_else(|| "DEFLATE".into()),
    };
    std::fs::create_dir_all(&out_dir).with_context(|| format!("cannot create {}", out_dir.display()))?;

    // Resolve output paths up front to detect collisions and in-place overwrites.
    let mut seen: HashMap<PathBuf, &PathBuf> = HashMap::new();
    let mut jobs = Vec::new();
    for f in &c.files {
        let ext = format.for_input(f).extension();
        let stem = f.file_stem().with_context(|| format!("invalid file name {}", f.display()))?;
        let dst = out_dir.join(stem).with_extension(ext);
        if let Some(prev) = seen.insert(dst.clone(), f) {
            bail!("{} and {} would both be written to {}", prev.display(), f.display(), dst.display());
        }
        if f.canonicalize().ok() == dst.canonicalize().ok() && dst.exists() {
            bail!("refusing to overwrite the input {}", f.display());
        }
        jobs.push((f, dst));
    }
    let seams = s.seams(&jobs.iter().map(|(f, _)| *f).collect::<Vec<_>>())?;
    let opts = s.options(output, args.debug.then(|| out_dir.join("debug")))?;

    let results: Vec<bool> = jobs
        .par_iter()
        .zip(&seams)
        .map(|((input, output), seams)| match process_file(input, output, &opts, seams.as_ref()) {
            Ok(r) => {
                if c.verbose {
                    let p = r.paper_lab;
                    eprintln!(
                        "{} -> {} [{}x{}{}] paper Lab ({:.1}, {:.1}, {:.1}), surface {} L {:.1}..{:.1}, {} seam(s), {:.3}s (thumbnail {:.3}, analysis {:.3}, pass {:.3}, finalise {:.3})",
                        input.display(),
                        output.display(),
                        r.width,
                        r.height,
                        if r.icc { ", ICC" } else { "" },
                        p[0], p[1], p[2],
                        r.degree.map_or("flat".into(), |d| format!("degree {d}")),
                        r.paper_l_range.0,
                        r.paper_l_range.1,
                        r.seam_edges,
                        r.timings.iter().sum::<f64>(),
                        r.timings[0], r.timings[1], r.timings[2], r.timings[3]
                    );
                }
                true
            }
            Err(e) => {
                eprintln!("error: {e:#}");
                false
            }
        })
        .collect();
    Ok(results.iter().all(|&ok| ok))
}

fn calibrate(args: CalibrateArgs) -> Result<bool> {
    let c = &args.common;
    setup(c)?;
    let s = Settings::new(c)?;
    let files: Vec<&PathBuf> = c.files.iter().collect();
    let seams = s.seams(&files)?;
    let opts = s.options(OutputOptions { jpeg_quality: 95, tiff_compress: "DEFLATE".into() }, None)?;
    let reports: Vec<_> = files
        .par_iter()
        .zip(&seams)
        .map(|(f, seams)| calibrate::calibrate_sheet(f, &opts, seams.as_ref()))
        .collect();
    let mut sheets = Vec::new();
    let mut ok = true;
    for r in reports {
        match r {
            Ok(s) => sheets.push(s),
            Err(e) => {
                eprintln!("error: {e:#}");
                ok = false;
            }
        }
    }
    let lighting = s.lighting.to_possible_value().map(|v| v.get_name().to_owned()).unwrap_or_default();
    let proposal = Proposal { target: s.target, lighting, icc: s.icc, gcp_dir: s.gcp_dir.clone(), paper_hints: s.paper.clone(), keep_hints: s.keep.clone() };
    calibrate::write(&args.out_dir, &sheets, &proposal)?;
    if c.verbose {
        eprintln!("{} sheet(s) analysed, report in {}", sheets.len(), Path::new(&args.out_dir).join("calibrate.txt").display());
    }
    Ok(ok)
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Some(Command::Calibrate(args)) => calibrate(args),
        None => process(cli.process),
    };
    match result {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
