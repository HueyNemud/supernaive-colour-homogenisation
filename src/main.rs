use std::collections::HashMap;
use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;

use homog::io::{init_gdal, Format, OutputOptions};
use homog::layout::Atlas;
use homog::{process_file, Options};

/// Homogenise the background of scanned map sheets: the paper colour of each
/// sheet (estimated locally) is mapped to a target white by Bradford chromatic
/// adaptation, preserving the other colours. Georeferencing is kept.
#[derive(Parser)]
#[command(version)]
struct Cli {
    /// Input images (RGB/RGBA, 8 or 16 bits: JPEG, TIFF/GeoTIFF, PNG...)
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Target background colour, as r,g,b
    #[arg(short, long, default_value = "255,255,255", value_parser = parse_rgb)]
    target: [u8; 3],

    /// Output directory
    #[arg(short, long, default_value = "output")]
    out_dir: PathBuf,

    /// Output format
    #[arg(short, long, value_enum, default_value_t = Format::Auto)]
    format: Format,

    /// JPEG quality (JPEG output only)
    #[arg(long, default_value_t = 95, value_parser = clap::value_parser!(u8).range(1..=100))]
    jpeg_quality: u8,

    /// TIFF compression (NONE, LZW, DEFLATE, ZSTD...)
    #[arg(long, default_value = "DEFLATE")]
    compress: String,

    /// Border fraction ignored for the global paper estimate
    #[arg(short, long, default_value_t = 0.05)]
    margin: f64,

    /// Use a single paper colour per sheet instead of a paper surface
    #[arg(long)]
    no_flat_field: bool,

    /// Degree of the polynomial paper surface (higher follows the paper more closely,
    /// but starts following large washes)
    #[arg(long, default_value_t = 3, value_parser = clap::value_parser!(u8).range(1..=6))]
    degree: u8,

    /// Directory of QGIS georeferencer GCP files (`<image file name>.points`) giving the layout
    /// of the sheets: seams with adjacent sheets are then corrected. Images are not resampled.
    #[arg(short, long)]
    gcp_dir: Option<PathBuf>,

    /// Ignore embedded ICC profiles (assume sRGB)
    #[arg(long)]
    ignore_icc: bool,

    /// Write debug images (paper surface, gains, seam measurements...) to <out-dir>/debug/<image>/
    #[arg(long)]
    debug: bool,

    /// Number of threads (default: all cores)
    #[arg(short, long)]
    jobs: Option<usize>,

    /// Print per-image details
    #[arg(short, long)]
    verbose: bool,
}

fn parse_rgb(s: &str) -> Result<[u8; 3], String> {
    let v: Vec<u8> = s
        .split(',')
        .map(|c| c.trim().parse::<u8>().map_err(|_| format!("invalid component '{c}' (expected 0-255)")))
        .collect::<Result<_, _>>()?;
    v.try_into().map_err(|_| format!("expected r,g,b, got '{s}'"))
}

fn run(cli: Cli) -> Result<bool> {
    if let Some(j) = cli.jobs {
        rayon::ThreadPoolBuilder::new().num_threads(j).build_global()?;
    }
    init_gdal();
    std::fs::create_dir_all(&cli.out_dir).with_context(|| format!("cannot create {}", cli.out_dir.display()))?;

    // Resolve output paths up front to detect collisions and in-place overwrites.
    let mut seen: HashMap<PathBuf, &PathBuf> = HashMap::new();
    let mut jobs = Vec::new();
    for f in &cli.files {
        let ext = match cli.format {
            Format::Auto => match f.extension().and_then(|e| e.to_str()).map(str::to_lowercase).as_deref() {
                Some("jpg" | "jpeg") => "jpg",
                Some("png") => "png",
                _ => "tif",
            },
            fmt => fmt.extension(),
        };
        let stem = f.file_stem().with_context(|| format!("invalid file name {}", f.display()))?;
        let out = cli.out_dir.join(stem).with_extension(ext);
        if let Some(prev) = seen.insert(out.clone(), f) {
            bail!("{} and {} would both be written to {}", prev.display(), f.display(), out.display());
        }
        if f.canonicalize().ok() == out.canonicalize().ok() && out.exists() {
            bail!("refusing to overwrite the input {}", f.display());
        }
        jobs.push((f, out));
    }

    let atlas = cli.gcp_dir.as_deref().map(Atlas::load).transpose()?;
    if let Some(a) = &atlas
        && a.is_empty()
    {
        bail!("no .points file in {}", cli.gcp_dir.as_ref().unwrap().display());
    }
    let seams: Vec<_> = jobs
        .iter()
        .map(|(f, _)| {
            let name = f.file_name().and_then(|n| n.to_str()).unwrap_or_default();
            let s = atlas.as_ref().and_then(|a| a.seams(name));
            if atlas.is_some() && s.is_none() {
                eprintln!("warning: no GCPs for {name}, its seams are not corrected");
            }
            s
        })
        .collect();

    let opts = Options {
        target: cli.target,
        margin: cli.margin,
        flat_field: !cli.no_flat_field,
        degree: cli.degree as usize,
        thumb_side: 400,
        use_icc: !cli.ignore_icc,
        output: OutputOptions { format: cli.format, jpeg_quality: cli.jpeg_quality, tiff_compress: cli.compress },
        debug: cli.debug.then(|| cli.out_dir.join("debug")),
    };

    let results: Vec<bool> = jobs
        .par_iter()
        .zip(&seams)
        .map(|((input, output), seams)| match process_file(input, output, &opts, seams.as_ref()) {
            Ok(r) => {
                if cli.verbose {
                    let p = r.paper_lab;
                    eprintln!(
                        "{} -> {} [{}x{}{}] paper Lab ({:.1}, {:.1}, {:.1}), local L {:.1}..{:.1}, {} seam(s), {:.3}s (thumbnail {:.3}, analysis {:.3}, pass {:.3}, finalise {:.3})",
                        input.display(),
                        output.display(),
                        r.width,
                        r.height,
                        if r.icc { ", ICC" } else { "" },
                        p[0], p[1], p[2],
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

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
