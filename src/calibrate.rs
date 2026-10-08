//! `homog calibrate`: analyse a batch of sheets without writing corrected images, and report
//! what was taken for paper so that hints can be added where the automatic estimate fails.
//!
//! Writes to the output directory:
//! - `contact_sheet.png`: one row per sheet (order of `calibrate.txt`): original, role of the
//!   pixels in the paper fit (green = bare paper, orange = lightness only, grey = ignored,
//!   dark red = outside the sheet), corrected result;
//! - `calibrate.txt`: paper colour of the batch and, per sheet, its paper colour, the share
//!   of bare paper found, the surface degree and warnings (paper far from the batch, little
//!   bare paper found);
//! - `homog.toml`: proposed configuration.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::background::median;
use crate::color::Vec3;
use crate::debug::{result_image, role_colour, srgb8};
use crate::io::{write_png, Input};
use crate::layout::Seams;
use crate::{analyse, Options};

/// Height of the contact sheet rows
const ROW_HEIGHT: usize = 220;
const GAP: usize = 6;

pub struct SheetReport {
    pub name: String,
    pub paper: Vec3,
    pub degree: Option<usize>,
    pub paper_fraction: f64,
    pub seam_edges: usize,
    row: (usize, usize, Vec<[u8; 3]>),
}

/// Bilinear resize of an 8-bit RGB image to the given height.
fn resize(w: usize, h: usize, px: &[[u8; 3]], nh: usize) -> (usize, Vec<[u8; 3]>) {
    let nw = (w * nh).div_ceil(h).max(1);
    let out = (0..nw * nh)
        .map(|i| {
            let x = ((i % nw) as f64 + 0.5) * w as f64 / nw as f64 - 0.5;
            let y = ((i / nw) as f64 + 0.5) * h as f64 / nh as f64 - 0.5;
            let (x, y) = (x.clamp(0.0, (w - 1) as f64), y.clamp(0.0, (h - 1) as f64));
            let (x0, y0) = (x as usize, y as usize);
            let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
            let (fx, fy) = (x - x0 as f64, y - y0 as f64);
            [0, 1, 2].map(|c| {
                let at = |x: usize, y: usize| px[y * w + x][c] as f64;
                let top = at(x0, y0) + (at(x1, y0) - at(x0, y0)) * fx;
                let bottom = at(x0, y1) + (at(x1, y1) - at(x0, y1)) * fx;
                (top + (bottom - top) * fy).round() as u8
            })
        })
        .collect();
    (nw, out)
}

/// Images side by side, separated by white gaps.
fn hstack(images: &[(usize, Vec<[u8; 3]>)], h: usize) -> (usize, usize, Vec<[u8; 3]>) {
    let w: usize = images.iter().map(|(w, _)| w).sum::<usize>() + GAP * (images.len() - 1);
    let mut out = vec![[255; 3]; w * h];
    let mut x0 = 0;
    for (iw, px) in images {
        for y in 0..h {
            out[y * w + x0..][..*iw].copy_from_slice(&px[y * iw..][..*iw]);
        }
        x0 += iw + GAP;
    }
    (w, h, out)
}

pub fn calibrate_sheet(input: &Path, opts: &Options, seams: Option<&Seams>) -> Result<SheetReport> {
    let src = Input::open(input)?;
    let a = analyse(&src, opts, seams).with_context(|| format!("while analysing {}", input.display()))?;
    let t = &a.thumb;
    let tiles = [
        t.pixels.iter().map(|p| srgb8(p.map(|c| c as f64))).collect::<Vec<_>>(),
        t.pixels.iter().zip(&a.surface.role).map(|(p, r)| role_colour(*r, *p)).collect(),
        result_image(&a, t),
    ];
    let tiles: Vec<_> = tiles.iter().map(|px| resize(t.width, t.height, px, ROW_HEIGHT)).collect();
    Ok(SheetReport {
        name: input.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default(),
        paper: a.paper.global,
        degree: a.surface.degree,
        paper_fraction: a.paper_fraction(),
        seam_edges: a.seams.as_ref().map_or(0, |s| s.edges()),
        row: hstack(&tiles, ROW_HEIGHT),
    })
}

/// Settings echoed in the proposed `homog.toml`.
pub struct Proposal {
    pub target: [u8; 3],
    pub icc: bool,
    pub gcp_dir: Option<PathBuf>,
    pub paper_hints: Option<PathBuf>,
    pub keep_hints: Option<PathBuf>,
}

fn absolute(p: &Path) -> String {
    p.canonicalize().unwrap_or_else(|_| p.to_owned()).display().to_string()
}

pub fn write(dir: &Path, sheets: &[SheetReport], proposal: &Proposal) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("cannot create {}", dir.display()))?;

    // Contact sheet
    let width = sheets.iter().map(|s| s.row.0).max().unwrap_or(1);
    let height = sheets.len() * (ROW_HEIGHT + GAP);
    let mut contact = vec![[255; 3]; width * height.max(1)];
    for (k, s) in sheets.iter().enumerate() {
        let (w, h, px) = &s.row;
        for y in 0..*h {
            contact[(k * (ROW_HEIGHT + GAP) + y) * width..][..*w].copy_from_slice(&px[y * w..][..*w]);
        }
    }
    write_png(&dir.join("contact_sheet.png"), width, height.max(1), &contact)?;

    // Report
    let batch = [0, 1, 2].map(|c| median(sheets.iter().map(|s| s.paper[c]).collect()));
    let mut txt = String::new();
    writeln!(txt, "{} sheet(s); batch paper colour (CIELAB): ({:.1}, {:.1}, {:.1})", sheets.len(), batch[0], batch[1], batch[2])?;
    writeln!(txt, "Rows of contact_sheet.png in this order: original | pixels used as paper | result.\n")?;
    writeln!(txt, "{:>4}  {:<44} {:>22} {:>8} {:>7} {:>7} {:>6}  warnings", "row", "sheet", "paper Lab", "ΔE batch", "paper", "degree", "seams")?;
    let mut flagged = 0;
    for (k, s) in sheets.iter().enumerate() {
        let de = ((s.paper[0] - batch[0]).powi(2) + (s.paper[1] - batch[1]).powi(2) + (s.paper[2] - batch[2]).powi(2)).sqrt();
        let mut warnings = Vec::new();
        if sheets.len() > 2 && de > 6.0 {
            warnings.push("paper far from the batch");
        }
        if s.paper_fraction < 0.15 {
            warnings.push("little bare paper found");
        }
        flagged += usize::from(!warnings.is_empty());
        writeln!(
            txt,
            "{:>4}  {:<44} {:>22} {:>8.1} {:>6.0}% {:>7} {:>6}  {}",
            k + 1,
            s.name,
            format!("({:.1}, {:.1}, {:.1})", s.paper[0], s.paper[1], s.paper[2]),
            de,
            s.paper_fraction * 100.0,
            s.degree.map_or("-".into(), |d| d.to_string()),
            s.seam_edges,
            warnings.join("; ")
        )?;
    }
    writeln!(txt)?;
    if flagged > 0 {
        writeln!(
            txt,
            "{flagged} sheet(s) with warnings: check their row of the contact sheet. If something other\n\
             than paper is green (e.g. a pale wash), crop a sample of it into the `keep` hint directory;\n\
             if the paper is not found, crop a sample of bare paper into the `paper` hint directory."
        )?;
    } else {
        writeln!(txt, "No warning: the automatic estimate looks consistent across the batch.")?;
    }
    std::fs::write(dir.join("calibrate.txt"), txt)?;

    // Proposed configuration
    let mut toml = String::new();
    writeln!(toml, "# Proposed by `homog calibrate` on {} sheet(s). See calibrate.txt and contact_sheet.png.", sheets.len())?;
    writeln!(toml, "# Batch paper colour (CIELAB): ({:.1}, {:.1}, {:.1})", batch[0], batch[1], batch[2])?;
    writeln!(toml, "target = [{}, {}, {}]", proposal.target[0], proposal.target[1], proposal.target[2])?;
    let mut hist = std::collections::BTreeMap::new();
    sheets.iter().filter_map(|s| s.degree).for_each(|d| *hist.entry(d).or_insert(0) += 1);
    let hist: Vec<String> = hist.iter().map(|(d, n)| format!("degree {d}: {n}")).collect();
    writeln!(toml, "# Paper surface degrees chosen: {}", if hist.is_empty() { "none".into() } else { hist.join(", ") })?;
    writeln!(toml, "lighting = \"auto\"   # auto | even | uneven | none")?;
    writeln!(toml, "icc = {}", proposal.icc)?;
    writeln!(toml, "\n[hints]")?;
    writeln!(toml, "# Imagettes: crops of the scans (at least 10 x 10 px), one example per image file.")?;
    match &proposal.paper_hints {
        Some(p) => writeln!(toml, "paper = \"{}\"", absolute(p))?,
        None => writeln!(toml, "# paper = \"hints/paper\"   # bare paper, when the automatic estimate picks something else")?,
    }
    match &proposal.keep_hints {
        Some(p) => writeln!(toml, "keep = \"{}\"", absolute(p))?,
        None => writeln!(toml, "# keep = \"hints/keep\"     # colours that must not become white (pale washes close to the paper)")?,
    }
    writeln!(toml, "\n[seams]")?;
    match &proposal.gcp_dir {
        Some(p) => writeln!(toml, "gcp_dir = \"{}\"", absolute(p))?,
        None => writeln!(toml, "# gcp_dir = \"GCP\"   # QGIS .points files, for maps made of several sheets")?,
    }
    std::fs::write(dir.join("homog.toml"), toml)?;
    Ok(())
}
