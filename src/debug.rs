//! Debug images: what the analysis estimated and what the correction does, per sheet.
//!
//! Written to `<dir>/`:
//! - `1_original.png`: input (reduced);
//! - `2_paper_surface.png`: estimated paper colour (flat-field);
//! - `3_paper_pixels.png`: role of each pixel in the paper fit (see [`role_colour`]);
//! - `4_gain_surface.png`: luminance gain of the surface correction (scale in `debug.txt`);
//! - `5_gain_seams.png`: luminance gain of the seam correction, blue = darker, white = 1,
//!   red = lighter (symmetric scale in `debug.txt`);
//! - `6_seam_samples.png`: result with the seam measurement points: green = paper, red = rejected;
//! - `7_result.png`: result (reduced);
//! - `debug.txt`: numbers behind the images.

use std::fmt::Write as _;
use std::path::Path;

use anyhow::{Context, Result};
use rayon::prelude::*;

use crate::background::PixelRole;
use crate::color::{lab_to_srgb, linear_rgb_to_lab, linear_to_srgb};
use crate::io::write_png;
use crate::seams::luminance_gain;
use crate::transform::GainMap;
use crate::{Analysis, Reduced};

pub fn srgb8(lin: [f64; 3]) -> [u8; 3] {
    lin.map(|c| (linear_to_srgb(c) * 255.0).round() as u8)
}

fn lerp_colors(stops: &[[f64; 3]], t: f64) -> [u8; 3] {
    let x = t.clamp(0.0, 1.0) * (stops.len() - 1) as f64;
    let i = (x as usize).min(stops.len() - 2);
    let f = x - i as f64;
    [0, 1, 2].map(|c| (stops[i][c] + (stops[i + 1][c] - stops[i][c]) * f).round() as u8)
}

/// Sequential colour map (magma-like) on [0, 1].
fn sequential(t: f64) -> [u8; 3] {
    lerp_colors(&[[0.0, 0.0, 4.0], [81.0, 18.0, 124.0], [183.0, 55.0, 121.0], [252.0, 137.0, 97.0], [252.0, 253.0, 191.0]], t)
}

/// Diverging colour map on [-1, 1]: blue, white, red.
fn diverging(t: f64) -> [u8; 3] {
    lerp_colors(&[[33.0, 102.0, 172.0], [247.0, 247.0, 247.0], [178.0, 24.0, 43.0]], (t + 1.0) / 2.0)
}

/// Luminance gain of a map of LMS gains.
fn luminance(g: &GainMap) -> Vec<f64> {
    g.data.iter().map(|v| luminance_gain(&v.map(|c| c as f64))).collect()
}

/// Colour of a thumbnail pixel by role in the paper fit, shaded by its lightness:
/// green = bare paper (L* and a*b*), orange = L* envelope only, grey = ignored (ink, strong
/// washes, `keep` colours), dark red = outside the sheet (scanner bed, colour charts).
pub fn role_colour(role: PixelRole, linear: [f32; 3]) -> [u8; 3] {
    let l = linear_rgb_to_lab(&linear.map(|c| c as f64))[0] / 100.0;
    let (tint, k) = match role {
        PixelRole::Chroma => ([40.0, 170.0, 60.0], 0.35 + 0.65 * l),
        PixelRole::Lightness => ([240.0, 150.0, 30.0], 0.35 + 0.65 * l),
        PixelRole::Other => ([255.0, 255.0, 255.0], 0.6 * l),
        PixelRole::Background => ([150.0, 20.0, 40.0], 0.3 + 0.4 * l),
    };
    tint.map(|c: f64| (c * k) as u8)
}

/// Rendering of a reduced image with the final correction.
pub fn result_image(a: &Analysis, img: &Reduced) -> Vec<[u8; 3]> {
    img.pixels
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let s = img.factor as f64;
            let px = [(((i % img.width) as f64 + 0.5) * s).min(a.width as f64), (((i / img.width) as f64 + 0.5) * s).min(a.height as f64)];
            srgb8(a.adapter.adapt_point(px, p.map(|c| c as f64)))
        })
        .collect()
}

pub fn write(dir: &Path, a: &Analysis) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("cannot create {}", dir.display()))?;
    let t = &a.thumb;
    let m = a.medium.as_ref().context("debug images need the reduced image")?;
    let mut txt = String::new();
    let g = a.paper.global;
    writeln!(txt, "global paper Lab: ({:.2}, {:.2}, {:.2}), a*b* spread {:.2}", g[0], g[1], g[2], a.paper.sd_ab())?;
    match a.surface.degree {
        Some(d) => writeln!(txt, "paper surface degree: {d}")?,
        None => writeln!(txt, "single paper colour (no flat-field)")?,
    }

    // 1. Original
    write_png(&dir.join("1_original.png"), m.width, m.height, &m.pixels.iter().map(|p| srgb8(p.map(|c| c as f64))).collect::<Vec<_>>())?;

    // 2. Paper surface
    let surf: Vec<[u8; 3]> = a.surface.lab.iter().map(|lab| lab_to_srgb(lab).map(|c| (c.clamp(0.0, 1.0) * 255.0).round() as u8)).collect();
    write_png(&dir.join("2_paper_surface.png"), t.width, t.height, &surf)?;
    let (lo, hi) = a.paper_l_range();
    writeln!(txt, "paper surface L*: {lo:.1} .. {hi:.1}")?;

    // 3. Pixels used by the fit
    let roles: Vec<[u8; 3]> = t.pixels.iter().zip(&a.surface.role).map(|(p, r)| role_colour(*r, *p)).collect();
    write_png(&dir.join("3_paper_pixels.png"), t.width, t.height, &roles)?;
    let count = |r: PixelRole| a.surface.role.iter().filter(|x| **x == r).count() as f64 / a.surface.role.len() as f64 * 100.0;
    writeln!(
        txt,
        "pixels: bare paper {:.1} %, L* envelope only {:.1} %, ignored {:.1} %, outside the sheet {:.1} %",
        count(PixelRole::Chroma),
        count(PixelRole::Lightness),
        count(PixelRole::Other),
        count(PixelRole::Background)
    )?;

    // 4. Surface gain
    let gs = luminance(&a.surface_gains);
    let (gmin, gmax) = gs.iter().fold((f64::MAX, f64::MIN), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
    let span = (gmax - 1.0).max(1e-6);
    let img: Vec<[u8; 3]> = gs.iter().map(|v| sequential((v - 1.0) / span)).collect();
    write_png(&dir.join("4_gain_surface.png"), a.surface_gains.width, a.surface_gains.height, &img)?;
    writeln!(txt, "surface luminance gain: {gmin:.3} .. {gmax:.3} (image: black = 1.0, light yellow = {gmax:.3})")?;

    // 5. Seam gain (ratio of the final gains to the surface gains)
    let total = luminance(a.adapter.gains());
    let ratio: Vec<f64> = total.iter().zip(&gs).map(|(x, y)| x / y).collect();
    let (rmin, rmax) = ratio.iter().fold((f64::MAX, f64::MIN), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
    // Symmetric scale covering the largest deviation, at least ±5 %
    let scale = (rmax - 1.0).max(1.0 - rmin).max(0.05);
    let img: Vec<[u8; 3]> = ratio.iter().map(|r| diverging((r - 1.0) / scale)).collect();
    write_png(&dir.join("5_gain_seams.png"), a.surface_gains.width, a.surface_gains.height, &img)?;
    writeln!(txt, "seam luminance gain: {rmin:.3} .. {rmax:.3} (image: blue = {:.3}, white = 1, red = {:.3})", 1.0 - scale, 1.0 + scale)?;

    // 7. Result (reduced), 6. with the seam samples
    let result = result_image(a, m);
    write_png(&dir.join("7_result.png"), m.width, m.height, &result)?;
    let mut with_samples = result;
    if let Some(field) = &a.seams {
        for &(px, paper) in &field.samples {
            let (x, y) = ((px[0] / m.factor as f64) as isize, (px[1] / m.factor as f64) as isize);
            if (0..m.width as isize).contains(&x) && (0..m.height as isize).contains(&y) {
                with_samples[y as usize * m.width + x as usize] = if paper { [0, 200, 0] } else { [230, 0, 0] };
            }
        }
        let n_paper = field.samples.iter().filter(|s| s.1).count();
        writeln!(txt, "seam band: {:.2} .. {:.2} map units from the edges", field.band.0, field.band.1)?;
        writeln!(txt, "seam samples: {} ({} taken as paper)", field.samples.len(), n_paper)?;
        for (name, gains) in field.describe() {
            let list: Vec<String> = gains.iter().map(|g| format!("{:+.1}", (g - 1.0) * 100.0)).collect();
            writeln!(txt, "  {name}: luminance gain per segment (%): {}", list.join(" "))?;
        }
    } else {
        writeln!(txt, "no seam correction")?;
    }
    write_png(&dir.join("6_seam_samples.png"), m.width, m.height, &with_samples)?;

    std::fs::write(dir.join("debug.txt"), txt)?;
    Ok(())
}
