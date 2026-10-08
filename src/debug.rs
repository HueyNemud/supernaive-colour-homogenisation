//! Debug images: what the analysis estimated and what the correction does, per sheet.
//!
//! Written to `<dir>/`:
//! - `1_original.png`: input (reduced);
//! - `2_paper_surface.png`: estimated paper colour (flat-field);
//! - `3_paper_pixels.png`: pixels used for the surface fit: green = bare paper (L* and a*b*),
//!   orange = L* envelope only, grey = ignored (ink, strong washes, background);
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

use crate::background::{PaperSurface, PixelRole};
use crate::color::{lab_to_srgb, linear_rgb_to_lab, linear_to_srgb, Vec3};
use crate::io::write_png;
use crate::seams::{luminance_gain, EdgeField};
use crate::transform::{Adapter, GainMap};

/// A reduced linear RGB image and its reduction factor w.r.t. the full resolution.
pub struct Reduced<'a> {
    pub width: usize,
    pub height: usize,
    pub factor: usize,
    pub pixels: &'a [[f32; 3]],
}

pub struct DebugInput<'a> {
    pub full_size: (usize, usize),
    pub global_paper: Vec3,
    pub thumb: Reduced<'a>,
    pub surface: &'a PaperSurface,
    pub surface_gains: &'a GainMap,
    pub seams: Option<&'a EdgeField>,
    pub medium: Reduced<'a>,
    pub adapter: &'a Adapter,
}

fn srgb8(lin: [f64; 3]) -> [u8; 3] {
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

pub fn write(dir: &Path, d: &DebugInput) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("cannot create {}", dir.display()))?;
    let (t, m) = (&d.thumb, &d.medium);
    let mut txt = String::new();
    let g = d.global_paper;
    writeln!(txt, "global paper Lab: ({:.2}, {:.2}, {:.2})", g[0], g[1], g[2])?;

    // 1. Original
    write_png(&dir.join("1_original.png"), m.width, m.height, &m.pixels.iter().map(|p| srgb8(p.map(|c| c as f64))).collect::<Vec<_>>())?;

    // 2. Paper surface
    let surf: Vec<[u8; 3]> = d.surface.lab.iter().map(|lab| lab_to_srgb(lab).map(|c| (c.clamp(0.0, 1.0) * 255.0).round() as u8)).collect();
    write_png(&dir.join("2_paper_surface.png"), t.width, t.height, &surf)?;
    let (lo, hi) = d.surface.lab.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])));
    writeln!(txt, "paper surface L*: {lo:.1} .. {hi:.1}")?;

    // 3. Pixels used by the fit
    let roles: Vec<[u8; 3]> = t
        .pixels
        .iter()
        .zip(&d.surface.role)
        .map(|(p, role)| {
            let l = linear_rgb_to_lab(&p.map(|c| c as f64))[0] / 100.0;
            let tint = match role {
                PixelRole::Chroma => [40.0, 170.0, 60.0],
                PixelRole::Lightness => [240.0, 150.0, 30.0],
                PixelRole::Other => [255.0, 255.0, 255.0],
            };
            let k = if *role == PixelRole::Other { 0.6 * l } else { 0.35 + 0.65 * l };
            tint.map(|c| (c * k) as u8)
        })
        .collect();
    write_png(&dir.join("3_paper_pixels.png"), t.width, t.height, &roles)?;
    let count = |r: PixelRole| d.surface.role.iter().filter(|x| **x == r).count() as f64 / d.surface.role.len() as f64 * 100.0;
    writeln!(
        txt,
        "fit pixels: bare paper {:.1} %, L* envelope only {:.1} %, ignored {:.1} %",
        count(PixelRole::Chroma),
        count(PixelRole::Lightness),
        count(PixelRole::Other)
    )?;

    // 4. Surface gain
    let gs = luminance(d.surface_gains);
    let (gmin, gmax) = gs.iter().fold((f64::MAX, f64::MIN), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
    let span = (gmax - 1.0).max(1e-6);
    let img: Vec<[u8; 3]> = gs.iter().map(|v| sequential((v - 1.0) / span)).collect();
    write_png(&dir.join("4_gain_surface.png"), d.surface_gains.width, d.surface_gains.height, &img)?;
    writeln!(txt, "surface luminance gain: {gmin:.3} .. {gmax:.3} (image: black = 1.0, light yellow = {gmax:.3})")?;

    // 5. Seam gain (ratio of the final gains to the surface gains)
    let total = luminance(d.adapter.gains());
    let ratio: Vec<f64> = total.iter().zip(&gs).map(|(a, b)| a / b).collect();
    let (rmin, rmax) = ratio.iter().fold((f64::MAX, f64::MIN), |(lo, hi), v| (lo.min(*v), hi.max(*v)));
    // Symmetric scale covering the largest deviation, at least ±5 %
    let scale = (rmax - 1.0).max(1.0 - rmin).max(0.05);
    let img: Vec<[u8; 3]> = ratio.iter().map(|r| diverging((r - 1.0) / scale)).collect();
    write_png(&dir.join("5_gain_seams.png"), d.surface_gains.width, d.surface_gains.height, &img)?;
    writeln!(txt, "seam luminance gain: {rmin:.3} .. {rmax:.3} (image: blue = {:.3}, white = 1, red = {:.3})", 1.0 - scale, 1.0 + scale)?;

    // 7. Result (reduced), 6. with the seam samples
    let (w, h) = d.full_size;
    let result: Vec<[u8; 3]> = m
        .pixels
        .par_iter()
        .enumerate()
        .map(|(i, p)| {
            let px = [((i % m.width) as f64 + 0.5) * m.factor as f64, ((i / m.width) as f64 + 0.5) * m.factor as f64];
            srgb8(d.adapter.adapt_point([px[0].min(w as f64), px[1].min(h as f64)], p.map(|c| c as f64)))
        })
        .collect();
    write_png(&dir.join("7_result.png"), m.width, m.height, &result)?;
    let mut with_samples = result;
    if let Some(field) = d.seams {
        for &(px, paper) in &field.samples {
            let (x, y) = ((px[0] / m.factor as f64) as isize, (px[1] / m.factor as f64) as isize);
            if (0..m.width as isize).contains(&x) && (0..m.height as isize).contains(&y) {
                with_samples[y as usize * m.width + x as usize] = if paper { [0, 200, 0] } else { [230, 0, 0] };
            }
        }
        let n_paper = field.samples.iter().filter(|s| s.1).count();
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
