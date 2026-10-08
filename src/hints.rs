//! Colour hints given as small images cropped from the scans ("imagettes"):
//! `paper/` shows what bare paper looks like, `keep/` shows colours that must never be taken
//! for paper (e.g. a pale yellow wash close to yellowed paper). Each imagette gives a colour
//! model: robust centre and spread, from which the tolerances are derived.

use std::path::Path;

use anyhow::{bail, ensure, Context, Result};

use crate::background::median;
use crate::color::{linear_rgb_to_lab, srgb_to_linear, Vec3};
use crate::io::Input;

/// Mahalanobis radius (in a*b*) within which a colour matches a model
pub const MATCH_RADIUS: f64 = 3.5;

/// Robust colour model in CIELAB: centre, a*b* and L* spreads.
#[derive(Clone, Debug)]
pub struct ColourModel {
    pub mean: Vec3,
    /// a*b* standard deviation (per axis, after flooring)
    pub sd_ab: f64,
    pub sd_l: f64,
}

impl ColourModel {
    /// Model of a set of Lab pixels: median centre, robust spreads (the median a*b* distance
    /// to the centre is 1.18 standard deviations for a 2-D normal spread), a*b* standard
    /// deviation floored at `min_sd_ab`.
    pub fn from_pixels(px: &[Vec3], min_sd_ab: f64) -> Option<Self> {
        if px.len() < 10 {
            return None;
        }
        let mean = [0, 1, 2].map(|c| median(px.iter().map(|p| p[c]).collect()));
        let sd_ab = median(px.iter().map(|p| (p[1] - mean[1]).hypot(p[2] - mean[2])).collect()) / 1.1774;
        let sd_l = 1.4826 * median(px.iter().map(|p| (p[0] - mean[0]).abs()).collect());
        Some(Self { mean, sd_ab: sd_ab.max(min_sd_ab), sd_l })
    }

    /// a*b* distance to the centre, in standard deviations.
    pub fn distance_ab(&self, lab: &Vec3) -> f64 {
        (lab[1] - self.mean[1]).hypot(lab[2] - self.mean[2]) / self.sd_ab
    }
}

/// Colour hints of a project.
#[derive(Clone, Debug, Default)]
pub struct Hints {
    pub paper: Vec<ColourModel>,
    pub keep: Vec<ColourModel>,
}

/// Floor of the a*b* standard deviation of an imagette model
const SAMPLE_MIN_SD_AB: f64 = 1.5;
/// Floor of the a*b* standard deviation of a `keep` model
const KEEP_MIN_SD_AB: f64 = 1.0;

impl Hints {
    pub fn load(paper: Option<&Path>, keep: Option<&Path>) -> Result<Self> {
        Ok(Self {
            paper: paper.map(|d| load_dir(d, SAMPLE_MIN_SD_AB)).transpose()?.unwrap_or_default(),
            keep: keep.map(|d| load_dir(d, KEEP_MIN_SD_AB)).transpose()?.unwrap_or_default(),
        })
    }

    /// Whether a colour (as scanned) matches one of the `keep` imagettes.
    pub fn keeps(&self, lab: &Vec3) -> bool {
        self.keep
            .iter()
            .any(|m| m.distance_ab(lab) < 3.0 && (lab[0] - m.mean[0]).abs() < (3.0 * m.sd_l).max(15.0))
    }
}

/// One colour model per image of a directory (PNG, JPEG or TIFF, read as sRGB).
fn load_dir(dir: &Path, min_sd_ab: f64) -> Result<Vec<ColourModel>> {
    let mut paths: Vec<_> = std::fs::read_dir(dir)
        .with_context(|| format!("cannot read the hint directory {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .is_some_and(|e| ["png", "jpg", "jpeg", "tif", "tiff"].contains(&e.to_lowercase().as_str()))
        })
        .collect();
    paths.sort();
    if paths.is_empty() {
        bail!("no image in the hint directory {}", dir.display());
    }
    paths.iter().map(|p| load_imagette(p, min_sd_ab).with_context(|| format!("hint {}", p.display()))).collect()
}

fn load_imagette(path: &Path, min_sd_ab: f64) -> Result<ColourModel> {
    let input = Input::open(path)?;
    let (_, _, _, px) = input.read_thumbnail(256)?;
    let lab: Vec<Vec3> = px.iter().map(|p| linear_rgb_to_lab(&p.map(|c| srgb_to_linear(c as f64)))).collect();
    ensure!(lab.len() >= 10, "imagette too small ({} pixels), crop at least 10 x 10 pixels", lab.len());
    ColourModel::from_pixels(&lab, min_sd_ab).context("cannot model the imagette colours")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn model_of_a_cloud() {
        // a*b* uniformly spread on a disc of radius 6 around (2, 15), plus 10 % outliers
        let mut px: Vec<Vec3> = (0..2000)
            .map(|i| {
                let r = 6.0 * ((i as f64 * 0.618_034).fract()).sqrt();
                let t = std::f64::consts::TAU * (i as f64 * 0.414_214).fract();
                [70.0, 2.0 + r * t.cos(), 15.0 + r * t.sin()]
            })
            .collect();
        px.extend((0..200).map(|_| [30.0, 40.0, -20.0]));
        let m = ColourModel::from_pixels(&px, 0.5).unwrap();
        assert!((m.mean[1] - 2.0).abs() < 0.5 && (m.mean[2] - 15.0).abs() < 0.5, "{:?}", m.mean);
        // Median radius of the disc: 6 / sqrt(2) = 4.24 -> sd 3.6
        assert!((m.sd_ab - 3.6).abs() < 0.4, "sd {}", m.sd_ab);
    }

    #[test]
    fn spread_is_floored() {
        let px = vec![[70.0, 1.0, 15.0]; 100];
        let m = ColourModel::from_pixels(&px, 2.0).unwrap();
        assert!((m.sd_ab - 2.0).abs() < 1e-9);
        assert!((m.distance_ab(&[70.0, 1.0, 19.0]) - 2.0).abs() < 1e-9);
    }
}
