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

/// Robust colour model in CIELAB: centre, a*b* covariance and L* spread.
#[derive(Clone, Debug)]
pub struct ColourModel {
    pub mean: Vec3,
    inv_ab: [[f64; 2]; 2],
    /// Largest a*b* standard deviation (after flooring)
    pub sd_ab: f64,
    pub sd_l: f64,
}

impl ColourModel {
    /// Model of a set of Lab pixels: median centre; covariance of the pixels within 3 robust
    /// standard deviations of it; a*b* standard deviations floored at `min_sd_ab`.
    pub fn from_pixels(px: &[Vec3], min_sd_ab: f64) -> Option<Self> {
        if px.len() < 10 {
            return None;
        }
        let mean = [0, 1, 2].map(|c| median(px.iter().map(|p| p[c]).collect()));
        let dist = |p: &Vec3| ((p[1] - mean[1]).powi(2) + (p[2] - mean[2]).powi(2)).sqrt();
        let mad = 1.4826 * median(px.iter().map(dist).collect());
        let core: Vec<&Vec3> = px.iter().filter(|p| dist(p) <= 3.0 * mad.max(0.5)).collect();
        let n = core.len() as f64;
        let cov = |i: usize, j: usize| core.iter().map(|p| (p[i] - mean[i]) * (p[j] - mean[j])).sum::<f64>() / n;
        let sd_l = 1.4826 * median(px.iter().map(|p| (p[0] - mean[0]).abs()).collect());
        // Eigen-decomposition of the 2x2 covariance, eigenvalues floored
        let (saa, sab, sbb) = (cov(1, 1), cov(1, 2), cov(2, 2));
        let tr = saa + sbb;
        let det = saa * sbb - sab * sab;
        let disc = ((tr * tr / 4.0) - det).max(0.0).sqrt();
        let (l1, l2) = (tr / 2.0 + disc, tr / 2.0 - disc);
        let (v1, v2) = if sab.abs() > 1e-12 {
            let v = [l1 - sbb, sab];
            let n = v[0].hypot(v[1]);
            ([v[0] / n, v[1] / n], [-v[1] / n, v[0] / n])
        } else if saa >= sbb {
            ([1.0, 0.0], [0.0, 1.0])
        } else {
            ([0.0, 1.0], [1.0, 0.0])
        };
        let floor = min_sd_ab * min_sd_ab;
        let (e1, e2) = (l1.max(floor), l2.max(floor));
        let inv = |i: usize, j: usize| v1[i] * v1[j] / e1 + v2[i] * v2[j] / e2;
        Some(Self { mean, inv_ab: [[inv(0, 0), inv(0, 1)], [inv(1, 0), inv(1, 1)]], sd_ab: e1.sqrt(), sd_l })
    }

    /// Mahalanobis distance in a*b*.
    pub fn distance_ab(&self, lab: &Vec3) -> f64 {
        let (da, db) = (lab[1] - self.mean[1], lab[2] - self.mean[2]);
        let m = &self.inv_ab;
        (da * (m[0][0] * da + m[0][1] * db) + db * (m[1][0] * da + m[1][1] * db)).max(0.0).sqrt()
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
    fn model_of_an_elongated_cloud() {
        // a* spread 4 along a*, b* spread 1: Mahalanobis distance scales accordingly
        let px: Vec<Vec3> = (0..2000)
            .map(|i| {
                let t = (i as f64 * 0.618).fract() * 2.0 - 1.0;
                let u = (i as f64 * 0.414).fract() * 2.0 - 1.0;
                [70.0, 2.0 + 4.0 * 3f64.sqrt() * t, 15.0 + 3f64.sqrt() * u]
            })
            .collect();
        let m = ColourModel::from_pixels(&px, 0.5).unwrap();
        assert!((m.mean[1] - 2.0).abs() < 0.3 && (m.mean[2] - 15.0).abs() < 0.3, "{:?}", m.mean);
        assert!((m.sd_ab - 4.0).abs() < 0.4, "sd {}", m.sd_ab);
        let along_a = m.distance_ab(&[70.0, 2.0 + 8.0, 15.0]);
        let along_b = m.distance_ab(&[70.0, 2.0, 15.0 + 2.0]);
        assert!((along_a - 2.0).abs() < 0.3 && (along_b - 2.0).abs() < 0.3, "{along_a} {along_b}");
    }

    #[test]
    fn spread_is_floored() {
        let px = vec![[70.0, 1.0, 15.0]; 100];
        let m = ColourModel::from_pixels(&px, 2.0).unwrap();
        assert!((m.sd_ab - 2.0).abs() < 1e-9);
        assert!((m.distance_ab(&[70.0, 1.0, 19.0]) - 2.0).abs() < 1e-9);
    }
}
