//! Seam correction between adjacent sheets. Mirrors `edge_gain_field` of the prototype.
//!
//! The paper surface cannot follow the darkening close to the edges of a sheet, which shows
//! as a step where two sheets meet. Along every edge shared with another sheet, the paper is
//! measured (after the surface correction) in a thin band inside the neatline, per segment,
//! and brought back to the target white by a Bradford-LMS gain that fades out inwards.

use crate::background::median;
use crate::color::{linear_rgb_to_lab, mat_vec, Vec3, M_BRADFORD, M_RGB2XYZ};
use crate::layout::{Edge, Layout};

/// Band where the paper is measured, from the edge (fractions of the sheet size)
const BAND: (f64, f64) = (0.003, 0.013);
/// Distance over which the correction fades out (fraction of the sheet size)
const REACH: f64 = 0.12;
const BINS: usize = 24;
const ALONG: usize = 400;
const ACROSS: usize = 8;

struct EdgeProfile {
    edge: Edge,
    /// Log LMS gain per segment along the edge
    gains: Vec<[f64; 3]>,
}

/// A measurement point: full-resolution pixel position and whether it was taken as paper.
pub type Sample = ([f64; 2], bool);

/// Log-gain field of one sheet.
pub struct EdgeField {
    layout: Layout,
    profiles: Vec<EdgeProfile>,
    reach: f64,
    /// All measurement points (for the debug images)
    pub samples: Vec<Sample>,
}

/// Factor applied to the luminance Y of a neutral colour by the Bradford-LMS gains `g`.
pub fn luminance_gain(g: &[f64; 3]) -> f64 {
    let w = crate::color::white_d65();
    let lms = mat_vec(&M_BRADFORD, &w);
    let scaled = [lms[0] * g[0], lms[1] * g[1], lms[2] * g[2]];
    mat_vec(&crate::color::mat_inv(&M_BRADFORD), &scaled)[1] / w[1]
}

/// np.interp on regularly spaced bin centres, clamped at both ends.
fn interp_bins(u: f64, values: &[[f64; 3]]) -> [f64; 3] {
    let n = values.len();
    let x = (u.clamp(0.0, 1.0) * n as f64 - 0.5).clamp(0.0, (n - 1) as f64);
    let i = (x.floor() as usize).min(n - 1);
    let j = (i + 1).min(n - 1);
    let f = x - i as f64;
    [0, 1, 2].map(|c| values[i][c] + (values[j][c] - values[i][c]) * f)
}

impl EdgeField {
    /// Measure the paper along `edges`. `sample` returns the linear RGB colour (after the
    /// surface correction) at a full-resolution pixel position, if inside the image.
    pub fn measure(layout: &Layout, edges: &[Edge], target_xyz: &Vec3, sample: impl Fn([f64; 2]) -> Option<[f64; 3]>) -> Self {
        let size = layout.size();
        let lms_t = mat_vec(&M_BRADFORD, target_xyz);
        let to_lms = |rgb: &[f64; 3]| mat_vec(&M_BRADFORD, &mat_vec(&M_RGB2XYZ, rgb));
        let mut samples = Vec::new();
        let profiles = edges
            .iter()
            .filter_map(|&edge| {
                let (u0, u1) = edge.span;
                let mut bins: Vec<Vec<[f64; 3]>> = vec![Vec::new(); BINS];
                for a in 0..ALONG {
                    let f = a as f64 / (ALONG - 1) as f64;
                    let u = u0 + (u1 - u0) * f;
                    for k in 0..ACROSS {
                        let t = size * (BAND.0 + (BAND.1 - BAND.0) * k as f64 / (ACROSS - 1) as f64);
                        let d = edge.pos + edge.inward * t;
                        let xy = if edge.axis == 0 { [d, u] } else { [u, d] };
                        let px = layout.pixel(xy);
                        let Some(rgb) = sample(px) else { continue };
                        let lab = linear_rgb_to_lab(&rgb);
                        let is_paper = lab[1].hypot(lab[2]) < 4.0 && lab[0] > 80.0;
                        if is_paper {
                            bins[((f * BINS as f64) as usize).min(BINS - 1)].push(to_lms(&rgb));
                        }
                        samples.push((px, is_paper));
                    }
                }
                // Gain per segment with enough paper, interpolated over the others
                let known: Vec<(f64, [f64; 3])> = bins
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| b.len() >= 10)
                    .map(|(k, b)| {
                        let g = [0, 1, 2].map(|c| (lms_t[c] / median(b.iter().map(|p| p[c]).collect())).ln());
                        ((k as f64 + 0.5) / BINS as f64, g)
                    })
                    .collect();
                if known.is_empty() {
                    return None;
                }
                let filled: Vec<[f64; 3]> = (0..BINS)
                    .map(|k| {
                        let c = (k as f64 + 0.5) / BINS as f64;
                        match known.iter().position(|(x, _)| *x >= c) {
                            None => known.last().unwrap().1,
                            Some(0) => known[0].1,
                            Some(i) => {
                                let ((x0, g0), (x1, g1)) = (known[i - 1], known[i]);
                                let f = (c - x0) / (x1 - x0);
                                [0, 1, 2].map(|ch| g0[ch] + (g1[ch] - g0[ch]) * f)
                            }
                        }
                    })
                    .collect();
                // 3-segment moving average (edges repeated)
                let gains = (0..BINS)
                    .map(|k| {
                        let idx = [k.saturating_sub(1), k, (k + 1).min(BINS - 1)];
                        [0, 1, 2].map(|c| idx.iter().map(|&i| filled[i][c]).sum::<f64>() / 3.0)
                    })
                    .collect();
                Some(EdgeProfile { edge, gains })
            })
            .collect();
        Self { layout: layout.clone(), profiles, reach: REACH * size, samples }
    }

    pub fn edges(&self) -> usize {
        self.profiles.len()
    }

    /// Per corrected edge: description and luminance gain (Y factor) per segment.
    pub fn describe(&self) -> Vec<(String, Vec<f64>)> {
        self.profiles
            .iter()
            .map(|p| {
                let e = &p.edge;
                let name = format!(
                    "{} edge {} = {:.1} (span {:.1}..{:.1})",
                    if e.axis == 0 { "vertical" } else { "horizontal" },
                    if e.axis == 0 { "x" } else { "y" },
                    e.pos,
                    e.span.0,
                    e.span.1
                );
                (name, p.gains.iter().map(|g| luminance_gain(&g.map(f64::exp))).collect())
            })
            .collect()
    }

    /// Log LMS gain at a full-resolution pixel position. Near a corner, the corrections of
    /// the two edges measure the same darkening: they are averaged rather than added
    /// (weights normalised when they sum to more than 1).
    pub fn log_gain(&self, px: [f64; 2]) -> [f64; 3] {
        let xy = self.layout.map(px);
        let mut g = [0.0; 3];
        let mut weight = 0.0;
        for p in &self.profiles {
            let (across, along) = if p.edge.axis == 0 { (xy[0], xy[1]) } else { (xy[1], xy[0]) };
            let d = p.edge.inward * (across - p.edge.pos);
            let phi = (1.0 - d / self.reach).clamp(0.0, 1.0).powi(2);
            if phi > 0.0 {
                let (u0, u1) = p.edge.span;
                let v = interp_bins((along - u0) / (u1 - u0), &p.gains);
                (0..3).for_each(|c| g[c] += phi * v[c]);
                weight += phi;
            }
        }
        g.map(|v| v / weight.max(1.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::color::{srgb_to_linear, white_d65};

    /// Sheet 600 x 400 map units, 1 px per unit, y up; paper darkening over the last 30
    /// units before the east edge (x = 600), shared with a neighbour.
    fn sheet() -> Layout {
        Layout::from_gcps(&[[0.0, 400.0, 0.0, 0.0], [600.0, 400.0, 600.0, 0.0], [0.0, 0.0, 0.0, 400.0], [600.0, 0.0, 600.0, 400.0]])
            .unwrap()
    }

    fn paper(x: f64) -> f64 {
        let dark = ((x - 570.0) / 30.0).clamp(0.0, 1.0) * 0.12;
        srgb_to_linear(1.0 - dark)
    }

    #[test]
    fn brings_edge_paper_to_white_and_fades_inwards() {
        let layout = sheet();
        let edge = Edge { axis: 0, pos: 600.0, inward: -1.0, span: (0.0, 400.0) };
        let field = EdgeField::measure(&layout, &[edge], &white_d65(), |px| {
            (px[0] >= 0.0 && px[0] < 600.0).then(|| [paper(px[0]); 3])
        });
        assert_eq!(field.edges(), 1);
        // In the measurement band (2 to 8 units from the edge) the corrected paper is ~white
        let corrected = |x: f64| paper(x) * field.log_gain([x, 200.0])[1].exp();
        let band: Vec<f64> = (592..598).map(|x| corrected(x as f64)).collect();
        let mean = band.iter().sum::<f64>() / band.len() as f64;
        assert!((mean - 1.0).abs() < 0.03, "band mean {mean}");
        // No correction beyond the reach (12 % of 600 = 72 units)
        assert_eq!(field.log_gain([500.0, 200.0]), [0.0; 3]);
    }

    #[test]
    fn corners_average_rather_than_add() {
        let layout = sheet();
        let east = Edge { axis: 0, pos: 600.0, inward: -1.0, span: (0.0, 400.0) };
        let north = Edge { axis: 1, pos: 400.0, inward: -1.0, span: (0.0, 600.0) };
        // Uniformly grey paper: both edges measure the same gain
        let grey = |_px: [f64; 2]| Some([0.5; 3]);
        let one = EdgeField::measure(&layout, &[east], &white_d65(), grey);
        let two = EdgeField::measure(&layout, &[east, north], &white_d65(), grey);
        // Pixel near the north-east corner (map y = 400 is pixel row 0)
        let (g1, g2) = (one.log_gain([599.0, 1.0]), two.log_gain([599.0, 1.0]));
        assert!((g1[1] - g2[1]).abs() < 1e-9, "{g1:?} vs {g2:?}");
    }
}
