//! Seam correction between adjacent sheets.
//!
//! The paper surface cannot follow the darkening close to the edges of a sheet, which shows
//! as a step where two sheets meet. Along every edge shared with another sheet, the paper is
//! measured (after the surface correction) in a thin band inside the neatline, per segment,
//! and brought back to the target by a Bradford-LMS gain that fades out inwards. The band
//! starts beyond the GCP inaccuracy and the fading distance is measured on the sheet.

use crate::background::median;
use crate::color::{mat_inv, mat_vec, white_d65, Vec3, M_BRADFORD, M_RGB2XYZ};
use crate::layout::{Edge, Layout};

/// Width of the measurement band (fraction of the sheet size), and its minimal start
const BAND_WIDTH: f64 = 0.01;
const BAND_MIN_START: f64 = 0.003;
/// Bounds of the fading distance (fractions of the sheet size)
const REACH_BOUNDS: (f64, f64) = (0.04, 0.20);
/// Profile of the darkening: distance tranches up to this fraction of the sheet size
const PROFILE_DEPTH: f64 = 0.25;
const PROFILE_TRANCHES: usize = 12;
const BINS: usize = 24;
const ALONG: usize = 400;
const ACROSS: usize = 8;

/// Factor applied to the luminance Y of a neutral colour by the Bradford-LMS gains `g`.
pub fn luminance_gain(g: &[f64; 3]) -> f64 {
    let w = white_d65();
    let lms = mat_vec(&M_BRADFORD, &w);
    let scaled = [lms[0] * g[0], lms[1] * g[1], lms[2] * g[2]];
    mat_vec(&mat_inv(&M_BRADFORD), &scaled)[1] / w[1]
}

/// A measurement point: full-resolution pixel position and whether it was taken as paper.
pub type Sample = ([f64; 2], bool);

/// Colours at a full-resolution pixel position: (as scanned, after the surface correction),
/// both linear RGB; None outside the image.
pub trait Sampler: Fn([f64; 2]) -> Option<([f64; 3], [f64; 3])> {}
impl<F: Fn([f64; 2]) -> Option<([f64; 3], [f64; 3])>> Sampler for F {}

struct EdgeProfile {
    edge: Edge,
    /// Log LMS gain per segment along the edge
    gains: Vec<[f64; 3]>,
    /// Distance over which the correction fades out
    reach: f64,
}

/// Log-gain field of one sheet.
pub struct EdgeField {
    layout: Layout,
    profiles: Vec<EdgeProfile>,
    /// Measurement band (start, end), in map units
    pub band: (f64, f64),
    /// All band measurement points (for the debug images)
    pub samples: Vec<Sample>,
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

fn to_lms(rgb: &[f64; 3]) -> Vec3 {
    mat_vec(&M_BRADFORD, &mat_vec(&M_RGB2XYZ, rgb))
}

/// Map position at distance `t` inside `edge`, at position `u` along it.
fn point(edge: &Edge, t: f64, u: f64) -> [f64; 2] {
    let d = edge.pos + edge.inward * t;
    if edge.axis == 0 { [d, u] } else { [u, d] }
}

impl EdgeField {
    /// Measure the paper along `edges`. `is_paper(scanned, corrected)` tells whether a pair of
    /// linear RGB colours is bare paper.
    pub fn measure(
        layout: &Layout,
        edges: &[Edge],
        target_xyz: &Vec3,
        is_paper: impl Fn(&[f64; 3], &[f64; 3]) -> bool,
        sample: impl Sampler,
    ) -> Self {
        let size = layout.size();
        let lms_t = mat_vec(&M_BRADFORD, target_xyz);
        let log_gain = |lms: &[Vec3]| [0, 1, 2].map(|c| (lms_t[c] / median(lms.iter().map(|p| p[c]).collect())).ln());
        let start = (BAND_MIN_START * size).max(2.0 * layout.rms);
        let band = (start, start + BAND_WIDTH * size);
        let mut samples = Vec::new();
        let profiles = edges
            .iter()
            .filter_map(|&edge| {
                let (u0, u1) = edge.span;
                // Paper in the band, per segment along the edge
                let mut bins: Vec<Vec<Vec3>> = vec![Vec::new(); BINS];
                for a in 0..ALONG {
                    let f = a as f64 / (ALONG - 1) as f64;
                    for k in 0..ACROSS {
                        let t = band.0 + (band.1 - band.0) * k as f64 / (ACROSS - 1) as f64;
                        let px = layout.pixel(point(&edge, t, u0 + (u1 - u0) * f));
                        let Some((raw, corrected)) = sample(px) else { continue };
                        let paper = is_paper(&raw, &corrected);
                        if paper {
                            bins[((f * BINS as f64) as usize).min(BINS - 1)].push(to_lms(&corrected));
                        }
                        samples.push((px, paper));
                    }
                }
                let known: Vec<(f64, [f64; 3])> = bins
                    .iter()
                    .enumerate()
                    .filter(|(_, b)| b.len() >= 10)
                    .map(|(k, b)| ((k as f64 + 0.5) / BINS as f64, log_gain(b)))
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
                let gains: Vec<[f64; 3]> = (0..BINS)
                    .map(|k| {
                        let idx = [k.saturating_sub(1), k, (k + 1).min(BINS - 1)];
                        [0, 1, 2].map(|c| idx.iter().map(|&i| filled[i][c]).sum::<f64>() / 3.0)
                    })
                    .collect();
                // Darkening profile: paper gain per distance tranche; the correction fades out
                // where the remaining gain falls below a quarter of the gain at the edge
                let at_edge = luminance_gain(&log_gain(&bins.concat()).map(f64::exp)).ln();
                let mut reach = REACH_BOUNDS.1 * size;
                if at_edge.abs() < 0.01 {
                    reach = REACH_BOUNDS.0 * size;
                } else {
                    for k in 0..PROFILE_TRANCHES {
                        let t = band.1 + (PROFILE_DEPTH * size - band.1) * (k as f64 + 0.5) / PROFILE_TRANCHES as f64;
                        let paper: Vec<Vec3> = (0..ALONG / 4)
                            .filter_map(|a| {
                                let u = u0 + (u1 - u0) * a as f64 / (ALONG / 4 - 1) as f64;
                                let (raw, corrected) = sample(layout.pixel(point(&edge, t, u)))?;
                                is_paper(&raw, &corrected).then(|| to_lms(&corrected))
                            })
                            .collect();
                        if paper.len() >= 10 {
                            let g = luminance_gain(&log_gain(&paper).map(f64::exp)).ln();
                            if g.abs() < 0.25 * at_edge.abs() || g.signum() != at_edge.signum() {
                                reach = t;
                                break;
                            }
                        }
                    }
                }
                let reach = reach.clamp(REACH_BOUNDS.0 * size, REACH_BOUNDS.1 * size);
                Some(EdgeProfile { edge, gains, reach })
            })
            .collect();
        Self { layout: layout.clone(), profiles, band, samples }
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
                    "{} edge {} = {:.1} (span {:.1}..{:.1}), fading over {:.1} map units",
                    if e.axis == 0 { "vertical" } else { "horizontal" },
                    if e.axis == 0 { "x" } else { "y" },
                    e.pos,
                    e.span.0,
                    e.span.1,
                    p.reach
                );
                (name, p.gains.iter().map(|g| luminance_gain(&g.map(f64::exp))).collect())
            })
            .collect()
    }

    /// Log LMS gain at a full-resolution pixel position: full correction up to the middle of
    /// the measurement band, fading out (squared ramp) at the edge's reach. Near a corner,
    /// the corrections of the two edges measure the same darkening: they are averaged
    /// rather than added (weights normalised when they sum to more than 1).
    pub fn log_gain(&self, px: [f64; 2]) -> [f64; 3] {
        let xy = self.layout.map(px);
        let mut g = [0.0; 3];
        let mut weight = 0.0;
        for p in &self.profiles {
            let (across, along) = if p.edge.axis == 0 { (xy[0], xy[1]) } else { (xy[1], xy[0]) };
            // Full correction up to the middle of the measurement band, then fading out
            let d = p.edge.inward * (across - p.edge.pos);
            let mid = (self.band.0 + self.band.1) / 2.0;
            let phi = if d <= mid { 1.0 } else { (1.0 - (d - mid) / (p.reach - mid)).clamp(0.0, 1.0).powi(2) };
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
    use crate::color::srgb_to_linear;

    /// Sheet 600 x 400 map units, 1 px per unit, y up.
    fn sheet() -> Layout {
        Layout::from_gcps(&[[0.0, 400.0, 0.0, 0.0], [600.0, 400.0, 600.0, 0.0], [0.0, 0.0, 0.0, 400.0], [600.0, 0.0, 600.0, 400.0]])
            .unwrap()
    }

    /// Paper darkening over the last `width` units before the east edge (x = 600).
    fn paper(x: f64, width: f64) -> f64 {
        let dark = ((x - (600.0 - width)) / width).clamp(0.0, 1.0) * 0.12;
        srgb_to_linear(1.0 - dark)
    }

    fn east() -> Edge {
        Edge { axis: 0, pos: 600.0, inward: -1.0, span: (0.0, 400.0) }
    }

    fn grey(rgb: &[f64; 3], _: &[f64; 3]) -> bool {
        rgb[0] > 0.5
    }

    fn field(width: f64) -> EdgeField {
        EdgeField::measure(&sheet(), &[east()], &white_d65(), grey, move |px: [f64; 2]| {
            (px[0] >= 0.0 && px[0] < 600.0).then(|| ([paper(px[0], width); 3], [paper(px[0], width); 3]))
        })
    }

    #[test]
    fn brings_edge_paper_to_white_and_fades_inwards() {
        let f = field(30.0);
        assert_eq!(f.edges(), 1);
        // In the measurement band the corrected paper is ~white
        let corrected = |x: f64| paper(x, 30.0) * f.log_gain([x, 200.0])[1].exp();
        let band: Vec<f64> = (592..597).map(|x| corrected(x as f64)).collect();
        let mean = band.iter().sum::<f64>() / band.len() as f64;
        assert!((mean - 1.0).abs() < 0.03, "band mean {mean}");
        // No correction far from the edge
        assert_eq!(f.log_gain([300.0, 200.0]), [0.0; 3]);
    }

    #[test]
    fn fading_distance_follows_the_darkening_width() {
        // Darkening over 60 then 100 units: the fading distance grows accordingly
        let (r60, r100) = (field(60.0).profiles[0].reach, field(100.0).profiles[0].reach);
        assert!(r100 > r60 * 1.3, "reach {r60} then {r100}");
        assert!((40.0..=90.0).contains(&r60), "reach {r60} for a 60-unit darkening");
    }

    #[test]
    fn band_starts_beyond_the_gcp_inaccuracy() {
        let mut layout = sheet();
        layout.rms = 5.0;
        let f = EdgeField::measure(&layout, &[east()], &white_d65(), grey, |_px: [f64; 2]| Some(([0.9; 3], [0.9; 3])));
        assert!((f.band.0 - 10.0).abs() < 1e-9, "band {:?}", f.band);
    }

    #[test]
    fn corners_average_rather_than_add() {
        let layout = sheet();
        let north = Edge { axis: 1, pos: 400.0, inward: -1.0, span: (0.0, 600.0) };
        // Uniformly grey paper: both edges measure the same gain
        let s = |_px: [f64; 2]| Some(([0.5; 3], [0.5; 3]));
        let paper = |_: &[f64; 3], _: &[f64; 3]| true;
        let one = EdgeField::measure(&layout, &[east()], &white_d65(), paper, s);
        let two = EdgeField::measure(&layout, &[east(), north], &white_d65(), paper, s);
        // Pixel near the north-east corner (map y = 400 is pixel row 0)
        let (g1, g2) = (one.log_gain([599.0, 1.0]), two.log_gain([599.0, 1.0]));
        assert!((g1[1] - g2[1]).abs() < 1e-9, "{g1:?} vs {g2:?}");
    }
}
