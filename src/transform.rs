//! Per-pixel Bradford adaptation with a spatially varying paper white (flat-field).
//!
//! For each pixel: linear RGB -> Bradford LMS (one 3x3 matrix), multiply by the
//! gain LMS_target / LMS_paper(x, y) bilinearly interpolated from a low-res map,
//! back to linear RGB (one 3x3 matrix), then sRGB encoding.

use crate::color::{lab_to_xyz, mat_inv, mat_mul, mat_vec, Mat3, Vec3, M_BRADFORD, M_RGB2XYZ};

/// Low-resolution map of per-channel LMS gains.
#[derive(Clone)]
pub struct GainMap {
    pub width: usize,
    pub height: usize,
    pub data: Vec<[f32; 3]>,
}

impl GainMap {
    /// Gains mapping each paper colour (Lab, row-major `width` x `height`) to the target XYZ.
    pub fn from_paper_lab(paper: &[Vec3], width: usize, height: usize, target_xyz: &Vec3) -> Self {
        let lms_t = mat_vec(&M_BRADFORD, target_xyz);
        let data = paper
            .iter()
            .map(|lab| {
                let lms_p = mat_vec(&M_BRADFORD, &lab_to_xyz(lab));
                [0, 1, 2].map(|c| (lms_t[c] / lms_p[c]) as f32)
            })
            .collect();
        Self { width, height, data }
    }
}

/// Linear interpolation coordinates of output pixels in the low-res map
/// (pixel-centre aligned, clamped at the edges like `map_coordinates(mode="nearest")`).
fn lerp_coords(n_out: usize, n_in: usize) -> Vec<(usize, usize, f32)> {
    (0..n_out)
        .map(|i| {
            let x = ((i as f64 + 0.5) * n_in as f64 / n_out as f64 - 0.5).clamp(0.0, (n_in - 1) as f64);
            let i0 = x.floor() as usize;
            let i1 = (i0 + 1).min(n_in - 1);
            (i0, i1, (x - i0 as f64) as f32)
        })
        .collect()
}

/// Everything needed to transform rows of a full-resolution image.
pub struct Adapter {
    to_lms: [[f32; 3]; 3],
    from_lms: [[f32; 3]; 3],
    gains: GainMap,
    xs: Vec<(usize, usize, f32)>,
    ys: Vec<(usize, usize, f32)>,
}

fn to_f32(m: &Mat3) -> [[f32; 3]; 3] {
    m.map(|r| r.map(|v| v as f32))
}

impl Adapter {
    pub fn new(gains: GainMap, width: usize, height: usize) -> Self {
        let to_lms = mat_mul(&M_BRADFORD, &M_RGB2XYZ);
        let from_lms = mat_inv(&to_lms);
        let xs = lerp_coords(width, gains.width);
        let ys = lerp_coords(height, gains.height);
        Self { to_lms: to_f32(&to_lms), from_lms: to_f32(&from_lms), gains, xs, ys }
    }

    pub fn gains(&self) -> &GainMap {
        &self.gains
    }

    /// Multiply the gains of the low-res map by `f(x, y)`, evaluated at each cell centre in
    /// full-resolution pixel coordinates.
    pub fn scale_gains(&mut self, width: usize, height: usize, f: impl Fn([f64; 2]) -> [f64; 3]) {
        let (w, h) = (self.gains.width, self.gains.height);
        for (i, g) in self.gains.data.iter_mut().enumerate() {
            let x = ((i % w) as f64 + 0.5) * width as f64 / w as f64;
            let y = ((i / w) as f64 + 0.5) * height as f64 / h as f64;
            let s = f([x, y]);
            (0..3).for_each(|c| g[c] *= s[c] as f32);
        }
    }

    /// Adapt one linear RGB colour located at full-resolution pixel position `px`.
    pub fn adapt_point(&self, px: [f64; 2], rgb: [f64; 3]) -> [f64; 3] {
        let (w, h) = (self.gains.width, self.gains.height);
        let (fx, fy) = (px[0] * w as f64 / self.xs.len() as f64 - 0.5, px[1] * h as f64 / self.ys.len() as f64 - 0.5);
        let (fx, fy) = (fx.clamp(0.0, (w - 1) as f64), fy.clamp(0.0, (h - 1) as f64));
        let (x0, y0) = (fx as usize, fy as usize);
        let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
        let (tx, ty) = ((fx - x0 as f64) as f32, (fy - y0 as f64) as f32);
        let g = |x: usize, y: usize| self.gains.data[y * w + x];
        let gain: [f32; 3] = [0, 1, 2].map(|c| {
            let top = g(x0, y0)[c] + (g(x1, y0)[c] - g(x0, y0)[c]) * tx;
            let bottom = g(x0, y1)[c] + (g(x1, y1)[c] - g(x0, y1)[c]) * tx;
            top + (bottom - top) * ty
        });
        self.adapt(rgb.map(|c| c as f32), gain).map(|c| c as f64)
    }

    /// Bradford adaptation of one linear RGB colour with the given LMS gains.
    #[inline(always)]
    fn adapt(&self, [r, g, b]: [f32; 3], gain: [f32; 3]) -> [f32; 3] {
        let (p, q) = (&self.to_lms, &self.from_lms);
        let l = (p[0][0] * r + p[0][1] * g + p[0][2] * b) * gain[0];
        let m = (p[1][0] * r + p[1][1] * g + p[1][2] * b) * gain[1];
        let s = (p[2][0] * r + p[2][1] * g + p[2][2] * b) * gain[2];
        [
            q[0][0] * l + q[0][1] * m + q[0][2] * s,
            q[1][0] * l + q[1][1] * m + q[1][2] * s,
            q[2][0] * l + q[2][1] * m + q[2][2] * s,
        ]
    }

    /// Gains of the low-res map, interpolated vertically for full-resolution row `y`.
    pub fn row_gains(&self, y: usize) -> Vec<[f32; 3]> {
        let (y0, y1, fy) = self.ys[y];
        let w = self.gains.width;
        let (r0, r1) = (&self.gains.data[y0 * w..][..w], &self.gains.data[y1 * w..][..w]);
        r0.iter().zip(r1).map(|(a, b)| [0, 1, 2].map(|c| a[c] + (b[c] - a[c]) * fy)).collect()
    }

    /// Adapt one row of linear RGB pixels in place. `row_gains` comes from [`Self::row_gains`].
    #[inline]
    pub fn apply_row(&self, row: &mut [[f32; 3]], row_gains: &[[f32; 3]]) {
        for (px, &(x0, x1, fx)) in row.iter_mut().zip(&self.xs) {
            let (g0, g1) = (row_gains[x0], row_gains[x1]);
            *px = self.adapt(*px, [0, 1, 2].map(|c| g0[c] + (g1[c] - g0[c]) * fx));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::color::{linear_rgb_to_lab, srgb_to_lab, srgb_to_linear, white_d65};

    fn adapt(paper_srgb: Vec3, px_srgb: Vec3) -> Vec3 {
        let paper = srgb_to_lab(&paper_srgb);
        let gains = GainMap::from_paper_lab(&[paper], 1, 1, &white_d65());
        let ad = Adapter::new(gains, 1, 1);
        let mut row = [px_srgb.map(|c| srgb_to_linear(c) as f32)];
        ad.apply_row(&mut row, &ad.row_gains(0));
        linear_rgb_to_lab(&row[0].map(|c| c as f64))
    }

    #[test]
    fn paper_becomes_white() {
        let paper = [0.75, 0.68, 0.55];
        let lab = adapt(paper, paper);
        assert!((lab[0] - 100.0).abs() < 1e-3 && lab[1].abs() < 1e-3 && lab[2].abs() < 1e-3, "{lab:?}");
    }

    #[test]
    fn black_stays_black_and_grey_stays_neutral() {
        let paper = [0.75, 0.68, 0.55];
        assert!(adapt(paper, [0.0; 3])[0].abs() < 1e-6);
        // A grey with the paper's chromaticity becomes neutral
        let ink = paper.map(|c| srgb_to_linear(c) * 0.1).map(crate::color::linear_to_srgb);
        let lab = adapt(paper, ink);
        assert!(lab[1].abs() < 1e-2 && lab[2].abs() < 1e-2, "{lab:?}");
    }

    #[test]
    fn interpolation_coords() {
        // 2 low-res cells -> 4 output pixels: centres at -0.25, 0.25, 0.75, 1.25
        let c = lerp_coords(4, 2);
        assert_eq!(c[0], (0, 1, 0.0));
        assert_eq!(c[1], (0, 1, 0.25));
        assert_eq!(c[2], (0, 1, 0.75));
        assert_eq!(c[3], (1, 1, 0.0));
    }
}
