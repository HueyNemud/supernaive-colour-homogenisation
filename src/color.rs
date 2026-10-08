//! Colour conversions (sRGB, D65). Mirrors `proto/homog_proto.py`.

pub type Mat3 = [[f64; 3]; 3];
pub type Vec3 = [f64; 3];

pub const M_RGB2XYZ: Mat3 = [
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041],
];

pub const M_BRADFORD: Mat3 = [
    [0.8951, 0.2664, -0.1614],
    [-0.7502, 1.7135, 0.0367],
    [0.0389, -0.0685, 1.0296],
];

const DELTA: f64 = 6.0 / 29.0;

pub fn mat_vec(m: &Mat3, v: &Vec3) -> Vec3 {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

pub fn mat_mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    out
}

pub fn mat_inv(m: &Mat3) -> Mat3 {
    let [[a, b, c], [d, e, f], [g, h, i]] = *m;
    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    [
        [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]
}

pub fn white_d65() -> Vec3 {
    mat_vec(&M_RGB2XYZ, &[1.0, 1.0, 1.0])
}

pub fn srgb_to_linear(c: f64) -> f64 {
    if c <= 0.04045 { c / 12.92 } else { ((c + 0.055) / 1.055).powf(2.4) }
}

pub fn linear_to_srgb(c: f64) -> f64 {
    let c = c.clamp(0.0, 1.0);
    if c <= 0.0031308 { 12.92 * c } else { 1.055 * c.powf(1.0 / 2.4) - 0.055 }
}

fn f(t: f64) -> f64 {
    if t > DELTA * DELTA * DELTA { t.cbrt() } else { t / (3.0 * DELTA * DELTA) + 4.0 / 29.0 }
}

fn finv(t: f64) -> f64 {
    if t > DELTA { t * t * t } else { 3.0 * DELTA * DELTA * (t - 4.0 / 29.0) }
}

pub fn xyz_to_lab(xyz: &Vec3) -> Vec3 {
    let w = white_d65();
    let fx = f(xyz[0] / w[0]);
    let fy = f(xyz[1] / w[1]);
    let fz = f(xyz[2] / w[2]);
    [116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)]
}

pub fn lab_to_xyz(lab: &Vec3) -> Vec3 {
    let w = white_d65();
    let fy = (lab[0] + 16.0) / 116.0;
    let fx = fy + lab[1] / 500.0;
    let fz = fy - lab[2] / 200.0;
    [finv(fx) * w[0], finv(fy) * w[1], finv(fz) * w[2]]
}

pub fn linear_rgb_to_lab(rgb: &Vec3) -> Vec3 {
    xyz_to_lab(&mat_vec(&M_RGB2XYZ, rgb))
}

pub fn srgb_to_lab(rgb: &Vec3) -> Vec3 {
    linear_rgb_to_lab(&rgb.map(srgb_to_linear))
}

pub fn srgb_u8_to_linear(rgb: [u8; 3]) -> Vec3 {
    rgb.map(|c| srgb_to_linear(c as f64 / 255.0))
}

pub fn srgb_to_xyz_u8(rgb: [u8; 3]) -> Vec3 {
    mat_vec(&M_RGB2XYZ, &srgb_u8_to_linear(rgb))
}

pub fn lab_to_srgb(lab: &Vec3) -> Vec3 {
    mat_vec(&mat_inv(&M_RGB2XYZ), &lab_to_xyz(lab)).map(linear_to_srgb)
}

/// Fast linear -> sRGB encoder: table indexed by sqrt(linear), linearly interpolated.
/// Max error < 1e-6 over [0, 1], i.e. well below one 16-bit code value.
pub struct SrgbEncoder {
    table: Vec<f32>,
}

impl SrgbEncoder {
    const N: usize = 4096;

    pub fn new() -> Self {
        let table = (0..=Self::N + 1)
            .map(|i| {
                let s = (i as f64 / Self::N as f64).min(1.0);
                linear_to_srgb(s * s) as f32
            })
            .collect();
        Self { table }
    }

    #[inline(always)]
    pub fn encode(&self, lin: f32) -> f32 {
        let x = lin.clamp(0.0, 1.0).sqrt() * Self::N as f32;
        let i = x as usize;
        let t = x - i as f32;
        let (a, b) = (self.table[i], self.table[i + 1]);
        a + (b - a) * t
    }
}

impl Default for SrgbEncoder {
    fn default() -> Self {
        Self::new()
    }
}

/// sRGB-encoded integer -> linear, for every code value of a given bit depth.
pub fn decode_table(max_code: u32) -> Vec<f32> {
    (0..=max_code)
        .map(|v| srgb_to_linear(v as f64 / max_code as f64) as f32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lab_roundtrip() {
        for rgb in [[0.2, 0.5, 0.7], [1.0, 1.0, 1.0], [0.01, 0.0, 0.3], [0.9, 0.8, 0.6]] {
            let back = lab_to_srgb(&srgb_to_lab(&rgb));
            for c in 0..3 {
                assert!((back[c] - rgb[c]).abs() < 1e-9, "{rgb:?} -> {back:?}");
            }
        }
    }

    #[test]
    fn white_is_l100() {
        let lab = srgb_to_lab(&[1.0, 1.0, 1.0]);
        assert!((lab[0] - 100.0).abs() < 1e-9 && lab[1].abs() < 1e-9 && lab[2].abs() < 1e-9);
    }

    #[test]
    fn encoder_accuracy() {
        let enc = SrgbEncoder::new();
        let mut max_err: f64 = 0.0;
        for i in 0..=1_000_000 {
            let lin = i as f64 / 1_000_000.0;
            max_err = max_err.max((enc.encode(lin as f32) as f64 - linear_to_srgb(lin)).abs());
        }
        assert!(max_err < 0.5 / 65535.0, "max error {max_err}");
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn inverse() {
        let p = mat_mul(&M_BRADFORD, &mat_inv(&M_BRADFORD));
        for i in 0..3 {
            for j in 0..3 {
                assert!((p[i][j] - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
    }
}
