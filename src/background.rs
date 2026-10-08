//! Paper (background) estimation on a low-resolution Lab thumbnail: paper colour model
//! (automatic, or from `paper` imagettes), content mask, and smooth paper surface.

use anyhow::{bail, Result};
use rayon::prelude::*;

use crate::color::Vec3;
use crate::hints::{ColourModel, Hints, MATCH_RADIUS};

/// Row-major Lab thumbnail.
pub struct LabImage {
    pub width: usize,
    pub height: usize,
    pub data: Vec<Vec3>,
}

/// numpy-compatible percentile (linear interpolation). Sorts `v` in place.
fn percentile(v: &mut [f64], q: f64) -> f64 {
    v.sort_unstable_by(f64::total_cmp);
    let pos = q / 100.0 * (v.len() - 1) as f64;
    let i = pos.floor() as usize;
    let j = (i + 1).min(v.len() - 1);
    v[i] + (v[j] - v[i]) * (pos - i as f64)
}

pub(crate) fn median(mut v: Vec<f64>) -> f64 {
    percentile(&mut v, 50.0)
}

fn median3(px: &[Vec3]) -> Vec3 {
    [0, 1, 2].map(|c| median(px.iter().map(|p| p[c]).collect()))
}

fn ab_dist(p: &Vec3, q: &Vec3) -> f64 {
    ((p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// How uneven the lighting of the scans is, i.e. how flexible the paper surface is.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Lighting {
    /// Degree chosen per sheet by spatial cross-validation
    Auto,
    /// Degree 2
    Even,
    /// Degree 4
    Uneven,
    /// Single paper colour per sheet (no flat-field)
    None,
}

/// Floor of the a*b* spread of the automatic paper model: paper within 8 Δab of the
/// estimated colour is always accepted.
const AUTO_MIN_SD_AB: f64 = 8.0 / MATCH_RADIUS;

/// Automatic paper colour: the bright pixels (75th-98th percentile of L*) with the dominant
/// chroma, among `pixels` (content, not matching a `keep` hint).
fn auto_paper(pixels: Vec<Vec3>) -> Option<ColourModel> {
    if pixels.len() < 50 {
        return None;
    }
    let mut ls: Vec<f64> = pixels.iter().map(|p| p[0]).collect();
    let lo = percentile(&mut ls, 75.0);
    let hi = percentile(&mut ls, 98.0);
    let sel: Vec<Vec3> = pixels.into_iter().filter(|p| p[0] >= lo && p[0] <= hi).collect();
    let ab_med = [0.0, median(sel.iter().map(|p| p[1]).collect()), median(sel.iter().map(|p| p[2]).collect())];
    let sel: Vec<Vec3> = sel.into_iter().filter(|p| ab_dist(p, &ab_med) < 6.0).collect();
    ColourModel::from_pixels(&sel, AUTO_MIN_SD_AB)
}

fn delta_e(p: &Vec3, q: &Vec3) -> f64 {
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// 4-connected flood fill from `seeds` through the pixels where `inside` holds.
fn flood(w: usize, h: usize, seeds: impl IntoIterator<Item = usize>, inside: impl Fn(usize) -> bool) -> Vec<bool> {
    let mut reached = vec![false; w * h];
    let mut stack: Vec<usize> = seeds.into_iter().collect();
    while let Some(i) = stack.pop() {
        if reached[i] || !inside(i) {
            continue;
        }
        reached[i] = true;
        let (x, y) = (i % w, i / w);
        if x > 0 {
            stack.push(i - 1);
        }
        if x + 1 < w {
            stack.push(i + 1);
        }
        if y > 0 {
            stack.push(i - w);
        }
        if y + 1 < h {
            stack.push(i + w);
        }
    }
    reached
}

/// Pixels of the sheet, and the colour of the scanner bed if one was found. The bed is the
/// uniform colour touching the image border: it is flood-filled from the border, then only
/// the largest remaining region (the sheet) is kept, which also drops colour charts or rulers
/// lying on the bed. If the "bed" covers most of the image, the scan is cropped to the map
/// and everything is content.
pub fn content_mask(img: &LabImage) -> (Vec<bool>, Option<Vec3>) {
    let (w, h) = (img.width, img.height);
    let border: Vec<usize> = (0..w).flat_map(|x| [x, (h - 1) * w + x]).chain((0..h).flat_map(|y| [y * w, y * w + w - 1])).collect();
    let bed = median3(&border.iter().map(|&i| img.data[i]).collect::<Vec<_>>());
    let background = flood(w, h, border, |i| delta_e(&img.data[i], &bed) < 15.0);
    if background.iter().filter(|&&b| b).count() * 10 > w * h * 6 {
        return (vec![true; w * h], None);
    }
    // Largest connected region of the rest (single-pass labelling)
    let mut label = vec![usize::MAX; w * h];
    let mut best = (0, usize::MAX);
    let mut stack = Vec::new();
    for start in 0..w * h {
        if background[start] || label[start] != usize::MAX {
            continue;
        }
        let mut size = 0;
        stack.push(start);
        label[start] = start;
        while let Some(i) = stack.pop() {
            size += 1;
            let (x, y) = (i % w, i / w);
            let neighbours = [(x > 0).then(|| i - 1), (x + 1 < w).then(|| i + 1), (y > 0).then(|| i - w), (y + 1 < h).then(|| i + w)];
            for j in neighbours.into_iter().flatten() {
                if !background[j] && label[j] == usize::MAX {
                    label[j] = start;
                    stack.push(j);
                }
            }
        }
        if size > best.0 {
            best = (size, start);
        }
    }
    (label.iter().map(|&l| l == best.1).collect(), Some(bed))
}

/// Paper of a sheet: reference colour models (automatic estimate or `paper` imagettes),
/// global paper colour and content mask.
pub struct PaperModel {
    pub global: Vec3,
    pub content: Vec<bool>,
    refs: Vec<ColourModel>,
}

impl PaperModel {
    pub fn estimate(img: &LabImage, hints: &Hints) -> Result<Self> {
        let usable = |mask: &dyn Fn(usize) -> bool| -> Vec<Vec3> {
            (0..img.data.len()).filter(|&i| mask(i) && !hints.keeps(&img.data[i])).map(|i| img.data[i]).collect()
        };
        let (mut content, bed) = content_mask(img);
        // A "bed" of the paper colour is paper touching the image border (scan cropped to the
        // map, cells closed by ink lines): no background then
        let is_paper = |bed: &Vec3, paper: &Vec3| delta_e(bed, paper) < 10.0;
        if hints.paper.is_empty() {
            let auto = |content: &[bool]| auto_paper(usable(&|i| content[i]));
            let Some(mut model) = auto(&content) else { bail!("no paper found: provide `paper` imagettes") };
            if bed.is_some_and(|b| is_paper(&b, &model.mean)) {
                content = vec![true; img.data.len()];
                model = auto(&content).expect("more pixels than before");
            }
            return Ok(Self { global: model.mean, content, refs: vec![model] });
        }
        if bed.is_some_and(|b| hints.paper.iter().any(|m| is_paper(&b, &m.mean))) {
            content = vec![true; img.data.len()];
        }
        // Imagettes: global colour = median of the matching pixels
        let mut model = Self { global: hints.paper[0].mean, content, refs: hints.paper.clone() };
        let matching: Vec<Vec3> = usable(&|i| model.content[i]).into_iter().filter(|p| model.accepts(p)).collect();
        if matching.len() < 50 {
            bail!("the `paper` imagettes match almost no pixel of the sheet");
        }
        model.global = median3(&matching);
        Ok(model)
    }

    /// Whether a colour is compatible with one of the paper references (any lighting).
    pub fn accepts(&self, lab: &Vec3) -> bool {
        self.refs.iter().any(|m| {
            let dl = lab[0] - m.mean[0];
            m.distance_ab(lab) < MATCH_RADIUS && dl > -(3.0 * m.sd_l).max(15.0) && dl < (2.0 * m.sd_l).max(10.0)
        })
    }

    /// Largest a*b* spread of the paper references.
    pub fn sd_ab(&self) -> f64 {
        self.refs.iter().map(|m| m.sd_ab).fold(0.0, f64::max)
    }

    /// Indices of the pixels used to fit the paper surface.
    fn candidates(&self, img: &LabImage, hints: &Hints) -> Vec<usize> {
        (0..img.data.len())
            .filter(|&i| self.content[i] && self.accepts(&img.data[i]) && !hints.keeps(&img.data[i]))
            .collect()
    }
}

/// Monomials x^i y^j (i + j <= degree) of the pixel centre, coordinates normalised to [-1, 1].
fn monomials(x: usize, y: usize, w: usize, h: usize, degree: usize, out: &mut [f64]) {
    let xn = (x as f64 + 0.5) / w as f64 * 2.0 - 1.0;
    let yn = (y as f64 + 0.5) / h as f64 * 2.0 - 1.0;
    let mut k = 0;
    for i in 0..=degree {
        for j in 0..=degree - i {
            out[k] = xn.powi(i as i32) * yn.powi(j as i32);
            k += 1;
        }
    }
}

/// Weighted least squares: coefficients minimising sum_i w_i (x_i . beta - y_i)^2 for each
/// of the `Y` targets, by normal equations and Cholesky decomposition.
pub(crate) fn weighted_lstsq<const Y: usize>(x: &[f64], k: usize, y: &[[f64; Y]], w: &[f64]) -> [Vec<f64>; Y] {
    let zero = || (vec![0.0; k * k], vec![[0.0; Y]; k]);
    let (mut a, b) = x
        .par_chunks_exact(k)
        .zip(y)
        .zip(w)
        .fold(zero, |(mut a, mut b), ((xi, yi), &wi)| {
            for r in 0..k {
                let wx = wi * xi[r];
                for c in 0..=r {
                    a[r * k + c] += wx * xi[c];
                }
                for t in 0..Y {
                    b[r][t] += wx * yi[t];
                }
            }
            (a, b)
        })
        .reduce(zero, |(mut a1, mut b1), (a2, b2)| {
            a1.iter_mut().zip(&a2).for_each(|(u, v)| *u += v);
            b1.iter_mut().zip(&b2).for_each(|(u, v)| (0..Y).for_each(|t| u[t] += v[t]));
            (a1, b1)
        });
    // Tiny scale-invariant ridge for numerical safety, then in-place Cholesky (lower triangle)
    for i in 0..k {
        a[i * k + i] *= 1.0 + 1e-12;
        a[i * k + i] += 1e-300;
    }
    for j in 0..k {
        let d = (a[j * k + j] - (0..j).map(|p| a[j * k + p].powi(2)).sum::<f64>()).sqrt();
        a[j * k + j] = d;
        for i in j + 1..k {
            a[i * k + j] = (a[i * k + j] - (0..j).map(|p| a[i * k + p] * a[j * k + p]).sum::<f64>()) / d;
        }
    }
    std::array::from_fn(|t| {
        let mut z = vec![0.0; k];
        for i in 0..k {
            z[i] = (b[i][t] - (0..i).map(|p| a[i * k + p] * z[p]).sum::<f64>()) / a[i * k + i];
        }
        for i in (0..k).rev() {
            z[i] = (z[i] - (i + 1..k).map(|p| a[p * k + i] * z[p]).sum::<f64>()) / a[i * k + i];
        }
        z
    })
}

pub(crate) fn eval(beta: &[f64], xi: &[f64]) -> f64 {
    beta.iter().zip(xi).map(|(b, x)| b * x).sum()
}

/// What a thumbnail pixel was used for in the paper surface fit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PixelRole {
    /// Outside the sheet (scanner bed, colour chart...)
    Background,
    /// Not compatible with the paper colour (ink, strong wash) or matching a `keep` hint
    Other,
    /// Used for the L* envelope only
    Lightness,
    /// Used for the L* envelope and for a*, b* (bare paper)
    Chroma,
}

/// Paper surface: one Lab value per thumbnail pixel, the role of each pixel in the fit and
/// the polynomial degree.
pub struct PaperSurface {
    pub lab: Vec<Vec3>,
    pub role: Vec<PixelRole>,
    pub degree: Option<usize>,
}

/// Envelope quantile and IRLS settings
const TAU: f64 = 0.8;
const NEAR: f64 = 4.0;
const C_AB: f64 = 3.0;

fn n_terms(degree: usize) -> usize {
    (degree + 1) * (degree + 2) / 2
}

/// Upper-envelope fit of `l` (quantile regression, quantile TAU) by IRLS.
fn fit_envelope(x: &[f64], k: usize, l: &[[f64; 1]], iterations: usize) -> Vec<f64> {
    let mut wt = vec![1.0; l.len()];
    let mut beta = Vec::new();
    for _ in 0..=iterations {
        [beta] = weighted_lstsq(x, k, l, &wt);
        for ((xi, li), wi) in x.chunks_exact(k).zip(l).zip(wt.iter_mut()) {
            let r = li[0] - eval(&beta, xi);
            *wi = if r > 0.0 { TAU } else { 1.0 - TAU } / r.abs().max(0.05);
        }
    }
    beta
}

fn pinball(r: f64) -> f64 {
    if r > 0.0 { TAU * r } else { (TAU - 1.0) * r }
}

/// Degree of the paper surface for `lighting = auto`: spatial 4-fold cross-validation of the
/// envelope fit (folds = interleaved blocks of a 4 x 4 grid), smallest degree whose held-out
/// quantile loss is within 1 % of the best. A stiffer surface protects the washes.
pub fn select_degree(img: &LabImage, cand: &[usize]) -> usize {
    const DEGREES: std::ops::RangeInclusive<usize> = 1..=5;
    let (w, h) = (img.width, img.height);
    let step = (cand.len() / 20_000).max(1);
    let pts: Vec<usize> = cand.iter().step_by(step).copied().collect();
    let fold = |i: usize| {
        let (bx, by) = ((i % w) * 4 / w, (i / w) * 4 / h);
        (bx % 2) + 2 * (by % 2)
    };
    let losses: Vec<f64> = DEGREES
        .into_par_iter()
        .map(|d| {
            let k = n_terms(d);
            let mut row = vec![0.0; k];
            let (mut loss, mut n) = (0.0, 0.0);
            for f in 0..4 {
                let (train, test): (Vec<usize>, Vec<usize>) = pts.iter().partition(|&&i| fold(i) != f);
                if train.len() < 20 * k || test.is_empty() {
                    continue;
                }
                let x: Vec<f64> = train
                    .iter()
                    .flat_map(|&i| {
                        monomials(i % w, i / w, w, h, d, &mut row);
                        row.clone()
                    })
                    .collect();
                let l: Vec<[f64; 1]> = train.iter().map(|&i| [img.data[i][0]]).collect();
                let beta = fit_envelope(&x, k, &l, 15);
                for &i in &test {
                    monomials(i % w, i / w, w, h, d, &mut row);
                    loss += pinball(img.data[i][0] - eval(&beta, &row));
                    n += 1.0;
                }
            }
            if n > 0.0 { loss / n } else { f64::INFINITY }
        })
        .collect();
    let best = losses.iter().copied().fold(f64::INFINITY, f64::min);
    DEGREES.zip(&losses).find(|(_, l)| **l <= best * 1.01).map_or(3, |(d, _)| d)
}

/// Paper surface: smooth polynomial surfaces of the given degree (None: single colour),
/// robust to large washes.
///
/// Bare paper is the lightest material (a transparent wash or ink can only darken it):
/// - L*: upper-envelope fit, i.e. quantile regression (80 %) by iteratively reweighted
///   least squares, over the pixels compatible with the paper model;
/// - a*, b*: least squares with Tukey reweighting (cut-off 3 Δab units, starting from the
///   global paper chroma) over the pixels at most 4 L* units below the envelope, so that
///   washes do not contribute.
///
/// A low-degree surface cannot follow a wash patch, however large.
pub fn paper_surface(img: &LabImage, paper: &PaperModel, hints: &Hints, degree: Option<usize>) -> PaperSurface {
    let (w, h) = (img.width, img.height);
    let global = paper.global;
    let cand = paper.candidates(img, hints);
    let mut role: Vec<PixelRole> =
        paper.content.iter().map(|&c| if c { PixelRole::Other } else { PixelRole::Background }).collect();
    cand.iter().for_each(|&i| role[i] = PixelRole::Lightness);
    let flat = |role| PaperSurface { lab: vec![global; w * h], role, degree: None };
    let Some(degree) = degree else { return flat(role) };
    let k = n_terms(degree);
    if cand.len() < 20 * k {
        return flat(role);
    }

    let mut design = vec![0.0; w * h * k];
    design.par_chunks_exact_mut(k).enumerate().for_each(|(i, row)| monomials(i % w, i / w, w, h, degree, row));
    let rows = |idx: &[usize]| -> Vec<f64> { idx.iter().flat_map(|&i| design[i * k..][..k].iter().copied()).collect() };

    let l: Vec<[f64; 1]> = cand.iter().map(|&i| [img.data[i][0]]).collect();
    let beta_l = fit_envelope(&rows(&cand), k, &l, 30);
    let l_surf: Vec<f64> = design.par_chunks_exact(k).map(|xi| eval(&beta_l, xi)).collect();

    let near: Vec<usize> = cand.into_iter().filter(|&i| img.data[i][0] - l_surf[i] > -NEAR).collect();
    let (beta_a, beta_b) = if near.len() < 20 * k {
        (None, None)
    } else {
        let xn = rows(&near);
        let ab: Vec<[f64; 2]> = near.iter().map(|&i| [img.data[i][1], img.data[i][2]]).collect();
        let mut pred: Vec<[f64; 2]> = vec![[global[1], global[2]]; near.len()];
        let mut wt = vec![0.0; near.len()];
        let mut beta = [Vec::new(), Vec::new()];
        for _ in 0..10 {
            for ((wi, p), q) in wt.iter_mut().zip(&ab).zip(&pred) {
                let u = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt() / C_AB;
                *wi = if u < 1.0 { (1.0 - u * u).powi(2) } else { 0.0 };
            }
            beta = weighted_lstsq(&xn, k, &ab, &wt);
            for (q, xi) in pred.iter_mut().zip(xn.chunks_exact(k)) {
                *q = [eval(&beta[0], xi), eval(&beta[1], xi)];
            }
        }
        near.iter().zip(&wt).filter(|(_, w)| **w > 0.0).for_each(|(&i, _)| role[i] = PixelRole::Chroma);
        let [a, b] = beta;
        (Some(a), Some(b))
    };

    let lab = design
        .par_chunks_exact(k)
        .zip(&l_surf)
        .map(|(xi, &ls)| {
            // Keep the extrapolation (outside the paper, e.g. the scanner background) sane
            let ls = ls.clamp((global[0] - 30.0).max(5.0), 100.0);
            let a = beta_a.as_ref().map_or(global[1], |b| eval(b, xi));
            let b = beta_b.as_ref().map_or(global[2], |b| eval(b, xi));
            [ls, a, b]
        })
        .collect();
    PaperSurface { lab, role, degree: Some(degree) }
}

/// Degree of the paper surface for a lighting setting (None: single paper colour).
pub fn surface_degree(lighting: Lighting, img: &LabImage, paper: &PaperModel, hints: &Hints) -> Option<usize> {
    match lighting {
        Lighting::Auto => Some(select_degree(img, &paper.candidates(img, hints))),
        Lighting::Even => Some(2),
        Lighting::Uneven => Some(4),
        Lighting::None => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percentile_matches_numpy() {
        // np.percentile([1, 2, 3, 4], [50, 75, 98]) == [2.5, 3.25, 3.94]
        let mut v = vec![4.0, 1.0, 3.0, 2.0];
        assert!((percentile(&mut v, 50.0) - 2.5).abs() < 1e-12);
        assert!((percentile(&mut v, 75.0) - 3.25).abs() < 1e-12);
        assert!((percentile(&mut v, 98.0) - 3.94).abs() < 1e-12);
    }

    fn surface(img: &LabImage, degree: usize) -> Vec<Vec3> {
        let hints = Hints::default();
        let paper = PaperModel::estimate(img, &hints).unwrap();
        paper_surface(img, &paper, &hints, Some(degree)).lab
    }

    #[test]
    fn paper_found_despite_ink() {
        // Yellowish paper with dark ink stripes and a pale blue wash
        let (w, h) = (100, 80);
        let paper = [70.0, 1.0, 16.0];
        let data = (0..w * h)
            .map(|i| match i % w {
                x if x % 10 == 0 => [20.0, 0.0, 2.0],
                x if (40..60).contains(&x) => [72.0, -8.0, -10.0],
                _ => paper,
            })
            .collect();
        let img = LabImage { width: w, height: h, data };
        let g = PaperModel::estimate(&img, &Hints::default()).unwrap().global;
        assert!((0..3).all(|c| (g[c] - paper[c]).abs() < 1e-9), "{g:?}");
        let est = surface(&img, 3);
        assert!(est.iter().all(|p| (0..3).all(|c| (p[c] - paper[c]).abs() < 0.5)), "{:?}", est[0]);
    }

    /// Paper darkening from left to right (illumination) and turning yellower at the
    /// bottom, with a large green wash (40 % of the sheet) slightly darker than the paper.
    fn sheet_with_wash(x: usize, y: usize, w: usize, h: usize) -> (Vec3, bool) {
        let paper = [75.0 - 15.0 * x as f64 / w as f64, 1.0, 14.0 + 4.0 * y as f64 / h as f64];
        let in_wash = (w / 5..w * 3 / 4).contains(&x) && (h / 6..h * 5 / 6).contains(&y);
        if in_wash { ([paper[0] - 4.0, paper[1] - 4.0, paper[2] - 2.0], true) } else { (paper, false) }
    }

    #[test]
    fn surface_follows_paper_not_washes() {
        let (w, h) = (160, 120);
        let data = (0..w * h).map(|i| sheet_with_wash(i % w, i / w, w, h).0).collect();
        let img = LabImage { width: w, height: h, data };
        let est = surface(&img, 3);
        for y in 0..h {
            for x in 0..w {
                let paper = sheet_with_wash(x, y, w, h).0;
                let (paper, wash) = if sheet_with_wash(x, y, w, h).1 {
                    ([paper[0] + 4.0, paper[1] + 4.0, paper[2] + 2.0], true)
                } else {
                    (paper, false)
                };
                let e = est[y * w + x];
                let tol = if wash { 0.8 } else { 0.5 };
                assert!((0..3).all(|c| (e[c] - paper[c]).abs() < tol), "({x}, {y}) wash={wash}: {e:?} vs {paper:?}");
            }
        }
    }

    #[test]
    fn least_squares_recovers_polynomial() {
        let (w, h, d) = (40, 30, 3);
        let k = 10;
        let mut x = vec![0.0; w * h * k];
        x.chunks_exact_mut(k).enumerate().for_each(|(i, r)| monomials(i % w, i / w, w, h, d, r));
        let truth: Vec<f64> = (0..k).map(|i| i as f64 - 4.5).collect();
        let y: Vec<[f64; 1]> = x.chunks_exact(k).map(|r| [eval(&truth, r)]).collect();
        let [beta] = weighted_lstsq(&x, k, &y, &vec![1.0; w * h]);
        assert!(beta.iter().zip(&truth).all(|(b, t)| (b - t).abs() < 1e-8), "{beta:?}");
    }

    /// Paper lightness following a polynomial of degree `d` with ±1 L* texture noise.
    fn lightness_sheet(d: usize) -> LabImage {
        let (w, h) = (200, 150);
        let data = (0..w * h)
            .map(|i| {
                let (x, y) = ((i % w) as f64 / w as f64 * 2.0 - 1.0, (i / w) as f64 / h as f64 * 2.0 - 1.0);
                let l = match d {
                    1 => 70.0 + 6.0 * x,
                    _ => 70.0 + 6.0 * x - 8.0 * x * x * y + 5.0 * y * y * y,
                };
                let noise = ((i as f64 * 0.754_877_666).fract() - 0.5) * 2.0;
                [l + noise, 1.0, 15.0]
            })
            .collect();
        LabImage { width: w, height: h, data }
    }

    #[test]
    fn degree_selection_finds_the_needed_flexibility() {
        for (d, expected) in [(1, 1), (3, 3)] {
            let img = lightness_sheet(d);
            let paper = PaperModel::estimate(&img, &Hints::default()).unwrap();
            let got = select_degree(&img, &paper.candidates(&img, &Hints::default()));
            assert_eq!(got, expected, "lighting of degree {d}");
        }
    }

    #[test]
    fn scanner_background_is_not_content() {
        // Paper in the centre, white scanner bed with a colour chart patch around it
        let (w, h) = (120, 90);
        let data = (0..w * h)
            .map(|i| {
                let (x, y) = (i % w, i / w);
                if (15..105).contains(&x) && (10..80).contains(&y) {
                    [72.0, 1.0, 16.0]
                } else if x < 10 && y < 10 {
                    [50.0, 60.0, 40.0]
                } else {
                    [97.0, 0.0, 0.0]
                }
            })
            .collect();
        let img = LabImage { width: w, height: h, data };
        let paper = PaperModel::estimate(&img, &Hints::default()).unwrap();
        assert!((paper.global[0] - 72.0).abs() < 1e-9, "{:?}", paper.global);
        assert!(!paper.content[0] && !paper.content[5 * w + 50] && paper.content[45 * w + 60]);
    }

    #[test]
    fn paper_touching_the_border_is_content() {
        // Scan cropped to the map: paper up to the border, cells closed by ink lines
        let (w, h) = (120, 90);
        let data = (0..w * h)
            .map(|i| if (i % w) % 20 == 10 || (i / w) % 20 == 10 { [20.0, 0.0, 2.0] } else { [72.0, 1.0, 16.0] })
            .collect();
        let img = LabImage { width: w, height: h, data };
        let paper = PaperModel::estimate(&img, &Hints::default()).unwrap();
        assert!(paper.content.iter().all(|&c| c));
    }
}
