//! Paper (background) estimation on a low-resolution Lab thumbnail.
//! Mirrors `paper_global` and `paper_poly` of `proto/homog_proto.py`.

use rayon::prelude::*;

use crate::color::Vec3;

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

/// Robust global paper colour: median Lab of the bright, paper-chroma pixels
/// of the central region (`margin` = border fraction ignored on each side).
pub fn paper_global(img: &LabImage, margin: f64) -> Vec3 {
    let (w, h) = (img.width, img.height);
    let (mw, mh) = ((w as f64 * margin) as usize, (h as f64 * margin) as usize);
    let central: Vec<Vec3> = (mh..h - mh)
        .flat_map(|y| (mw..w - mw).map(move |x| y * w + x))
        .map(|i| img.data[i])
        .collect();
    let mut ls: Vec<f64> = central.iter().map(|p| p[0]).collect();
    let lo = percentile(&mut ls, 75.0);
    let hi = percentile(&mut ls, 98.0);
    let sel: Vec<Vec3> = central.into_iter().filter(|p| p[0] >= lo && p[0] <= hi).collect();
    let ab_med = [0.0, median(sel.iter().map(|p| p[1]).collect()), median(sel.iter().map(|p| p[2]).collect())];
    let sel: Vec<Vec3> = sel.into_iter().filter(|p| ab_dist(p, &ab_med) < 6.0).collect();
    median3(&sel)
}

/// Pixels loosely compatible with the global paper colour (excludes ink, strong
/// washes and the scanner background).
fn paper_candidates(img: &LabImage, global: &Vec3) -> Vec<usize> {
    (0..img.data.len())
        .filter(|&i| {
            let p = &img.data[i];
            let dl = p[0] - global[0];
            ab_dist(p, global) < 8.0 && dl > -15.0 && dl < 10.0
        })
        .collect()
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
    /// Not compatible with the paper colour (ink, strong wash, background)
    Other,
    /// Used for the L* envelope only
    Lightness,
    /// Used for the L* envelope and for a*, b* (bare paper)
    Chroma,
}

/// Paper surface: one Lab value per thumbnail pixel, and the role of each pixel in the fit.
pub struct PaperSurface {
    pub lab: Vec<Vec3>,
    pub role: Vec<PixelRole>,
}

/// Paper map as smooth polynomial surfaces of the given degree, robust to large washes.
/// Mirrors `paper_poly` of the prototype.
///
/// Bare paper is the lightest material (a transparent wash or ink can only darken it):
/// - L*: upper-envelope fit, i.e. quantile regression (80 %) by iteratively reweighted
///   least squares, over the pixels loosely compatible with the global paper colour;
/// - a*, b*: least squares with Tukey reweighting (cut-off 3 Δab units, starting from the
///   global paper chroma) over the pixels at most 4 L* units below the envelope, so that
///   washes do not contribute.
///
/// A low-degree surface cannot follow a wash patch, however large.
pub fn paper_surface(img: &LabImage, global: Vec3, degree: usize) -> PaperSurface {
    const TAU: f64 = 0.8;
    const NEAR: f64 = 4.0;
    const C_AB: f64 = 3.0;
    let (w, h) = (img.width, img.height);
    let k = (degree + 1) * (degree + 2) / 2;
    let mut role = vec![PixelRole::Other; w * h];
    let fallback = |role| PaperSurface { lab: vec![global; w * h], role };

    let mut design = vec![0.0; w * h * k];
    design.par_chunks_exact_mut(k).enumerate().for_each(|(i, row)| monomials(i % w, i / w, w, h, degree, row));
    let rows = |idx: &[usize]| -> Vec<f64> { idx.iter().flat_map(|&i| design[i * k..][..k].iter().copied()).collect() };

    let cand = paper_candidates(img, &global);
    cand.iter().for_each(|&i| role[i] = PixelRole::Lightness);
    if cand.len() < 20 * k {
        return fallback(role);
    }
    let xc = rows(&cand);
    let l: Vec<[f64; 1]> = cand.iter().map(|&i| [img.data[i][0]]).collect();
    let mut wt = vec![1.0; cand.len()];
    let mut beta_l = Vec::new();
    for _ in 0..=30 {
        [beta_l] = weighted_lstsq(&xc, k, &l, &wt);
        for ((xi, li), wi) in xc.chunks_exact(k).zip(&l).zip(wt.iter_mut()) {
            let r = li[0] - eval(&beta_l, xi);
            *wi = if r > 0.0 { TAU } else { 1.0 - TAU } / r.abs().max(0.05);
        }
    }
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
    PaperSurface { lab, role }
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
        let g = paper_global(&img, 0.05);
        assert!((0..3).all(|c| (g[c] - paper[c]).abs() < 1e-9), "{g:?}");
        let est = paper_surface(&img, g, 3).lab;
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
        let est = paper_surface(&img, paper_global(&img, 0.05), 3).lab;
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
}
