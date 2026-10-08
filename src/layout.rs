//! Layout of the sheets in map space, from their ground control points (GCPs).
//! No resampling is involved: the layout only tells which sheets touch and where.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{ensure, Context, Result};

use crate::background::{eval, weighted_lstsq};

/// Enabled GCPs of a QGIS georeferencer `.points` file, as (pixel, line, X, Y).
pub fn read_points(path: &Path) -> Result<Vec<[f64; 4]>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let mut gcps = Vec::new();
    for line in text.lines() {
        let v: Vec<&str> = line.trim().split(',').collect();
        if line.starts_with('#') || line.starts_with("mapX") || v.len() < 5 || v[4].trim() != "1" {
            continue;
        }
        let num = |i: usize| -> Result<f64> {
            v[i].trim().parse().with_context(|| format!("{}: invalid number '{}'", path.display(), v[i]))
        };
        // QGIS stores image rows as negative sourceY
        gcps.push([num(2)?, -num(3)?, num(0)?, num(1)?]);
    }
    Ok(gcps)
}

/// Map extent (x0, x1, y0, y1).
pub type Rect = [f64; 4];

/// Position of a sheet in map space: affine transforms map <-> pixel (least squares on the
/// GCPs) and extent of its content, taken as the bounding box of the GCPs (which sit on the
/// neatline).
#[derive(Clone, Debug)]
pub struct Layout {
    to_pixel: Affine,
    to_map: Affine,
    pub rect: Rect,
    /// RMS residual of the GCPs, in map units: how precisely the neatline is located
    pub rms: f64,
}

/// 2-D affine transform fitted by least squares on centred coordinates (well conditioned).
#[derive(Clone, Debug)]
struct Affine {
    centre: [f64; 2],
    coefs: [Vec<f64>; 2],
}

impl Affine {
    fn fit(src: &[[f64; 2]], dst: &[[f64; 2]]) -> Self {
        let n = src.len() as f64;
        let centre = [0, 1].map(|c| src.iter().map(|p| p[c]).sum::<f64>() / n);
        let x: Vec<f64> = src.iter().flat_map(|p| [p[0] - centre[0], p[1] - centre[1], 1.0]).collect();
        Self { centre, coefs: weighted_lstsq(&x, 3, dst, &vec![1.0; src.len()]) }
    }

    fn apply(&self, p: [f64; 2]) -> [f64; 2] {
        let v = [p[0] - self.centre[0], p[1] - self.centre[1], 1.0];
        [eval(&self.coefs[0], &v), eval(&self.coefs[1], &v)]
    }
}

impl Layout {
    pub fn from_gcps(gcps: &[[f64; 4]]) -> Result<Self> {
        ensure!(gcps.len() >= 3, "at least 3 enabled GCPs are needed, found {}", gcps.len());
        let pixels: Vec<[f64; 2]> = gcps.iter().map(|g| [g[0], g[1]]).collect();
        let points: Vec<[f64; 2]> = gcps.iter().map(|g| [g[2], g[3]]).collect();
        let (xs, ys) = (gcps.iter().map(|g| g[2]), gcps.iter().map(|g| g[3]));
        let rect = [
            xs.clone().fold(f64::INFINITY, f64::min),
            xs.fold(f64::NEG_INFINITY, f64::max),
            ys.clone().fold(f64::INFINITY, f64::min),
            ys.fold(f64::NEG_INFINITY, f64::max),
        ];
        let to_map = Affine::fit(&pixels, &points);
        let rms = (pixels.iter().zip(&points).map(|(p, q)| {
            let m = to_map.apply(*p);
            (m[0] - q[0]).powi(2) + (m[1] - q[1]).powi(2)
        }).sum::<f64>() / gcps.len() as f64).sqrt();
        Ok(Self { to_pixel: Affine::fit(&points, &pixels), to_map, rect, rms })
    }

    /// Map coordinates -> full-resolution pixel coordinates (continuous, 0 = left edge).
    pub fn pixel(&self, xy: [f64; 2]) -> [f64; 2] {
        self.to_pixel.apply(xy)
    }

    /// Full-resolution pixel coordinates -> map coordinates.
    pub fn map(&self, px: [f64; 2]) -> [f64; 2] {
        self.to_map.apply(px)
    }

    /// Largest side of the extent, the unit of the seam distances.
    pub fn size(&self) -> f64 {
        (self.rect[1] - self.rect[0]).max(self.rect[3] - self.rect[2])
    }
}

/// An edge of a sheet's extent shared with another sheet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Edge {
    /// 0: vertical edge x = `pos`; 1: horizontal edge y = `pos`
    pub axis: usize,
    pub pos: f64,
    /// +1 if the sheet lies on the increasing side of the edge, -1 otherwise
    pub inward: f64,
    /// Extent of the shared part along the edge
    pub span: (f64, f64),
}

/// Edges of extent `ri` shared with extent `rj` (tolerance: 1 % of the sheet size).
pub fn shared_edges(ri: &Rect, rj: &Rect) -> Vec<Edge> {
    let tol = 0.01 * (ri[1] - ri[0]).max(ri[3] - ri[2]);
    let (ox0, ox1) = (ri[0].max(rj[0]), ri[1].min(rj[1]));
    let (oy0, oy1) = (ri[2].max(rj[2]), ri[3].min(rj[3]));
    let mut edges = Vec::new();
    if oy1 - oy0 > tol {
        if (ri[1] - rj[0]).abs() < tol {
            edges.push(Edge { axis: 0, pos: ri[1], inward: -1.0, span: (oy0, oy1) });
        }
        if (ri[0] - rj[1]).abs() < tol {
            edges.push(Edge { axis: 0, pos: ri[0], inward: 1.0, span: (oy0, oy1) });
        }
    }
    if ox1 - ox0 > tol {
        if (ri[3] - rj[2]).abs() < tol {
            edges.push(Edge { axis: 1, pos: ri[3], inward: -1.0, span: (ox0, ox1) });
        }
        if (ri[2] - rj[3]).abs() < tol {
            edges.push(Edge { axis: 1, pos: ri[2], inward: 1.0, span: (ox0, ox1) });
        }
    }
    edges
}

/// Layout of a sheet and the edges it shares with adjacent sheets.
pub type Seams = (Layout, Vec<Edge>);

/// Layout of a whole atlas: one `<image file name>.points` file per sheet in a directory.
pub struct Atlas {
    sheets: HashMap<String, Layout>,
}

impl Atlas {
    pub fn load(dir: &Path) -> Result<Self> {
        let mut sheets = HashMap::new();
        for entry in std::fs::read_dir(dir).with_context(|| format!("cannot read {}", dir.display()))? {
            let path = entry?.path();
            if let Some(name) = path.file_name().and_then(|n| n.to_str()).and_then(|n| n.strip_suffix(".points")) {
                let layout = Layout::from_gcps(&read_points(&path)?).with_context(|| format!("{}", path.display()))?;
                sheets.insert(name.to_owned(), layout);
            }
        }
        Ok(Self { sheets })
    }

    pub fn is_empty(&self) -> bool {
        self.sheets.is_empty()
    }

    /// Layout of the sheet stored as `file_name` and its edges shared with any other sheet of
    /// the atlas (whether or not it is processed in the same run).
    pub fn seams(&self, file_name: &str) -> Option<Seams> {
        let layout = self.sheets.get(file_name)?;
        let edges = self
            .sheets
            .iter()
            .filter(|(name, _)| name.as_str() != file_name)
            .flat_map(|(_, other)| shared_edges(&layout.rect, &other.rect))
            .collect();
        Some((layout.clone(), edges))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // pixel = 8.75 (X + 1600) + 286, line = 8.9 (2600 - Y) + 208, plus a disabled GCP
    const POINTS: &str = "#CRS: LOCAL\nmapX,mapY,sourceX,sourceY,enable,dX,dY,residual\n\
        -1600.0,2600.0,286.0,-208.0,1,0,0,0\n-1000.0,2600.0,5536.0,-208.0,1,0,0,0\n\
        -1000.0,2200.0,5536.0,-3768.0,1,0,0,0\n-1600.0,2200.0,286.0,-3768.0,1,0,0,0\n\
        -1300.0,2400.0,1.0,-1.0,0,0,0,0\n";

    #[test]
    fn parses_qgis_points() {
        let path = std::env::temp_dir().join(format!("homog-points-{}.points", std::process::id()));
        std::fs::write(&path, POINTS).unwrap();
        let gcps = read_points(&path).unwrap();
        assert_eq!(gcps.len(), 4, "disabled GCP skipped");
        assert_eq!(gcps[0], [286.0, 208.0, -1600.0, 2600.0]);
        let layout = Layout::from_gcps(&gcps).unwrap();
        assert_eq!(layout.rect, [-1600.0, -1000.0, 2200.0, 2600.0]);
        assert!(layout.rms < 1e-6, "exactly affine GCPs: rms {}", layout.rms);
        for g in &gcps {
            let px = layout.pixel([g[2], g[3]]);
            assert!((px[0] - g[0]).abs() < 1e-6 && (px[1] - g[1]).abs() < 1e-6, "{px:?} vs {g:?}");
            let xy = layout.map(px);
            assert!((xy[0] - g[2]).abs() < 1e-6 && (xy[1] - g[3]).abs() < 1e-6, "round trip {xy:?}");
        }
    }

    #[test]
    fn finds_shared_edges() {
        let a = [-1600.0, -1000.0, 2200.0, 2600.0];
        let b = [-1000.0, -400.0, 2200.0, 2600.0];
        let c = [-1600.0, -1000.0, 2600.0, 3000.0];
        let far = [0.0, 600.0, 0.0, 400.0];
        assert_eq!(shared_edges(&a, &b), vec![Edge { axis: 0, pos: -1000.0, inward: -1.0, span: (2200.0, 2600.0) }]);
        assert_eq!(shared_edges(&b, &a), vec![Edge { axis: 0, pos: -1000.0, inward: 1.0, span: (2200.0, 2600.0) }]);
        assert_eq!(shared_edges(&a, &c), vec![Edge { axis: 1, pos: 2600.0, inward: -1.0, span: (-1600.0, -1000.0) }]);
        assert!(shared_edges(&a, &far).is_empty());
        // Corner contact only: no shared edge
        assert!(shared_edges(&b, &c).is_empty());
    }
}
