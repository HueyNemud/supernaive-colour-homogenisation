# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "pillow"]
# ///
"""
Prototype for the colour homogenisation of scanned map sheets.

Reference implementation of the model ported to Rust (`homog`), and test bed for
its evolutions. Variants:

  bradford-global     Bradford chromatic adaptation of a single paper colour to the target
  bradford-surface    same, with the paper colour given by a robust polynomial surface
                      (flat-field: corrects uneven lighting, keeps the washes) -- the `homog` model
  bradford-seams      bradford-surface plus a seam correction: along the edges shared with
                      adjacent sheets (layout given by their GCPs), the paper is measured in a
                      thin band and brought back to the target, the correction fading inwards

Usage:
  uv run proto/homog_proto.py example/*.jpg [-g example/GCP] [-s 0.25] [-o proto/out]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

# ---------------------------------------------------------------------------
# Colour conversions (sRGB D65). Kept explicit so the Rust port is 1:1.
# ---------------------------------------------------------------------------

M_RGB2XYZ = np.array([
    [0.4124564, 0.3575761, 0.1804375],
    [0.2126729, 0.7151522, 0.0721750],
    [0.0193339, 0.1191920, 0.9503041],
])
M_XYZ2RGB = np.linalg.inv(M_RGB2XYZ)
WHITE_D65 = M_RGB2XYZ @ np.ones(3)

M_BRADFORD = np.array([
    [0.8951, 0.2664, -0.1614],
    [-0.7502, 1.7135, 0.0367],
    [0.0389, -0.0685, 1.0296],
])
M_BRADFORD_INV = np.linalg.inv(M_BRADFORD)

DELTA = 6 / 29


def srgb_to_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def srgb_to_xyz(rgb):
    return srgb_to_linear(rgb) @ M_RGB2XYZ.T


def xyz_to_srgb(xyz):
    return linear_to_srgb(xyz @ M_XYZ2RGB.T)


def _f(t):
    return np.where(t > DELTA ** 3, np.cbrt(t), t / (3 * DELTA ** 2) + 4 / 29)


def _finv(t):
    return np.where(t > DELTA, t ** 3, 3 * DELTA ** 2 * (t - 4 / 29))


def xyz_to_lab(xyz):
    f = _f(xyz / WHITE_D65)
    L = 116 * f[..., 1] - 16
    a = 500 * (f[..., 0] - f[..., 1])
    b = 200 * (f[..., 1] - f[..., 2])
    return np.stack([L, a, b], axis=-1)


def lab_to_xyz(lab):
    fy = (lab[..., 0] + 16) / 116
    fx = fy + lab[..., 1] / 500
    fz = fy - lab[..., 2] / 200
    return np.stack([_finv(fx), _finv(fy), _finv(fz)], axis=-1) * WHITE_D65


def srgb_to_lab(rgb):
    return xyz_to_lab(srgb_to_xyz(rgb))


def lab_to_srgb(lab):
    return xyz_to_srgb(lab_to_xyz(lab))


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def load(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("RGB"), dtype=np.float64) / 255.0


def save(rgb: np.ndarray, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(np.round(np.clip(rgb, 0, 1) * 255).astype(np.uint8)).save(path)


def thumbnail(rgb: np.ndarray, max_side: int = 400) -> np.ndarray:
    """Box downscale done in linear light, returned as sRGB."""
    h, w, _ = rgb.shape
    s = max(1, int(np.ceil(max(h, w) / max_side)))
    lin = srgb_to_linear(rgb)
    hh, ww = h // s * s, w // s * s
    lin = lin[:hh, :ww].reshape(hh // s, s, ww // s, s, 3).mean(axis=(1, 3))
    return linear_to_srgb(lin)


def upsample(img: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Bilinear upsampling of a low-res (h, w, c) map to `shape`, pixel-centre aligned."""
    h, w = img.shape[:2]
    H, W = shape
    yy = (np.arange(H) + 0.5) * h / H - 0.5
    xx = (np.arange(W) + 0.5) * w / W - 0.5
    gy, gx = np.meshgrid(yy, xx, indexing="ij")
    return np.stack(
        [ndimage.map_coordinates(img[..., c], [gy, gx], order=1, mode="nearest")
         for c in range(img.shape[2])],
        axis=-1,
    )


# ---------------------------------------------------------------------------
# Background (paper) estimation
# ---------------------------------------------------------------------------


def paper_global(lab_thumb: np.ndarray, margin: float) -> np.ndarray:
    """Robust global paper colour: median Lab of the bright, paper-chroma pixels
    of the central region (margins excluded)."""
    h, w, _ = lab_thumb.shape
    mh, mw = int(h * margin), int(w * margin)
    c = lab_thumb[mh:h - mh, mw:w - mw].reshape(-1, 3)
    L = c[:, 0]
    lo, hi = np.percentile(L, [75, 98])
    sel = c[(L >= lo) & (L <= hi)]
    # Keep the dominant chroma among bright pixels (rejects pale washes)
    ab_med = np.median(sel[:, 1:], axis=0)
    sel = sel[np.linalg.norm(sel[:, 1:] - ab_med, axis=1) < 6]
    return np.median(sel, axis=0)


def paper_mask(lab: np.ndarray, ref: np.ndarray, d_ab=8.0, d_dark=15.0, d_light=10.0):
    """Pixels loosely compatible with the paper colour `ref` (excludes ink, strong washes, scanner background)."""
    dab = np.linalg.norm(lab[..., 1:] - ref[..., 1:], axis=-1)
    dL = lab[..., 0] - ref[..., 0]
    return (dab < d_ab) & (dL > -d_dark) & (dL < d_light)


def poly_design(h, w, degree):
    """Monomials x^i y^j (i + j <= degree) on coordinates normalised to [-1, 1]."""
    y, x = np.mgrid[0:h, 0:w]
    x = (x + 0.5) / w * 2 - 1
    y = (y + 0.5) / h * 2 - 1
    terms = [x ** i * y ** j for i in range(degree + 1) for j in range(degree + 1 - i)]
    return np.stack([t.ravel() for t in terms], axis=1)


def weighted_lstsq(X, y, w):
    sw = np.sqrt(w)
    return np.linalg.lstsq(X * sw[:, None], y * sw[..., None] if y.ndim == 2 else y * sw, rcond=None)[0]


def paper_poly(lab_thumb, global_paper, degree=3, tau=0.8, near=4.0, c_ab=3.0, iterations=30):
    """Paper map as smooth polynomial surfaces, robust to large washes.

    Bare paper is the lightest material (a transparent wash or ink can only darken it):
    - L*: upper-envelope fit, i.e. quantile regression (quantile `tau`) by IRLS, over the
      pixels loosely compatible with the global paper colour;
    - a*, b*: least squares with Tukey reweighting (cut-off `c_ab` Δab units, starting from
      the global paper chroma) over the pixels within `near` L* units below the envelope,
      so that washes do not contribute.
    A low-degree surface cannot follow a wash patch, however large."""
    h, w, _ = lab_thumb.shape
    X = poly_design(h, w, degree)
    Y = lab_thumb.reshape(-1, 3)
    w0 = paper_mask(lab_thumb, np.broadcast_to(global_paper, lab_thumb.shape)).ravel()
    Xs, L = X[w0], Y[w0, 0]

    # Quantile regression: asymmetric absolute loss, IRLS weights tau/|r| or (1-tau)/|r|
    beta_L = weighted_lstsq(Xs, L, np.ones(len(L)))
    for _ in range(iterations):
        r = L - Xs @ beta_L
        beta_L = weighted_lstsq(Xs, L, np.where(r > 0, tau, 1 - tau) / np.maximum(np.abs(r), 0.05))
    L_surf = X @ beta_L

    sel = w0 & (Y[:, 0] - L_surf > -near)
    Xc, AB = X[sel], Y[sel, 1:]
    pred = np.broadcast_to(global_paper[1:], AB.shape)
    for _ in range(10):
        u = np.linalg.norm(AB - pred, axis=1) / c_ab
        wt = np.where(u < 1, (1 - u ** 2) ** 2, 0.0)
        beta_ab = weighted_lstsq(Xc, AB, wt)
        pred = Xc @ beta_ab
    L_surf = np.clip(L_surf, max(global_paper[0] - 30, 5), 100)  # sane extrapolation
    surf = np.column_stack([L_surf, X @ beta_ab]).reshape(h, w, 3)
    mask = np.zeros(h * w, bool)
    mask[np.flatnonzero(sel)[wt > 0]] = True
    return surf, mask.reshape(h, w)


# ---------------------------------------------------------------------------
# Colour models
# ---------------------------------------------------------------------------


def model_bradford(rgb, paper_xyz, target_xyz):
    """von Kries adaptation in Bradford cone space. paper_xyz: (3,) or (H, W, 3)."""
    gain = (M_BRADFORD @ target_xyz) / (paper_xyz @ M_BRADFORD.T)
    lms = srgb_to_xyz(rgb) @ M_BRADFORD.T
    return xyz_to_srgb((lms * gain) @ M_BRADFORD_INV.T)


# ---------------------------------------------------------------------------
# Inter-sheet harmonisation along the seams of adjacent sheets (GCP layout)
# ---------------------------------------------------------------------------


def read_points(path: Path) -> np.ndarray:
    """Enabled GCPs of a QGIS georeferencer `.points` file, as rows (pixel, line, X, Y)."""
    rows = []
    for line in path.read_text().splitlines():
        v = line.strip().split(",")
        if line.startswith(("#", "mapX")) or len(v) < 5 or v[4].strip() != "1":
            continue
        rows.append([float(v[2]), -float(v[3]), float(v[0]), float(v[1])])
    return np.array(rows)


class Layout:
    """Position of a sheet in map space: map -> pixel affine transform (least squares on the
    GCPs) and map extent of its content, taken as the bounding box of the GCPs (which sit
    on the neatline)."""

    def __init__(self, gcps: np.ndarray):
        m = np.column_stack([gcps[:, 2:], np.ones(len(gcps))])
        px = np.column_stack([gcps[:, :2], np.ones(len(gcps))])
        self.to_pixel = np.linalg.lstsq(m, gcps[:, :2], rcond=None)[0]
        self.to_map = np.linalg.lstsq(px, gcps[:, 2:], rcond=None)[0]
        self.rect = (gcps[:, 2].min(), gcps[:, 2].max(), gcps[:, 3].min(), gcps[:, 3].max())

    def pixel(self, xy: np.ndarray) -> np.ndarray:
        return np.column_stack([xy, np.ones(len(xy))]) @ self.to_pixel

    def map(self, px: np.ndarray) -> np.ndarray:
        return np.column_stack([px, np.ones(len(px))]) @ self.to_map


def shared_edges(ri, rj, tol=0.01):
    """Edges of extent ri shared with extent rj, as (axis, position, inward sign, (u0, u1)):
    axis 0 = vertical edge x = position, axis 1 = horizontal edge y = position."""
    size = max(ri[1] - ri[0], ri[3] - ri[2])
    tol = tol * size
    edges = []
    oy0, oy1 = max(ri[2], rj[2]), min(ri[3], rj[3])
    ox0, ox1 = max(ri[0], rj[0]), min(ri[1], rj[1])
    if oy1 - oy0 > tol:
        if abs(ri[1] - rj[0]) < tol:
            edges.append((0, ri[1], -1, (oy0, oy1)))
        if abs(ri[0] - rj[1]) < tol:
            edges.append((0, ri[0], 1, (oy0, oy1)))
    if ox1 - ox0 > tol:
        if abs(ri[3] - rj[2]) < tol:
            edges.append((1, ri[3], -1, (ox0, ox1)))
        if abs(ri[2] - rj[3]) < tol:
            edges.append((1, ri[2], 1, (ox0, ox1)))
    return edges


def edge_gain_field(lin_rgb, layout, edges, scale, target_xyz, d0=0.003, d1=0.013, reach=0.12, n_bins=24):
    """Log Bradford-LMS gain field (H, W, 3) bringing the paper back to the target along the
    given seam edges: paper is measured in a band [d0, d1] inside each edge, per segment, and
    the correction fades out (squared linear ramp) over `reach` (fractions of the sheet size)."""
    h, w, _ = lin_rgb.shape
    r = layout.rect
    size = max(r[1] - r[0], r[3] - r[2])
    d0, d1, reach = d0 * size, d1 * size, reach * size
    lms_t = M_BRADFORD @ target_xyz
    yy, xx = np.mgrid[0:h, 0:w]
    xy = layout.map(np.column_stack([(xx.ravel() + 0.5) * scale, (yy.ravel() + 0.5) * scale]))
    field = np.zeros((h * w, 3))
    weight = np.zeros(h * w)
    for axis, pos, sign, (u0, u1) in edges:
        u = np.linspace(u0, u1, 400)
        t = np.linspace(d0, d1, 8)
        U, T = np.meshgrid(u, t)
        pts = np.column_stack([pos + sign * T.ravel(), U.ravel()]) if axis == 0 else np.column_stack([U.ravel(), pos + sign * T.ravel()])
        vals, ok = sample(lin_rgb, layout.pixel(pts), scale)
        lab = xyz_to_lab(vals @ M_RGB2XYZ.T)
        paper = ok & (np.linalg.norm(lab[:, 1:], axis=1) < 4) & (lab[:, 0] > 80)
        lms = vals @ M_RGB2XYZ.T @ M_BRADFORD.T
        bins = np.minimum(((U.ravel() - u0) / (u1 - u0) * n_bins).astype(int), n_bins - 1)
        gain = np.full((n_bins, 3), np.nan)
        for k in range(n_bins):
            sel = paper & (bins == k)
            if sel.sum() >= 10:
                gain[k] = np.log(lms_t / np.median(lms[sel], axis=0))
        known = ~np.isnan(gain[:, 0])
        if known.sum() == 0:
            continue
        centres = (np.arange(n_bins) + 0.5) / n_bins
        gain = np.column_stack([np.interp(centres, centres[known], gain[known, c]) for c in range(3)])
        gain = ndimage.uniform_filter1d(gain, 3, axis=0, mode="nearest")
        # Spread inside the sheet
        d = sign * ((xy[:, 0] if axis == 0 else xy[:, 1]) - pos)
        uu = np.clip(((xy[:, 1] if axis == 0 else xy[:, 0]) - u0) / (u1 - u0), 0, 1)
        phi = np.clip(1 - d / reach, 0, 1) ** 2
        g = np.column_stack([np.interp(uu, centres, gain[:, c]) for c in range(3)])
        field += phi[:, None] * g
        weight += phi
    # Near a corner both edges measure the same darkening: average rather than add
    return (field / np.maximum(weight, 1)[:, None]).reshape(h, w, 3)


def apply_log_gain(lin_rgb, field):
    lms = lin_rgb @ M_RGB2XYZ.T @ M_BRADFORD.T
    return (lms * np.exp(field)) @ M_BRADFORD_INV.T @ np.linalg.inv(M_RGB2XYZ).T


def sample(lin_rgb: np.ndarray, pixels: np.ndarray, scale: float):
    """Bilinear sampling of a (reduced by `scale`) image at full-resolution pixel coordinates.
    Returns values and a validity mask (inside the image)."""
    h, w, _ = lin_rgb.shape
    x, y = pixels[:, 0] / scale - 0.5, pixels[:, 1] / scale - 0.5
    ok = (x >= 0) & (y >= 0) & (x <= w - 1) & (y <= h - 1)
    vals = np.stack([ndimage.map_coordinates(lin_rgb[..., c], [y, x], order=1, mode="nearest") for c in range(3)], -1)
    return vals, ok


def seam_points(ri, rj, d0=0.003, d1=0.013, n_along=400, n_across=8):
    """Map points mirrored across the edge shared by extents ri and rj, between d0 and d1
    from it (fractions of the sheet size): (points in i, points in j), or None."""
    edges = shared_edges(ri, rj)
    if not edges:
        return None
    axis, pos, sign, (u0, u1) = edges[0]
    size = max(ri[1] - ri[0], ri[3] - ri[2])
    U, T = np.meshgrid(np.linspace(u0, u1, n_along), np.linspace(d0 * size, d1 * size, n_across))
    U, T = U.ravel(), T.ravel()
    side = lambda sgn: np.column_stack([pos + sgn * T, U]) if axis == 0 else np.column_stack([U, pos + sgn * T])
    return side(sign), side(-sign)


def seam_report(pairs, names, label_):
    """Per seam: median Lab on each side of paper pairs and of light wash pairs."""
    for i, j, a, b in pairs:
        la, lb = xyz_to_lab(a @ M_RGB2XYZ.T), xyz_to_lab(b @ M_RGB2XYZ.T)
        same = np.linalg.norm(la - lb, axis=1) < 10
        paper = same & (np.linalg.norm(la[:, 1:], axis=1) < 4) & (np.linalg.norm(lb[:, 1:], axis=1) < 4) & (la[:, 0] > 80) & (lb[:, 0] > 80)
        wash = same & ~paper & (la[:, 0] > 60) & (lb[:, 0] > 60)
        line = f"  {label_:<10}{names[i][-4:]}|{names[j][-4:]}"
        for cls, m in [("paper", paper), ("washes", wash)]:
            ma, mb = np.median(la[m], 0), np.median(lb[m], 0)
            line += f"   {cls} {np.round(ma, 1).tolist()} | {np.round(mb, 1).tolist()} (ΔE {np.linalg.norm(ma - mb):.1f})"
        print(line)


def seam_crops(images, layouts, names, scales, width=0.08, height=0.25):
    """Side-by-side crops (pixel space, no resampling) on both sides of each seam."""
    tiles = []
    for i in range(len(images)):
        for j in range(i + 1, len(images)):
            if layouts[i] is None or layouts[j] is None:
                continue
            ri, rj = layouts[i].rect, layouts[j].rect
            pts = seam_points(ri, rj)
            if pts is None or abs(ri[1] - rj[0]) > 0.01 * (ri[1] - ri[0]):
                continue  # only west -> east seams are shown
            x, ym = (ri[1] + rj[0]) / 2, (max(ri[2], rj[2]) + min(ri[3], rj[3])) / 2
            dx, dy = width * (ri[1] - ri[0]), height * (ri[3] - ri[2])
            crops = []
            for k, (x0, x1) in [(i, (x - dx, x)), (j, (x, x + dx))]:
                px = layouts[k].pixel(np.array([[x0, ym + dy / 2], [x1, ym - dy / 2]])) / scales[k]
                (c0, r0), (c1, r1) = np.floor(px.min(0)).astype(int), np.ceil(px.max(0)).astype(int)
                crops.append(images[k][max(r0, 0):r1, max(c0, 0):c1])
            hh = min(c.shape[0] for c in crops)
            tiles.append((f"{names[i][-4:]} | {names[j][-4:]}", np.concatenate([c[:hh] for c in crops], 1)))
    return tiles


# ---------------------------------------------------------------------------
# Metrics and montage
# ---------------------------------------------------------------------------


def metrics(out_rgb, orig_lab_thumb, global_paper, global_out_lab, n=8, margin=0.06):
    """Paper whiteness and wash preservation of one rendering (on its thumbnail).

    - paper: L* 90th percentile of the paper-hued pixels (Δab < 3 from the global paper) of each
      block of an n x n grid (margins excluded); min and 10th percentile over the blocks;
    - washes: mean chroma of the light, non-paper-hued pixels relative to `global_out_lab` (the
      bradford-global rendering, which cannot bleach washes)."""
    lab = srgb_to_lab(thumbnail(out_rgb))
    dab = np.linalg.norm(orig_lab_thumb[..., 1:] - global_paper[1:], axis=-1)
    paper_hue = dab < 3
    wash = (dab > 3) & (dab < 15) & (orig_lab_thumb[..., 0] > global_paper[0] - 15)
    h, w = paper_hue.shape
    edges_y = np.linspace(h * margin, h * (1 - margin), n + 1).astype(int)
    edges_x = np.linspace(w * margin, w * (1 - margin), n + 1).astype(int)
    blocks = []
    for y0, y1 in zip(edges_y, edges_y[1:]):
        for x0, x1 in zip(edges_x, edges_x[1:]):
            sel = paper_hue[y0:y1, x0:x1]
            if sel.sum() > 50:
                blocks.append(np.percentile(lab[y0:y1, x0:x1][sel][:, 0], 90))
    chroma = lambda l: np.linalg.norm(l[wash][:, 1:], axis=1).mean()
    return np.min(blocks), np.percentile(blocks, 10), chroma(lab) / chroma(global_out_lab)


def label(img: Image.Image, text: str) -> Image.Image:
    bar = Image.new("RGB", (img.width, 28), "white")
    ImageDraw.Draw(bar).text((8, 6), text, fill="black")
    out = Image.new("RGB", (img.width, img.height + 28), "white")
    out.paste(bar, (0, 0))
    out.paste(img, (0, 28))
    return out


def montage(rows: dict[str, list[np.ndarray]], path: Path, height=360):
    imgs = []
    for name, sheets in rows.items():
        tiles = []
        for s in sheets:
            im = Image.fromarray(np.round(np.clip(s, 0, 1) * 255).astype(np.uint8))
            tiles.append(im.resize((round(im.width * height / im.height), height), Image.LANCZOS))
        row = Image.new("RGB", (sum(t.width for t in tiles) + 8 * (len(tiles) - 1), height), "white")
        x = 0
        for t in tiles:
            row.paste(t, (x, 0))
            x += t.width + 8
        imgs.append(label(row, name))
    W = max(i.width for i in imgs)
    out = Image.new("RGB", (W, sum(i.height for i in imgs)), "white")
    y = 0
    for i in imgs:
        out.paste(i, (0, y))
        y += i.height
    out.save(path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("-o", "--out", type=Path, default=Path("proto/out"))
    ap.add_argument("-t", "--target", default="255,255,255")
    ap.add_argument("-m", "--margin", type=float, default=0.05, help="border fraction ignored for the global paper estimate")
    ap.add_argument("-d", "--degree", type=int, default=3, help="degree of the paper surface")
    ap.add_argument("-g", "--gcp-dir", type=Path, help="directory of QGIS .points files (<image name>.points) giving the sheet layout")
    ap.add_argument("-s", "--scale", type=float, default=0.25, help="processing scale (1 = full resolution)")
    args = ap.parse_args()

    target_rgb = np.array([float(v) for v in args.target.split(",")]) / 255
    target_xyz = srgb_to_xyz(target_rgb)
    target_lab = xyz_to_lab(target_xyz)

    names = [f.stem for f in args.files]
    images = [load(f) for f in args.files]
    widths = [im.shape[1] for im in images]
    if args.scale < 1:
        images = [thumbnail(im, round(max(im.shape[:2]) * args.scale)) for im in images]
    scales = [w / im.shape[1] for w, im in zip(widths, images)]  # full-resolution px per processed px
    thumbs_lab = [srgb_to_lab(thumbnail(im)) for im in images]
    papers, surfaces = [], []
    for name, lab in zip(names, thumbs_lab):
        g = paper_global(lab, args.margin)
        surf, _ = paper_poly(lab, g, degree=args.degree)
        papers.append(g)
        surfaces.append(surf)
        print(f"{name}: global paper Lab = {np.round(g, 2)}, "
              f"surface L* range = [{surf[..., 0].min():.1f}, {surf[..., 0].max():.1f}]")

    results: dict[str, list[np.ndarray]] = {"original": images}
    results["bradford-global"] = [model_bradford(im, lab_to_xyz(g), target_xyz) for im, g in zip(images, papers)]
    results["bradford-surface"] = [
        model_bradford(im, lab_to_xyz(upsample(s, im.shape[:2])), target_xyz) for im, s in zip(images, surfaces)
    ]

    # Inter-sheet harmonisation along the seams (sheets laid out by their GCPs)
    layouts = []
    for f in args.files:
        pts = args.gcp_dir / f"{f.name}.points" if args.gcp_dir else None
        layouts.append(Layout(read_points(pts)) if pts and pts.exists() else None)
    lin = [srgb_to_linear(o) for o in results["bradford-surface"]]

    def seam_pairs(lin):
        pairs = []
        for i in range(len(images)):
            for j in range(i + 1, len(images)):
                if layouts[i] is None or layouts[j] is None:
                    continue
                sp = seam_points(layouts[i].rect, layouts[j].rect)
                if sp is None:
                    continue
                a, ok_a = sample(lin[i], layouts[i].pixel(sp[0]), scales[i])
                b, ok_b = sample(lin[j], layouts[j].pixel(sp[1]), scales[j])
                pairs.append((i, j, a[ok_a & ok_b], b[ok_a & ok_b]))
        return pairs

    pairs = seam_pairs(lin)
    if pairs:
        print("seams (medians per side):")
        seam_report(pairs, names, "surface")
    if pairs:
        all_rects = [Layout(read_points(f)).rect for f in sorted(args.gcp_dir.glob("*.points"))]
        for i in range(len(images)):
            if layouts[i] is None:
                continue
            edges = [e for r in all_rects if r != layouts[i].rect for e in shared_edges(layouts[i].rect, r)]
            lin[i] = apply_log_gain(lin[i], edge_gain_field(lin[i], layouts[i], edges, scales[i], target_xyz))
        results["bradford-seams"] = [linear_to_srgb(l) for l in lin]
        seam_report(seam_pairs(lin), names, "seams")

    # Outputs
    args.out.mkdir(parents=True, exist_ok=True)
    for variant, outs in results.items():
        if variant != "original":
            for name, o in zip(names, outs):
                save(o, args.out / variant / f"{name}.png")
    montage(results, args.out / "comparison.png")
    if pairs:
        rows = {}
        for variant in ["original", "bradford-surface", "bradford-seams"]:
            for label_, crop in seam_crops(results[variant], layouts, names, scales):
                rows.setdefault(label_, []).append(crop)
        montage({f"{k}  (original / surface / seams)": v for k, v in rows.items()}, args.out / "seams.png", height=500)

    # Metrics
    print(f"\n{'variant':<20}{'sheet':<12}{'paper L* min / p10':>20}{'wash chroma':>13}")
    globals_lab = [srgb_to_lab(thumbnail(o)) for o in results["bradford-global"]]
    for variant, outs in results.items():
        for name, o, tl, g, gl in zip(names, outs, thumbs_lab, papers, globals_lab):
            mn, p10, wc = metrics(o, tl, g, gl)
            print(f"{variant:<20}{name[-10:]:<12}{mn:>12.1f} / {p10:.1f}{wc:>13.2f}")
    print(f"\nOutputs in {args.out}/ (comparison.png = montage, seams.png = seams)")


if __name__ == "__main__":
    main()
