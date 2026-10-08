//! Background colour homogenisation of scanned map sheets.
//!
//! Pipeline per image:
//! 1. downsampled read -> Lab thumbnail;
//! 2. paper estimation: robust global colour, then a smooth paper surface (flat-field);
//! 3. seams (optional, layout from GCPs): along the edges shared with adjacent sheets, the
//!    remaining paper colour is measured on a reduced image and corrected, fading inwards;
//! 4. one streaming pass over full-resolution strips: decode (sRGB LUT or ICC),
//!    Bradford adaptation of the local paper colour to the target white, sRGB encode.

pub mod background;
pub mod color;
pub mod debug;
pub mod io;
pub mod layout;
pub mod seams;
pub mod transform;

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use lcms2::{CIExyY, CIExyYTRIPLE, DisallowCache, Flags, GlobalContext, Intent, PixelFormat, Profile, ToneCurve};
use rayon::prelude::*;

use background::{paper_global, paper_surface, LabImage, PaperSurface, PixelRole};
use color::{decode_table, linear_rgb_to_lab, srgb_to_linear, srgb_to_xyz_u8, SrgbEncoder, Vec3};
use io::{Input, Output, OutputOptions, Sample, SampleType};
use layout::{Edge, Layout};
use seams::EdgeField;
use transform::{Adapter, GainMap};

pub struct Options {
    /// Target background colour (sRGB, 8-bit)
    pub target: [u8; 3],
    /// Border fraction ignored for the global paper estimate
    pub margin: f64,
    /// Per-pixel paper map (polynomial surface) instead of a single paper colour
    pub flat_field: bool,
    /// Degree of the paper surface
    pub degree: usize,
    /// Largest side of the analysis thumbnail
    pub thumb_side: usize,
    /// Honour embedded ICC profiles (otherwise data are assumed sRGB)
    pub use_icc: bool,
    pub output: OutputOptions,
    /// Write debug images to `<debug>/<image name>/`
    pub debug: Option<PathBuf>,
}

pub struct Report {
    pub width: usize,
    pub height: usize,
    pub paper_lab: Vec3,
    pub paper_l_range: (f64, f64),
    pub icc: bool,
    /// Number of seam edges corrected
    pub seam_edges: usize,
    /// Wall-clock seconds: thumbnail read, paper analysis, full-resolution pass, file finalisation
    pub timings: [f64; 4],
}

/// Process one sheet. `seams`: its layout and the edges it shares with adjacent sheets.
pub fn process_file(input: &Path, output: &Path, opts: &Options, seams: Option<&(Layout, Vec<Edge>)>) -> Result<Report> {
    let src = Input::open(input)?;
    match src.sample_type {
        SampleType::U8 => run::<u8>(&src, output, opts, seams),
        SampleType::U16 => run::<u16>(&src, output, opts, seams),
    }
    .with_context(|| format!("while processing {}", input.display()))
}

type IccTransform<T> = lcms2::Transform<[T; 3], [f32; 3], GlobalContext, DisallowCache>;

/// Linear-light sRGB (D65) profile, the working space after ICC conversion.
fn linear_srgb_profile() -> Result<Profile> {
    let xy = |x, y| CIExyY { x, y, Y: 1.0 };
    let primaries = CIExyYTRIPLE { Red: xy(0.64, 0.33), Green: xy(0.30, 0.60), Blue: xy(0.15, 0.06) };
    let linear = ToneCurve::new(1.0);
    Ok(Profile::new_rgb(&xy(0.3127, 0.3290), &primaries, &[&linear, &linear, &linear])?)
}

fn icc_transform<I: Copy + lcms2::Pod>(icc: &[u8], format: PixelFormat) -> Result<lcms2::Transform<I, [f32; 3], GlobalContext, DisallowCache>> {
    let src = Profile::new_icc(icc).context("unreadable embedded ICC profile")?;
    Ok(lcms2::Transform::new_flags_context(
        GlobalContext::new(),
        &src,
        format,
        &linear_srgb_profile()?,
        PixelFormat::RGB_FLT,
        Intent::RelativeColorimetric,
        Flags::NO_CACHE,
    )?)
}

/// Integer samples -> linear RGB.
enum Decoder<T: Sample> {
    Lut(Vec<f32>),
    Icc(IccTransform<T>),
}

impl<T: Sample> Decoder<T> {
    fn decode_row(&self, row: &[T], bands: usize, out: &mut [[f32; 3]]) {
        match self {
            Decoder::Lut(lut) => {
                for (o, px) in out.iter_mut().zip(row.chunks_exact(bands)) {
                    *o = [lut[px[0].index()], lut[px[1].index()], lut[px[2].index()]];
                }
            }
            Decoder::Icc(t) => {
                let rgb: Vec<[T; 3]> = row.chunks_exact(bands).map(|px| [px[0], px[1], px[2]]).collect();
                t.transform_pixels(&rgb, out);
            }
        }
    }
}

/// Encoded [0, 1] values of a reduced read -> linear RGB.
fn to_linear(px: &[[f32; 3]], icc: Option<&[u8]>) -> Result<Vec<[f32; 3]>> {
    Ok(match icc {
        Some(p) => {
            let t = icc_transform::<[f32; 3]>(p, PixelFormat::RGB_FLT)?;
            let mut out = vec![[0f32; 3]; px.len()];
            t.transform_pixels(px, &mut out);
            out
        }
        None => px.iter().map(|p| p.map(|c| srgb_to_linear(c as f64) as f32)).collect(),
    })
}

/// Largest side of the reduced image used to measure the paper along the seams
const SEAM_IMAGE_SIDE: usize = 1600;

fn run<T: Sample>(src: &Input, output: &Path, opts: &Options, seams: Option<&(Layout, Vec<Edge>)>) -> Result<Report> {
    let t0 = Instant::now();
    let icc = src.icc.as_deref().filter(|_| opts.use_icc);

    // 1. Analysis thumbnail in Lab
    let (tw, th, ts, thumb) = src.read_thumbnail(opts.thumb_side)?;
    let t_thumb = t0.elapsed().as_secs_f64();
    let linear = to_linear(&thumb, icc)?;
    let lab = LabImage {
        width: tw,
        height: th,
        data: linear.iter().map(|p| linear_rgb_to_lab(&p.map(|c| c as f64))).collect(),
    };

    // 2. Paper estimation
    let global = paper_global(&lab, opts.margin);
    let paper = if opts.flat_field {
        paper_surface(&lab, global, opts.degree)
    } else {
        PaperSurface { lab: vec![global; tw * th], role: vec![PixelRole::Other; tw * th] }
    };
    let l_range = paper.lab.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])));
    let target_xyz = srgb_to_xyz_u8(opts.target);
    let surface_gains = GainMap::from_paper_lab(&paper.lab, tw, th, &target_xyz);
    let mut adapter = Adapter::new(surface_gains.clone(), src.width, src.height);

    // 3. Seams with adjacent sheets
    let seams = seams.filter(|(_, e)| !e.is_empty());
    let medium = if seams.is_some() || opts.debug.is_some() {
        let (mw, mh, s, px) = src.read_thumbnail(SEAM_IMAGE_SIDE)?;
        Some((mw, mh, s, to_linear(&px, icc)?))
    } else {
        None
    };
    let mut field = None;
    if let (Some((layout, edges)), Some((mw, mh, s, medium))) = (seams, &medium) {
        let (mw, mh, s) = (*mw, *mh, *s);
        let f = EdgeField::measure(layout, edges, &target_xyz, |px| {
            let (x, y) = (px[0] / s as f64 - 0.5, px[1] / s as f64 - 0.5);
            if !(0.0..=(mw - 1) as f64).contains(&x) || !(0.0..=(mh - 1) as f64).contains(&y) {
                return None;
            }
            let (x0, y0) = (x as usize, y as usize);
            let (x1, y1) = ((x0 + 1).min(mw - 1), (y0 + 1).min(mh - 1));
            let (fx, fy) = (x - x0 as f64, y - y0 as f64);
            let at = |x: usize, y: usize, c: usize| medium[y * mw + x][c] as f64;
            let rgb = [0, 1, 2].map(|c| {
                let top = at(x0, y0, c) + (at(x1, y0, c) - at(x0, y0, c)) * fx;
                let bottom = at(x0, y1, c) + (at(x1, y1, c) - at(x0, y1, c)) * fx;
                top + (bottom - top) * fy
            });
            Some(adapter.adapt_point(px, rgb))
        });
        adapter.scale_gains(src.width, src.height, |px| f.log_gain(px).map(f64::exp));
        field = Some(f);
    }
    if let (Some(dir), Some((mw, mh, s, medium))) = (&opts.debug, &medium) {
        let name = output.file_stem().unwrap_or_default();
        debug::write(
            &dir.join(name),
            &debug::DebugInput {
                full_size: (src.width, src.height),
                global_paper: global,
                thumb: debug::Reduced { width: tw, height: th, factor: ts, pixels: &linear },
                surface: &paper,
                surface_gains: &surface_gains,
                seams: field.as_ref(),
                medium: debug::Reduced { width: *mw, height: *mh, factor: *s, pixels: medium },
                adapter: &adapter,
            },
        )?;
    }
    let seam_edges = field.as_ref().map_or(0, EdgeField::edges);
    let t_analysis = t0.elapsed().as_secs_f64();

    // 4. Streaming full-resolution pass
    let decoder: Decoder<T> = match icc {
        Some(p) => Decoder::Icc(icc_transform::<[T; 3]>(p, T::ICC_FORMAT)?),
        None => Decoder::Lut(decode_table(T::MAX)),
    };
    let encoder = SrgbEncoder::new();
    let (w, bands) = (src.width, src.bands);
    let row_len = w * bands;
    let strip_rows = ((64 << 20) / (row_len * std::mem::size_of::<T>())).clamp(1, src.height);
    let mut buf = vec![T::default(); strip_rows * row_len];
    let mut dst = Output::create::<T>(src, output, &opts.output)?;

    for y0 in (0..src.height).step_by(strip_rows) {
        let rows = strip_rows.min(src.height - y0);
        let strip = &mut buf[..rows * row_len];
        src.read_rows(y0, rows, strip)?;
        strip.par_chunks_mut(row_len).enumerate().for_each_init(
            || vec![[0f32; 3]; w],
            |lin, (dy, row)| {
                decoder.decode_row(row, bands, lin);
                adapter.apply_row(lin, &adapter.row_gains(y0 + dy));
                for (px, l) in row.chunks_exact_mut(bands).zip(lin.iter()) {
                    for c in 0..3 {
                        px[c] = T::from_unit(encoder.encode(l[c]));
                    }
                }
            },
        );
        dst.write_rows(y0, rows, w, strip)?;
    }
    let t_pass = t0.elapsed().as_secs_f64();
    dst.finish()?;
    let t_end = t0.elapsed().as_secs_f64();

    Ok(Report {
        width: src.width,
        height: src.height,
        paper_lab: global,
        paper_l_range: l_range,
        icc: icc.is_some(),
        seam_edges,
        timings: [t_thumb, t_analysis - t_thumb, t_pass - t_analysis, t_end - t_pass],
    })
}
