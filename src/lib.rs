//! Not-So-Naive Map Colour Homogenisation: white, seamless paper for scanned map sheets.
//! See `docs/technical-report.md` for the model.
//!
//! Pipeline per image:
//! 1. downsampled read -> Lab thumbnail;
//! 2. paper estimation: content mask, paper colour model (automatic or from `paper`
//!    imagettes, `keep` imagettes excluded), then a smooth paper surface (flat-field);
//! 3. seams (optional, layout from GCPs): along the edges shared with adjacent sheets, the
//!    remaining paper colour is measured on a reduced image and corrected, fading inwards;
//! 4. one streaming pass over full-resolution strips: decode (sRGB LUT or ICC),
//!    Bradford adaptation of the local paper colour to the target white, sRGB encode.

pub mod background;
pub mod calibrate;
pub mod color;
pub mod config;
pub mod debug;
pub mod hints;
pub mod io;
pub mod layout;
pub mod seams;
pub mod transform;

use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use lcms2::{CIExyY, CIExyYTRIPLE, DisallowCache, Flags, GlobalContext, Intent, PixelFormat, Profile, ToneCurve};
use rayon::prelude::*;

use background::{paper_surface, LabImage, Lighting, PaperModel, PaperSurface, PixelRole};
use color::{decode_table, linear_rgb_to_lab, srgb_to_linear, srgb_to_xyz_u8, srgb_u8_to_linear, SrgbEncoder, Vec3};
use hints::Hints;
use io::{Input, Output, OutputOptions, Sample, SampleType};
use layout::Seams;
use seams::EdgeField;
use transform::{Adapter, GainMap};

pub struct Options {
    /// Target background colour (sRGB, 8-bit)
    pub target: [u8; 3],
    pub lighting: Lighting,
    /// Honour embedded ICC profiles (otherwise data are assumed sRGB)
    pub use_icc: bool,
    pub hints: Hints,
    pub output: OutputOptions,
    /// Write debug images to `<debug>/<image name>/`
    pub debug: Option<PathBuf>,
}

/// Largest side of the analysis thumbnail
const THUMB_SIDE: usize = 400;
/// Largest side of the reduced image used to measure the paper along the seams
const SEAM_IMAGE_SIDE: usize = 1600;

/// A reduced image in linear RGB and its reduction factor w.r.t. the full resolution.
pub struct Reduced {
    pub width: usize,
    pub height: usize,
    pub factor: usize,
    pub pixels: Vec<[f32; 3]>,
}

impl Reduced {
    /// Bilinear sample at a full-resolution pixel position, None outside.
    pub fn sample(&self, px: [f64; 2]) -> Option<[f64; 3]> {
        let (w, h, s) = (self.width, self.height, self.factor as f64);
        let (x, y) = (px[0] / s - 0.5, px[1] / s - 0.5);
        if !(0.0..=(w - 1) as f64).contains(&x) || !(0.0..=(h - 1) as f64).contains(&y) {
            return None;
        }
        let (x0, y0) = (x as usize, y as usize);
        let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
        let (fx, fy) = (x - x0 as f64, y - y0 as f64);
        let at = |x: usize, y: usize, c: usize| self.pixels[y * w + x][c] as f64;
        Some([0, 1, 2].map(|c| {
            let top = at(x0, y0, c) + (at(x1, y0, c) - at(x0, y0, c)) * fx;
            let bottom = at(x0, y1, c) + (at(x1, y1, c) - at(x0, y1, c)) * fx;
            top + (bottom - top) * fy
        }))
    }
}

/// Everything estimated on a sheet before the full-resolution pass.
pub struct Analysis {
    pub width: usize,
    pub height: usize,
    pub icc: bool,
    pub thumb: Reduced,
    pub paper: PaperModel,
    pub surface: PaperSurface,
    /// Gains of the surface correction alone
    pub surface_gains: GainMap,
    /// Final adaptation (surface and seams)
    pub adapter: Adapter,
    pub seams: Option<EdgeField>,
    /// Reduced image, read when seams or debug images are needed
    pub medium: Option<Reduced>,
    /// Seconds: thumbnail read, analysis
    pub timings: [f64; 2],
}

impl Analysis {
    /// Share of the content pixels taken as bare paper.
    pub fn paper_fraction(&self) -> f64 {
        let content = self.paper.content.iter().filter(|&&c| c).count().max(1);
        self.surface.role.iter().filter(|&&r| r == PixelRole::Chroma).count() as f64 / content as f64
    }

    pub fn paper_l_range(&self) -> (f64, f64) {
        self.surface.lab.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])))
    }
}

pub struct Report {
    pub width: usize,
    pub height: usize,
    pub paper_lab: Vec3,
    pub paper_l_range: (f64, f64),
    pub degree: Option<usize>,
    pub icc: bool,
    /// Number of seam edges corrected
    pub seam_edges: usize,
    /// Wall-clock seconds: thumbnail read, paper analysis, full-resolution pass, file finalisation
    pub timings: [f64; 4],
}

/// Process one sheet. `seams`: its layout and the edges it shares with adjacent sheets.
pub fn process_file(input: &Path, output: &Path, opts: &Options, seams: Option<&Seams>) -> Result<Report> {
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

/// Analyse a sheet: thumbnail, paper model and surface, seams with adjacent sheets.
pub fn analyse(src: &Input, opts: &Options, seams: Option<&Seams>) -> Result<Analysis> {
    let t0 = Instant::now();
    let icc = src.icc.as_deref().filter(|_| opts.use_icc);

    // 1. Analysis thumbnail in Lab
    let (tw, th, ts, thumb) = src.read_thumbnail(THUMB_SIDE)?;
    let t_thumb = t0.elapsed().as_secs_f64();
    let thumb = Reduced { width: tw, height: th, factor: ts, pixels: to_linear(&thumb, icc)? };
    let lab = LabImage {
        width: tw,
        height: th,
        data: thumb.pixels.iter().map(|p| linear_rgb_to_lab(&p.map(|c| c as f64))).collect(),
    };

    // 2. Paper model and surface
    let paper = PaperModel::estimate(&lab, &opts.hints)?;
    let surface = paper_surface(&lab, &paper, &opts.hints, opts.lighting.degree());
    let target_xyz = srgb_to_xyz_u8(opts.target);
    let surface_gains = GainMap::from_paper_lab(&surface.lab, tw, th, &target_xyz);
    let mut adapter = Adapter::new(surface_gains.clone(), src.width, src.height);

    // 3. Seams with adjacent sheets
    let seams = seams.filter(|(_, e)| !e.is_empty());
    let medium = if seams.is_some() || opts.debug.is_some() {
        let (mw, mh, s, px) = src.read_thumbnail(SEAM_IMAGE_SIDE)?;
        Some(Reduced { width: mw, height: mh, factor: s, pixels: to_linear(&px, icc)? })
    } else {
        None
    };
    let mut field = None;
    if let (Some((layout, edges)), Some(medium)) = (seams, &medium) {
        let target_lab = linear_rgb_to_lab(&srgb_u8_to_linear(opts.target));
        let tolerance = (2.0 * paper.sd_ab()).max(4.0);
        let is_paper = |raw: &[f64; 3], corrected: &[f64; 3]| {
            let lab = linear_rgb_to_lab(corrected);
            let dab = (lab[1] - target_lab[1]).hypot(lab[2] - target_lab[2]);
            dab < tolerance && lab[0] > target_lab[0] - 20.0 && !opts.hints.keeps(&linear_rgb_to_lab(raw))
        };
        let f = EdgeField::measure(layout, edges, &target_xyz, is_paper, |px: [f64; 2]| {
            medium.sample(px).map(|raw| (raw, adapter.adapt_point(px, raw)))
        });
        adapter.scale_gains(src.width, src.height, |px| f.log_gain(px).map(f64::exp));
        field = Some(f);
    }

    Ok(Analysis {
        width: src.width,
        height: src.height,
        icc: icc.is_some(),
        thumb,
        paper,
        surface,
        surface_gains,
        adapter,
        seams: field,
        medium,
        timings: [t_thumb, t0.elapsed().as_secs_f64() - t_thumb],
    })
}

fn run<T: Sample>(src: &Input, output: &Path, opts: &Options, seams: Option<&Seams>) -> Result<Report> {
    let t0 = Instant::now();
    let icc = src.icc.as_deref().filter(|_| opts.use_icc);
    let a = analyse(src, opts, seams)?;
    if let Some(dir) = &opts.debug {
        debug::write(&dir.join(output.file_stem().unwrap_or_default()), &a)?;
    }
    let adapter = &a.adapter;
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
        paper_lab: a.paper.global,
        paper_l_range: a.paper_l_range(),
        degree: a.surface.degree,
        icc: a.icc,
        seam_edges: a.seams.as_ref().map_or(0, EdgeField::edges),
        timings: [a.timings[0], t_analysis - a.timings[0], t_pass - t_analysis, t_end - t_pass],
    })
}
