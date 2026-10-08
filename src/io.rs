//! Raster I/O through GDAL: pixel-interleaved strip reads/writes, downsampled
//! thumbnail reads, ICC profile extraction and georeferencing copy.

use std::ffi::{c_int, c_void};
use std::path::{Path, PathBuf};

use anyhow::{bail, ensure, Context, Result};
use base64::Engine;
use gdal::raster::{GdalDataType, GdalType, RasterCreationOptions, ResampleAlg};
use gdal::{Dataset, DriverManager, Metadata};
use gdal_sys::{CPLErr, GDALRWFlag};

/// Integer sample types supported (8 and 16 bits).
pub trait Sample: Copy + Default + Send + Sync + GdalType + lcms2::Pod {
    const MAX: u32;
    const ICC_FORMAT: lcms2::PixelFormat;
    fn index(self) -> usize;
    fn from_unit(v: f32) -> Self;
}

impl Sample for u8 {
    const MAX: u32 = 255;
    const ICC_FORMAT: lcms2::PixelFormat = lcms2::PixelFormat::RGB_8;
    #[inline(always)]
    fn index(self) -> usize {
        self as usize
    }
    #[inline(always)]
    fn from_unit(v: f32) -> Self {
        (v * 255.0 + 0.5) as u8
    }
}

impl Sample for u16 {
    const MAX: u32 = 65535;
    const ICC_FORMAT: lcms2::PixelFormat = lcms2::PixelFormat::RGB_16;
    #[inline(always)]
    fn index(self) -> usize {
        self as usize
    }
    #[inline(always)]
    fn from_unit(v: f32) -> Self {
        (v * 65535.0 + 0.5) as u16
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleType {
    U8,
    U16,
}

pub struct Input {
    pub ds: Dataset,
    pub width: usize,
    pub height: usize,
    /// 3 (RGB) or 4 (RGB + alpha)
    pub bands: usize,
    pub sample_type: SampleType,
    pub icc: Option<Vec<u8>>,
}

impl Input {
    pub fn open(path: &Path) -> Result<Self> {
        let ds = Dataset::open(path).with_context(|| format!("cannot open {}", path.display()))?;
        let bands = ds.raster_count();
        ensure!(
            bands == 3 || bands == 4,
            "{}: {bands} band(s), only RGB or RGBA images are supported",
            path.display()
        );
        let types: Vec<GdalDataType> =
            (1..=bands).map(|b| ds.rasterband(b).map(|rb| rb.band_type())).collect::<Result<_, _>>()?;
        let sample_type = match types[0] {
            GdalDataType::UInt8 => SampleType::U8,
            GdalDataType::UInt16 => SampleType::U16,
            t => bail!("{}: unsupported sample type {t:?} (8 or 16-bit integers only)", path.display()),
        };
        ensure!(types.iter().all(|t| *t == types[0]), "{}: bands have different types", path.display());
        let icc = ds
            .metadata_item("SOURCE_ICC_PROFILE", "COLOR_PROFILE")
            .map(|b64| base64::engine::general_purpose::STANDARD.decode(b64.trim()))
            .transpose()
            .context("invalid embedded ICC profile")?;
        let (width, height) = ds.raster_size();
        Ok(Self { ds, width, height, bands, sample_type, icc })
    }

    /// Read rows `y0..y0+rows` of all bands into a pixel-interleaved buffer.
    pub fn read_rows<T: Sample>(&self, y0: usize, rows: usize, buf: &mut [T]) -> Result<()> {
        raster_io(&self.ds, GDALRWFlag::GF_Read, (0, y0), (self.width, rows), buf, (self.width, rows), self.bands, None)
    }

    /// Downsampled (box average) read of the RGB bands, with values in [0, 1].
    /// Returns (width, height, reduction factor, pixels).
    pub fn read_thumbnail(&self, max_side: usize) -> Result<(usize, usize, usize, Vec<[f32; 3]>)> {
        let s = self.width.max(self.height).div_ceil(max_side).max(1);
        let (tw, th) = (self.width / s, self.height / s);
        let mut buf = vec![0f32; tw * th * 3];
        raster_io(
            &self.ds,
            GDALRWFlag::GF_Read,
            (0, 0),
            (tw * s, th * s),
            &mut buf,
            (tw, th),
            3,
            Some(ResampleAlg::Average),
        )?;
        let k = 1.0 / match self.sample_type {
            SampleType::U8 => u8::MAX as f32,
            SampleType::U16 => u16::MAX as f32,
        };
        Ok((tw, th, s, buf.chunks_exact(3).map(|p| [p[0] * k, p[1] * k, p[2] * k]).collect()))
    }
}

/// Pixel-interleaved dataset RasterIO over the first `bands` bands.
#[allow(clippy::too_many_arguments)]
fn raster_io<T: GdalType>(
    ds: &Dataset,
    flag: GDALRWFlag::Type,
    (x0, y0): (usize, usize),
    (w, h): (usize, usize),
    buf: &mut [T],
    (bw, bh): (usize, usize),
    bands: usize,
    resample: Option<ResampleAlg>,
) -> Result<()> {
    ensure!(buf.len() == bw * bh * bands, "buffer size mismatch");
    let size = std::mem::size_of::<T>() as i64;
    let mut band_map: Vec<c_int> = (1..=bands as c_int).collect();
    let mut extra = gdal_sys::GDALRasterIOExtraArg {
        nVersion: 1,
        eResampleAlg: resample.unwrap_or(ResampleAlg::NearestNeighbour).to_gdal(),
        pfnProgress: None,
        pProgressData: std::ptr::null_mut(),
        bFloatingPointWindowValidity: 0,
        dfXOff: 0.0,
        dfYOff: 0.0,
        dfXSize: 0.0,
        dfYSize: 0.0,
        bUseOnlyThisScale: 0,
    };
    let rv = unsafe {
        gdal_sys::GDALDatasetRasterIOEx(
            ds.c_dataset(),
            flag,
            x0 as c_int,
            y0 as c_int,
            w as c_int,
            h as c_int,
            buf.as_mut_ptr() as *mut c_void,
            bw as c_int,
            bh as c_int,
            T::gdal_ordinal(),
            bands as c_int,
            band_map.as_mut_ptr(),
            size * bands as i64,
            size * bands as i64 * bw as i64,
            size,
            &mut extra,
        )
    };
    ensure!(rv == CPLErr::CE_None, "GDAL RasterIO failed: {}", last_gdal_error());
    Ok(())
}

fn last_gdal_error() -> String {
    unsafe {
        let msg = gdal_sys::CPLGetLastErrorMsg();
        if msg.is_null() { String::new() } else { std::ffi::CStr::from_ptr(msg).to_string_lossy().into_owned() }
    }
}

/// Output raster format.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    /// Same family as the input (TIFF -> TIFF, JPEG -> JPEG, PNG -> PNG, anything else -> TIFF)
    Auto,
    Tif,
    Png,
    Jpg,
}

impl Format {
    /// Format of a file, from its extension (anything else than PNG or JPEG is TIFF).
    pub fn of_path(path: &Path) -> Format {
        match path.extension().and_then(|e| e.to_str()).map(str::to_lowercase).as_deref() {
            Some("jpg" | "jpeg") => Format::Jpg,
            Some("png") => Format::Png,
            _ => Format::Tif,
        }
    }

    /// Format to write for `input` (Auto: same family as the input).
    pub fn for_input(self, input: &Path) -> Format {
        if self == Format::Auto { Format::of_path(input) } else { self }
    }

    pub fn extension(self) -> &'static str {
        match self {
            Format::Png => "png",
            Format::Jpg => "jpg",
            _ => "tif",
        }
    }

    fn driver(self) -> &'static str {
        match self {
            Format::Png => "PNG",
            Format::Jpg => "JPEG",
            _ => "GTiff",
        }
    }
}

pub struct OutputOptions {
    pub jpeg_quality: u8,
    pub tiff_compress: String,
}

/// Destination raster, its format given by the extension of its path. GTiff is written strip
/// by strip; JPEG/PNG only support CreateCopy, so they go through an in-memory dataset first.
pub struct Output {
    ds: Dataset,
    path: PathBuf,
    bands: usize,
    copy_to: Option<(&'static str, Vec<String>)>,
}

impl Output {
    pub fn create<T: Sample>(input: &Input, path: &Path, opts: &OutputOptions) -> Result<Self> {
        let format = Format::of_path(path);
        let (w, h, bands) = (input.width, input.height, input.bands);
        ensure!(
            !(format == Format::Jpg && T::MAX > 255),
            "16-bit images cannot be written as JPEG, use --format tif or png"
        );
        let mut options = Vec::new();
        let (ds, copy_to) = match format {
            Format::Tif => {
                options.extend([
                    format!("COMPRESS={}", opts.tiff_compress),
                    "BIGTIFF=IF_SAFER".into(),
                    "PHOTOMETRIC=RGB".into(),
                    "INTERLEAVE=PIXEL".into(),
                    "TILED=YES".into(),
                    "BLOCKXSIZE=512".into(),
                    "BLOCKYSIZE=512".into(),
                    "NUM_THREADS=ALL_CPUS".into(),
                ]);
                if ["DEFLATE", "LZW", "ZSTD"].contains(&opts.tiff_compress.to_uppercase().as_str()) {
                    options.push("PREDICTOR=2".into());
                }
                if bands == 4 {
                    options.push("ALPHA=YES".into());
                }
                let co = RasterCreationOptions::from_iter(options.iter().map(String::as_str));
                let ds = DriverManager::get_driver_by_name("GTiff")?
                    .create_with_band_type_with_options::<T, _>(path, w, h, bands, &co)
                    .with_context(|| format!("cannot create {}", path.display()))?;
                (ds, None)
            }
            _ => {
                if format == Format::Jpg {
                    options.push(format!("QUALITY={}", opts.jpeg_quality));
                }
                let ds = DriverManager::get_driver_by_name("MEM")?.create_with_band_type::<T, _>("", w, h, bands)?;
                (ds, Some((format.driver(), options)))
            }
        };
        let mut out = Self { ds, path: path.to_owned(), bands, copy_to };
        out.copy_georeferencing(&input.ds)?;
        Ok(out)
    }

    fn copy_georeferencing(&mut self, src: &Dataset) -> Result<()> {
        if let Ok(gt) = src.geo_transform() {
            self.ds.set_geo_transform(&gt)?;
        }
        let proj = src.projection();
        if !proj.is_empty() {
            self.ds.set_projection(&proj)?;
        }
        unsafe {
            let (s, d) = (src.c_dataset(), self.ds.c_dataset());
            let n = gdal_sys::GDALGetGCPCount(s);
            if n > 0 {
                let rv = gdal_sys::GDALSetGCPs(d, n, gdal_sys::GDALGetGCPs(s), gdal_sys::GDALGetGCPProjection(s));
                ensure!(rv == CPLErr::CE_None, "cannot copy GCPs: {}", last_gdal_error());
            }
        }
        // EXIF tags describe the original file (dimensions, software...): not copied
        for entry in src.metadata() {
            if entry.domain.is_empty() && !entry.key.starts_with("EXIF_") {
                self.ds.set_metadata_item(&entry.key, &entry.value, "")?;
            }
        }
        for b in 1..=self.bands {
            let sb = src.rasterband(b)?;
            let mut db = self.ds.rasterband(b)?;
            db.set_color_interpretation(sb.color_interpretation())?;
            if let Some(nd) = sb.no_data_value() {
                db.set_no_data_value(Some(nd))?;
            }
        }
        Ok(())
    }

    pub fn write_rows<T: Sample>(&mut self, y0: usize, rows: usize, width: usize, buf: &mut [T]) -> Result<()> {
        raster_io(&self.ds, GDALRWFlag::GF_Write, (0, y0), (width, rows), buf, (width, rows), self.bands, None)
    }

    pub fn finish(self) -> Result<()> {
        if let Some((driver, options)) = &self.copy_to {
            let co = RasterCreationOptions::from_iter(options.iter().map(String::as_str));
            let driver = DriverManager::get_driver_by_name(driver)?;
            self.ds
                .create_copy(&driver, &self.path, &co)
                .with_context(|| format!("cannot write {}", self.path.display()))?
                .close()?;
        }
        self.ds.close()?;
        Ok(())
    }
}

/// Write an 8-bit RGB image as PNG (debug output).
pub fn write_png(path: &Path, width: usize, height: usize, pixels: &[[u8; 3]]) -> Result<()> {
    let mem = DriverManager::get_driver_by_name("MEM")?.create_with_band_type::<u8, _>("", width, height, 3)?;
    let mut buf: Vec<u8> = pixels.iter().flatten().copied().collect();
    raster_io(&mem, GDALRWFlag::GF_Write, (0, 0), (width, height), &mut buf, (width, height), 3, None)?;
    for b in 1..=3 {
        mem.rasterband(b)?.set_color_interpretation(match b {
            1 => gdal::raster::ColorInterpretation::RedBand,
            2 => gdal::raster::ColorInterpretation::GreenBand,
            _ => gdal::raster::ColorInterpretation::BlueBand,
        })?;
    }
    let png = DriverManager::get_driver_by_name("PNG")?;
    mem.create_copy(&png, path, &RasterCreationOptions::new())
        .with_context(|| format!("cannot write {}", path.display()))?
        .close()?;
    Ok(())
}
