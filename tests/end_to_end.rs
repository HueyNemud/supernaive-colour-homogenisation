//! End-to-end check on a synthetic 16-bit GeoTIFF sheet: yellowish paper with a
//! strong illumination gradient, neutral ink lines and a blue wash.

use std::path::Path;

use gdal::raster::RasterCreationOptions;
use gdal::spatial_ref::SpatialRef;
use gdal::{Dataset, DriverManager};

use homog::color::{linear_rgb_to_lab, linear_to_srgb, srgb_to_linear, Vec3};
use homog::io::{Format, OutputOptions};
use homog::{process_file, Options};

const W: usize = 600;
const H: usize = 400;
const GT: [f64; 6] = [650000.0, 0.5, 0.0, 6865000.0, 0.0, -0.5];

/// Linear RGB of the synthetic sheet at (x, y), and whether it is bare paper.
fn sheet(x: usize, y: usize) -> ([f64; 3], bool) {
    let paper = [0.62, 0.52, 0.33]; // yellowish, L* ~ 78
    let light = 1.0 - 0.35 * (x as f64 / W as f64); // darker towards the right edge
    let (reflectance, is_paper) = if x % 50 < 3 || y % 50 < 3 {
        ([0.05, 0.05, 0.05], false) // neutral ink (relative to white)
    } else if (200..260).contains(&x) && (100..300).contains(&y) {
        ([0.55, 0.75, 0.95], false) // blue wash
    } else {
        ([1.0, 1.0, 1.0], true)
    };
    ([0, 1, 2].map(|c| paper[c] * light * reflectance[c]), is_paper)
}

fn write_input(path: &Path) {
    let mut buf = vec![0u16; W * H * 3];
    for y in 0..H {
        for x in 0..W {
            let (lin, _) = sheet(x, y);
            for c in 0..3 {
                buf[(y * W + x) * 3 + c] = (linear_to_srgb(lin[c]) * 65535.0).round() as u16;
            }
        }
    }
    let co = RasterCreationOptions::from_iter(["INTERLEAVE=PIXEL"]);
    let mut ds = DriverManager::get_driver_by_name("GTiff")
        .unwrap()
        .create_with_band_type_with_options::<u16, _>(path, W, H, 3, &co)
        .unwrap();
    ds.set_geo_transform(&GT).unwrap();
    ds.set_spatial_ref(&SpatialRef::from_epsg(2154).unwrap()).unwrap();
    for c in 0..3 {
        let mut band_buf = gdal::raster::Buffer::new((W, H), buf.iter().skip(c).step_by(3).copied().collect());
        ds.rasterband(c + 1).unwrap().write((0, 0), (W, H), &mut band_buf).unwrap();
    }
}

fn read_lab(path: &Path) -> Vec<Vec3> {
    let ds = Dataset::open(path).unwrap();
    assert_eq!(ds.rasterband(1).unwrap().band_type(), gdal::raster::GdalDataType::UInt16);
    let bands: Vec<Vec<u16>> = (1..=3)
        .map(|b| ds.rasterband(b).unwrap().read_band_as::<u16>().unwrap().into_shape_and_vec().1)
        .collect();
    (0..W * H)
        .map(|i| linear_rgb_to_lab(&[0, 1, 2].map(|c| srgb_to_linear(bands[c][i] as f64 / 65535.0))))
        .collect()
}

fn run(flat_field: bool) -> (Vec<Vec3>, Dataset) {
    let dir = std::env::temp_dir().join(format!("homog-e2e-{}-{flat_field}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (input, output) = (dir.join("in.tif"), dir.join("out.tif"));
    write_input(&input);
    let opts = Options {
        target: [255, 255, 255],
        margin: 0.05,
        flat_field,
        degree: 3,
        thumb_side: 400,
        use_icc: true,
        output: OutputOptions { format: Format::Tif, jpeg_quality: 95, tiff_compress: "DEFLATE".into() },
        debug: None,
    };
    process_file(&input, &output, &opts, None).unwrap();
    (read_lab(&output), Dataset::open(&output).unwrap())
}

fn paper_l_range(lab: &[Vec3]) -> (f64, f64) {
    // Bare paper away from the sheet border, where the estimate is well supported
    let mut lo: f64 = 100.0;
    let mut hi: f64 = 0.0;
    for y in 20..H - 20 {
        for x in 20..W - 20 {
            if sheet(x, y).1 && x % 50 > 8 && y % 50 > 8 {
                lo = lo.min(lab[y * W + x][0]);
                hi = hi.max(lab[y * W + x][0]);
            }
        }
    }
    (lo, hi)
}

#[test]
fn flat_field_whitens_the_whole_sheet() {
    let (lab, ds) = run(true);
    let (lo, hi) = paper_l_range(&lab);
    assert!(lo > 97.0 && hi <= 100.0 + 1e-6, "paper L* in [{lo:.2}, {hi:.2}]");

    // Paper chroma removed, ink stays neutral, blue wash stays blue
    let at = |x: usize, y: usize| lab[y * W + x];
    let p = at(120, 120);
    assert!(p[1].abs() < 1.0 && p[2].abs() < 1.0, "paper {p:?}");
    let ink = at(301, 200);
    assert!(ink[1].abs() < 1.5 && ink[2].abs() < 1.5, "ink {ink:?}");
    // The wash should look as it would on white paper (its reflectance), within ~2 ΔE
    let wash = at(230, 220);
    let expected = linear_rgb_to_lab(&[0.55, 0.75, 0.95]);
    let de = (0..3).map(|c| (wash[c] - expected[c]).powi(2)).sum::<f64>().sqrt();
    assert!(de < 2.0, "wash {wash:?}, expected {expected:?}, ΔE {de:.2}");

    // Georeferencing preserved
    assert_eq!(ds.geo_transform().unwrap(), GT);
    assert_eq!(ds.spatial_ref().unwrap().auth_code().unwrap(), 2154);
}

#[test]
fn global_mode_leaves_illumination_gradient() {
    let (lab, _) = run(false);
    let (lo, hi) = paper_l_range(&lab);
    // The 35 % light fall-off is not corrected by a single paper colour
    assert!(hi - lo > 10.0, "paper L* in [{lo:.2}, {hi:.2}]");
}

/// Two sheets side by side (map x 0..600 and 600..1200, y 0..400, 1 px per unit). The west
/// sheet's paper darkens sharply over its last 20 px before the seam, which the paper surface
/// cannot follow; the seam correction must bring it back to white.
mod seams {
    use super::*;
    use homog::layout::Atlas;

    const SW: usize = 600;
    const SH: usize = 400;

    fn west_paper(x: usize) -> f64 {
        let dark = ((x as f64 - 580.0) / 20.0).clamp(0.0, 1.0) * 0.25;
        0.55 * (1.0 - dark)
    }

    fn write_sheet(path: &Path, paper: impl Fn(usize) -> f64) {
        let ds = DriverManager::get_driver_by_name("GTiff").unwrap().create_with_band_type::<u8, _>(path, SW, SH, 3).unwrap();
        for (b, tint) in [1.0, 0.92, 0.7].iter().enumerate() {
            // Ink lines every 50 px so that the sheet is not blank
            let data: Vec<u8> = (0..SW * SH)
                .map(|i| {
                    let (x, y) = (i % SW, i / SW);
                    let v = if x % 50 == 25 || y % 50 == 25 { 0.03 } else { paper(x) * tint };
                    (linear_to_srgb(v) * 255.0).round() as u8
                })
                .collect();
            let mut buf = gdal::raster::Buffer::new((SW, SH), data);
            ds.rasterband(b + 1).unwrap().write((0, 0), (SW, SH), &mut buf).unwrap();
        }
    }

    fn write_points(path: &Path, x0: f64) {
        let mut s = String::from("#CRS: LOCAL\nmapX,mapY,sourceX,sourceY,enable,dX,dY,residual\n");
        for (px, line) in [(0.0, 0.0), (600.0, 0.0), (0.0, 400.0), (600.0, 400.0)] {
            s += &format!("{},{},{px},{},1,0,0,0\n", x0 + px, 400.0 - line, -line);
        }
        std::fs::write(path, s).unwrap();
    }

    /// Mean L* of the west sheet's output paper in a column band [x0, x1).
    fn paper_l(path: &Path, x0: usize, x1: usize) -> f64 {
        let ds = Dataset::open(path).unwrap();
        let bands: Vec<Vec<u8>> = (1..=3).map(|b| ds.rasterband(b).unwrap().read_band_as::<u8>().unwrap().into_shape_and_vec().1).collect();
        let mut sum = 0.0;
        let mut n = 0.0;
        for y in (0..SH).filter(|y| y % 50 != 25 && y % 50 != 24 && y % 50 != 26) {
            for x in x0..x1 {
                let i = y * SW + x;
                sum += linear_rgb_to_lab(&[0, 1, 2].map(|c| srgb_to_linear(bands[c][i] as f64 / 255.0)))[0];
                n += 1.0;
            }
        }
        sum / n
    }

    #[test]
    fn seam_correction_whitens_the_paper_along_the_shared_edge() {
        let dir = std::env::temp_dir().join(format!("homog-seams-{}", std::process::id()));
        let gcp = dir.join("gcp");
        std::fs::create_dir_all(&gcp).unwrap();
        let (west, east) = (dir.join("west.tif"), dir.join("east.tif"));
        write_sheet(&west, west_paper);
        write_sheet(&east, |_| 0.55);
        write_points(&gcp.join("west.tif.points"), 0.0);
        write_points(&gcp.join("east.tif.points"), 600.0);

        let atlas = Atlas::load(&gcp).unwrap();
        let seams = atlas.seams("west.tif").unwrap();
        assert_eq!(seams.1.len(), 1, "one shared edge");

        let opts = |out: &str| (dir.join(out), Options {
            target: [255, 255, 255],
            margin: 0.05,
            flat_field: true,
            degree: 3,
            thumb_side: 400,
            use_icc: true,
            output: OutputOptions { format: Format::Tif, jpeg_quality: 95, tiff_compress: "NONE".into() },
            debug: None,
        });
        let (plain, o1) = opts("plain.tif");
        let (seamed, o2) = opts("seamed.tif");
        process_file(&west, &plain, &o1, None).unwrap();
        let report = process_file(&west, &seamed, &o2, Some(&seams)).unwrap();
        assert_eq!(report.seam_edges, 1);

        // Measurement band: 2 to 8 px from the seam
        let (before, after) = (paper_l(&plain, 592, 598), paper_l(&seamed, 592, 598));
        assert!(after > 97.5, "paper L* along the seam: {before:.1} without, {after:.1} with seam correction");
        assert!(after - before > 3.0, "paper L* along the seam: {before:.1} without, {after:.1} with seam correction");
        // Far from the seam, nothing changes
        let (mid_before, mid_after) = (paper_l(&plain, 200, 300), paper_l(&seamed, 200, 300));
        assert!((mid_after - mid_before).abs() < 0.05, "centre: {mid_before:.2} vs {mid_after:.2}");
    }
}
