# Not-So-Naive Map Colour Homogenisation

[![Licence: AGPL v3](https://img.shields.io/badge/licence-AGPL%20v3-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-%E2%89%A5%201.85-orange.svg)](https://www.rust-lang.org)
[![GDAL](https://img.shields.io/badge/GDAL-3.x-green.svg)](https://gdal.org)

`homog` makes the paper of scanned map sheets white, evenly across each sheet and across sheets,
while keeping the inks and watercolour washes as they were drawn. It is meant as a pre-processing
step before georeferencing and mosaicking maps made of several sheets, so that the assembled map
does not look like a patchwork.

![Three adjacent sheets, as scanned (top) and corrected (bottom)](docs/img/before-after.jpg)

## Features

- **White paper, kept colours.** The paper is white-balanced like the light of a photograph
  (Bradford chromatic adaptation): black stays black, pale washes are not bleached.
- **Uneven lighting corrected.** The paper colour is estimated as a smooth surface over the sheet,
  which follows vignetting and yellowed edges but not the washes.
- **Seamless mosaics.** Given the ground control points of the sheets (QGIS georeferencer files),
  the paper is matched along the edges shared by adjacent sheets. Images are not resampled.
- **Automatic, with an escape hatch.** Everything that can be measured is; when the tool is fooled,
  show it a small crop of paper, or of a colour that is not paper.
- **Explains itself.** Debug images per sheet, and a calibration report that flags suspicious sheets
  of a batch.
- **Fast.** About half a second for a 5800 × 4000 sheet; 8 or 16 bits, JPEG / TIFF / GeoTIFF / PNG;
  ICC profiles honoured; georeferencing kept.

## Install

Requirements: Rust ≥ 1.85, GDAL 3.x with its development files, and `libclang`. On Debian/Ubuntu:

```bash
sudo apt install libgdal-dev libclang-dev
cargo build --release          # binary: target/release/homog
```

## Quick start

The repository comes with three adjacent sheets of an atlas of Paris (Service historique de la
Défense, GR 6 M J10 C 1188) and the ground control points of the whole atlas:

```text
example/
├── SHDGR__GR_6_M_J10_C_1188_001__0017.jpg    three adjacent sheets
├── SHDGR__GR_6_M_J10_C_1188_001__0019.jpg
├── SHDGR__GR_6_M_J10_C_1188_001__0021.jpg
├── GCP/                                       QGIS .points files of the 59 sheets of the atlas
└── homog.toml                                 settings of the example
```

```bash
target/release/homog -v -c example/homog.toml example/*.jpg
```

The corrected sheets are written to `output/`. Add `--debug` to see what the tool did (see
[Debug images](#debug-images)).

## Use it on your own scans

**A single sheet, or sheets that are not assembled:**

```bash
homog scans/*.tif                    # corrected images in output/, same format as the input
homog -o corrected -f png scans/*.jpg
```

**Several sheets of one map:** georeference them in QGIS first (or at least place their ground
control points on the neatline), and point `homog` at the directory of `.points` files. The sheets
that touch are found from the GCPs, and the paper is matched along their shared edges:

```bash
homog -g gcp/ scans/*.jpg            # gcp/<image file name>.points
```

The GCPs only tell where the edges are: the output images keep the geometry of the scans.

**Check a whole batch first:**

```bash
homog calibrate -g gcp/ scans/*.jpg  # writes calib/
```

`calib/contact_sheet.png` shows, for every sheet, the original, the pixels taken for paper and the
result; `calib/calibrate.txt` flags the sheets whose paper looks odd; `calib/homog.toml` is a
ready-to-use configuration.

![Calibration: a map fragment on a white board, without and with an imagette](docs/img/calibrate-imagette.jpg)

**When the paper is not found: imagettes.** The automatic estimate assumes that bare paper is the
lightest material of the map. When this fails (map pasted on a board, white body colour, dark
paper…), crop a small sample (at least 10 × 10 px) from a scan with any image viewer and drop it in
a directory:

| Directory | Put there | Example |
| --- | --- | --- |
| `hints/paper/` | crops of bare paper | the paper of a map pasted on a board |
| `hints/keep/` | crops of colours that must not become white | a pale yellow tint close to the paper colour |

```bash
homog --paper hints/paper --keep hints/keep scans/*.jpg
```

**Lighting.** The paper surface is flexible enough for usual scans. For very even or very uneven
lighting, use `--lighting even` or `--lighting uneven`; `--lighting none` uses a single paper colour
per sheet.

### Settings file

All settings can go in a `homog.toml`, every key optional, each overridden by the command line.
Relative paths are relative to the file.

```toml
target = [255, 255, 255]     # target paper colour (sRGB)
lighting = "normal"          # even | normal | uneven | none
icc = true                   # honour embedded ICC profiles
[hints]
paper = "hints/paper"
keep = "hints/keep"
[seams]
gcp_dir = "gcp"
[output]
dir = "output"
format = "auto"              # auto (same as input) | tif | png | jpg
jpeg_quality = 95
compress = "DEFLATE"         # TIFF compression
```

### All options

```text
homog [OPTIONS] <FILES>...              correct the sheets
homog calibrate [OPTIONS] <FILES>...    analyse them, write a report and a configuration

  -c, --config <FILE>      settings file (homog.toml)
  -t, --target <r,g,b>     target paper colour [default: 255,255,255]
  -l, --lighting <L>       even | normal | uneven | none [default: normal]
      --paper <DIR>        imagettes of bare paper
      --keep <DIR>         imagettes of colours that must not become white
  -g, --gcp-dir <DIR>      GCP files (<image file name>.points) giving the layout of the sheets
      --ignore-icc         ignore embedded ICC profiles (assume sRGB)
  -j, --jobs <N>           number of threads [default: all cores]
  -v, --verbose            per-image details and timings

  correction only:
  -o, --out-dir <DIR>      output directory [default: output]
  -f, --format <F>         auto | tif | png | jpg [default: auto]
      --jpeg-quality <Q>   [default: 95]
      --compress <C>       TIFF compression: NONE, LZW, DEFLATE, ZSTD... [default: DEFLATE]
      --debug              write debug images to <out-dir>/debug/<image>/

  calibration only:
  -o, --out-dir <DIR>      report directory [default: calib]
```

Inputs: RGB or RGBA, 8 or 16 bits, any format GDAL reads. Outputs keep the bit depth, the alpha
band and the georeferencing (geotransform, CRS, GCPs); TIFF outputs are tiled and compressed.

### Debug images

With `--debug`, each sheet gets a directory `<out-dir>/debug/<image>/`:

| File | Content |
| --- | --- |
| `1_original.png` | input, reduced (~1600 px) |
| `2_paper_surface.png` | estimated paper colour |
| `3_paper_pixels.png` | pixels used to estimate it: green = bare paper, orange = washes (lightness only), grey = ignored, dark red = outside the sheet |
| `4_gain_surface.png` | brightening of the paper correction (black = none, light yellow = maximum) |
| `5_gain_seams.png` | extra correction along the seams (blue = darker, red = lighter) |
| `6_seam_samples.png` | result with the seam measurement points: green = paper, red = rejected |
| `7_result.png` | result, reduced |
| `debug.txt` | the numbers behind the images |

![Debug images: pixels used, paper surface, gain](docs/img/debug-surface.jpg)

## How it works

1. **Find the sheet** on the scan, leaving out the scanner bed and colour charts.
2. **Find the paper colour**: the dominant colour of the lightest pixels (or the imagettes).
3. **Fit a smooth paper surface** that follows the lightest pixels, i.e. the paper, under uneven
   lighting, without sinking into the washes.
4. **White-balance** every pixel so that the local paper becomes white (Bradford chromatic adaptation).
5. **Match the seams**: along the edges shared with adjacent sheets, measure the remaining paper
   colour and correct it, fading inwards.
6. **Apply** all this in one streaming pass over the full-resolution image.

The [technical report](docs/technical-report.md) explains every step, each with a plain-language
summary, the model, the results on the example atlas and what was tried and abandoned.

## Development

```bash
cargo test --release                  # unit tests and end-to-end tests on synthetic sheets
cargo clippy --release --all-targets
```

The code is in `src/`: `background.rs` (sheet, paper model, paper surface), `transform.rs`
(adaptation), `layout.rs` and `seams.rs` (seams), `io.rs` (GDAL), `hints.rs` (imagettes),
`calibrate.rs`, `debug.rs`, `config.rs`, `main.rs` (command line).

## Licence

The code is licensed under the [GNU Affero General Public License v3.0](LICENSE) or later. The
example scans belong to the Service historique de la Défense (Vincennes).
