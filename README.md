# Super naive colour homogenisation

A fast command-line tool (`homog`, written in Rust) for background colour homogenisation of a series of images while preserving colour differences in each image.
Designed as a pre-processing for georeferencing old maps consisting of several scanned sheets.

## Motivations

Many large historical maps are composed of several individual sheets.
Because they are digitised separately, differences in lighting conditions and variations in the colour of the paper cause the background colour of the images to vary, resulting in an ugly patchwork effect when they are stitched to form the complete map.
As the paper tends to yellow and the colours fade over time, it may also be desirable to lighten and brighten the background evenly across all the sheets so that the full map can be used in historical GIS.

## How it works

The paper is treated as the illuminant of the scene: each sheet is "white-balanced" so that its paper becomes the target colour $t$ (white by default), using a von Kries chromatic adaptation in the Bradford cone space.
Since the paper colour also varies *within* a sheet (uneven lighting, vignetting, yellowed edges), it is estimated locally (flat-field).

Given an input image $I$:

1. **Paper estimation** on a ~400 px thumbnail, in CIELAB:
   - content mask: the scanner bed is the uniform colour touching the image border; it is flood-filled from the border and only the largest remaining region (the sheet) is kept, which drops colour charts or rulers lying on the bed;
   - global paper colour $p$: median of the bright pixels (75th–98th percentile of $L^*$) whose chroma is close to the dominant one — or, if `paper` imagettes are given, of the pixels matching them. Pixels matching `keep` imagettes are never taken for paper;
   - paper surface $p(x, y)$: smooth polynomial surfaces for $L^*$, $a^*$ and $b^*$, fitted on the pixels compatible with the paper model (Mahalanobis distance in $a^*b^*$, the tolerance following the spread of the paper colour). Bare paper is the lightest material, since a transparent wash or ink can only darken it, so $L^*$ is fitted as an *upper envelope* (80 % quantile regression), and $a^*, b^*$ only on the pixels close to that envelope, with a robust (Tukey) rejection of other hues. A low-degree surface follows uneven lighting and yellowed edges but cannot follow a wash patch, however large: pale washes keep their colour instead of being bleached. The degree (1 to 5) is chosen per sheet by spatial cross-validation of the envelope fit: the smallest degree within 1 % of the best held-out loss.
2. **Adaptation**, per pixel, in linear light:
   $`\mathrm{RGB}' = M^{-1}\,\mathrm{diag}\!\left(\frac{M\,t_{XYZ}}{M\,p_{XYZ}(x,y)}\right) M\,\mathrm{RGB}`$, with $M = M_{\mathrm{Bradford}} \cdot M_{\mathrm{sRGB \to XYZ}}$.

   The paper becomes $t$, black stays black, greys of the paper's tint become neutral, and the other colours are adapted consistently (no hue shift, no clipping of pale washes darker than the paper).

3. **Seams** (optional, `--gcp-dir`): the paper surface cannot follow the darkening close to the edges of a sheet, which shows as a step where two sheets meet. The layout of the sheets is read from their ground control points (QGIS georeferencer `.points` files): a map ↔ pixel affine transform per sheet, and the sheet extent, taken as the bounding box of its GCPs (which sit on the neatline). Along every edge shared with another sheet of the GCP directory, the paper is measured on the corrected image, in a band 1 % of the sheet size wide starting beyond twice the RMS residual of the GCPs (at least 0.3 %), per segment (24 along the edge), then brought back to $t$ by an extra Bradford gain. The gain is full up to the middle of the band and fades out where the measured darkening profile has decayed to a quarter (4–20 % of the sheet size). **Images are not resampled**: the GCPs only tell which sheets touch and where.

The gains are computed on the low-resolution surface and bilinearly interpolated, so the full-resolution pass is a single streaming loop (decode → 2 matrix products → encode), parallelised over rows.

## Installation

Requirements: Rust ≥ 1.85, GDAL 3.x (tested with 3.12) with its development files and `libclang` (for the bindings). On Debian/Ubuntu:

```bash
sudo apt install libgdal-dev libclang-dev
cargo build --release   # binary in target/release/homog
```

## Usage

```text
homog [OPTIONS] <FILES>...              correct the sheets
homog calibrate [OPTIONS] <FILES>...    analyse them and propose a configuration

Common options:
  -c, --config <FILE>           Configuration file (homog.toml)
  -t, --target <r,g,b>          Target background colour [default: 255,255,255]
  -l, --lighting <L>            auto | even | uneven | none [default: auto]
      --paper <DIR>             Imagettes of bare paper
      --keep <DIR>              Imagettes of colours that must not become white
  -g, --gcp-dir <DIR>           GCP files (<image file name>.points, QGIS georeferencer) giving the
                                sheet layout: seams with adjacent sheets are corrected
      --ignore-icc              Ignore embedded ICC profiles (assume sRGB)
  -j, --jobs <N>                Number of threads [default: all cores]
  -v, --verbose                 Per-image details and timings

Correction only:
  -o, --out-dir <DIR>           Output directory [default: output]
  -f, --format <FORMAT>         auto | tif | png | jpg [default: auto = same family as the input]
      --jpeg-quality <Q>        [default: 95]
      --compress <C>            TIFF compression: NONE, LZW, DEFLATE, ZSTD... [default: DEFLATE]
      --debug                   Write debug images to <out-dir>/debug/<image>/

Calibration only:
  -o, --out-dir <DIR>           Report directory [default: calib]
```

Examples:

```bash
homog -v example/*.jpg                                         # single sheets
homog -v -g example/GCP --debug example/*__00{17,19,21}.jpg    # adjacent sheets, with seams
homog calibrate -g example/GCP example/*.jpg                   # batch analysis
homog -c calib/homog.toml example/*.jpg                        # correction with the proposed settings
```

The sheets that touch are found among all the `.points` files of the GCP directory, so a sheet gets the same correction whether its neighbours are processed in the same run or not.

## Settings

Everything that can be measured is automatic: where the sheet is on the scan, the paper colour and its tolerance, the flexibility of the paper surface, the position and reach of the seam correction. What only the user knows is given as **imagettes**: small crops (at least 10 × 10 px, PNG/JPEG/TIFF) cut from the scans with any image viewer, one example per file.

| Directory | Content | When |
| --- | --- | --- |
| `paper` | bare paper | the automatic estimate picks something else: map pasted on a board, body colour (gouache) lighter than the paper, paper darker than large flat tints |
| `keep` | colours that must not become white | a pale wash or tint close to the paper colour and taken for it |

Each imagette gives a colour model (robust centre and spread); several `paper` imagettes describe several papers (e.g. clean and stained). The imagettes apply to all the sheets of the project.

Settings can be gathered in a `homog.toml` (relative paths are relative to this file), each overridden by the corresponding option:

```toml
target = [255, 255, 255]
lighting = "auto"            # auto | even (degree 2) | uneven (degree 4) | none (single colour)
icc = true
[hints]
paper = "hints/paper"
keep = "hints/keep"
[seams]
gcp_dir = "GCP"
[output]
dir = "output"
format = "auto"              # auto | tif | png | jpg
jpeg_quality = 95
compress = "DEFLATE"
```

### Calibrating a batch

`homog calibrate` analyses the sheets without writing corrected images and writes to `calib/`:

- `contact_sheet.png`: one row per sheet: original | pixels taken for paper (green = bare paper, orange = lightness only, grey = ignored, dark red = outside the sheet) | result;
- `calibrate.txt`: batch paper colour and, per sheet, paper colour, distance to the batch, share of bare paper found, surface degree, seams and warnings (*paper far from the batch*, *little bare paper found*);
- `homog.toml`: proposed configuration.

Typical workflow: run `calibrate`, look at the rows of the sheets with warnings; if something other than paper is green, crop it into `hints/keep/`; if the paper of the map is not green, crop it into `hints/paper/`; rerun `calibrate` with the hints, then `homog -c calib/homog.toml`.

### Debug images

With `--debug`, each sheet gets a `<out-dir>/debug/<image>/` directory:

| File | Content |
| --- | --- |
| `1_original.png` | input (reduced, ~1600 px) |
| `2_paper_surface.png` | estimated paper colour (flat-field) |
| `3_paper_pixels.png` | pixels used by the surface fit: green = bare paper ($L^*$ and $a^*b^*$), orange = $L^*$ envelope only, grey = ignored (ink, strong washes, `keep` colours), dark red = outside the sheet |
| `4_gain_surface.png` | luminance gain of the surface correction (black = ×1, light yellow = maximum, given in `debug.txt`) |
| `5_gain_seams.png` | luminance gain of the seam correction (blue = darker, white = unchanged, red = lighter; symmetric scale given in `debug.txt`) |
| `6_seam_samples.png` | result with the seam measurement points: green = paper, red = rejected |
| `7_result.png` | result (reduced) |
| `debug.txt` | global paper colour and spread, surface degree, gain ranges, seam band, reach and gain per segment |

Supported inputs: RGB or RGBA images, 8 or 16 bits (JPEG, TIFF/GeoTIFF, PNG, and any format GDAL reads).

- **Georeferencing** (geotransform, CRS, GCPs), nodata and colour interpretation are copied to the output. GeoTIFF outputs are tiled (512×512) and compressed.
- **ICC profiles**: an embedded profile is honoured (conversion through Little CMS); the output is sRGB.
- **Bit depth** is preserved (16-bit inputs cannot be written as JPEG).

## Performance

5800 × 4000 JPEG sheet (`example/…0019.jpg`), 20 threads:

| Input | Time |
| --- | --- |
| Single sheet | 0.5 s |
| With seams (`-g example/GCP`) | 0.6 s |

Most of the time is spent in JPEG decoding/encoding; the analysis (content mask, paper model, degree selection, surface fit) takes ~0.15 s per sheet whatever its size, the seam correction ~0.1 s more (a second, reduced read of the image). The 59 GCP-referenced sheets of the example atlas are calibrated in 8 s.

## Tests and prototype

```bash
cargo test --release          # unit tests + end-to-end tests on synthetic sheets
uv run proto/homog_proto.py -g example/GCP example/*__00{17,19,21}.jpg   # Python prototype
```

`proto/homog_proto.py` is the prototype used to design the colour model (Bradford adaptation, paper surface, seams; it writes `proto/out/comparison.png`, `proto/out/seams.png` and metrics; `-s` sets the processing scale, 0.25 by default). It does not include the automatic settings nor the imagettes of the Rust tool. `proto/compare_de.py` reports ΔE statistics between two renderings.

## Additional thoughts

1. The adaptation is exact for the paper and for neutral colours; for saturated inks, von Kries/Bradford is an approximation of how they would look on white paper (~2 ΔE in the end-to-end test).
2. Colours lighter than the local paper (rare on old maps) are clipped to the target.
3. The seam correction makes the paper continuous across sheets. Washes and inks may still differ slightly (hand-painted sheets; ~1–2.5 ΔE across the seams of the examples): a per-sheet colour matrix fitted on colour pairs across the seams was tried and brought no measurable improvement.
4. The extent of a sheet is the bounding box of its GCPs: sheets whose GCPs are not placed on the neatline get an approximate extent, hence no or partial seam correction.
5. Imagettes apply to the whole project: a `paper` imagette taken on one sheet is also used for the others, which is what is wanted when the sheets share the same paper.
