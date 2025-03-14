# Super naive colour homogenisation with Imagemagick

A bash script based on Imagemagick for background colour homogenisation of a series of images while preserving colour differences in each image.
Designed as a pre-processing for georeferencing old maps consisting of several scanned sheets.

## Motivations

Many large historical maps are composed of several individual sheets.
Because they are digitised separately, differences in lighting conditions and variations in the colour of the paper cause the background colour of the images to vary, resulting in an ugly patchwork effect when they are stitched to form the complete map.
As the paper tends to yellow and the colours fade over time, it may also be desirable to lighten and brighten the background evenly across all the sheets so that the full map can be used in historical GIS.

## How it works

The main idea is to transform the colors of an image so that its background matches a target color, while preserving the differences between colors within the image as much as possible.

The process consists of three stages, summarised in the diagram below:

![Pipeline](figs/process.png)

Given an input image $I$ and a target background colour $t$:

1. $I$ is quantised to extract a set $C$ of dominant colours, then $`c^* = \arg\min_{c \in C} \| c - t \|_2`$ is chosen as the background colour of the image;
2. Colours of the image $I$ are rescaled by a factor $`(\frac{t_{l}}{c^*_{l}},\frac{t_{a}}{c^*_{a}},\frac{t_{b}}{c^*_{b}})`$ ;
3. A sigmoïdal contrast enhancement is eventually applied to produce $I_h$ the final homogenized image. See [https://imagemagick.org/script/command-line-options.php?#sigmoidal-contrast](https://imagemagick.org/script/command-line-options.php?#sigmoidal-contrast) for more info.

Notes:

- calculations are made in the CIELAB colour space. The difference between two colours is the Euclidean distance between them: $`d_E = \| .\|_2`$.
- Step 1 assumes that the background colour is dominant in the image. If not, a background colour can be set with the `--background-reference` option.

## Usage

[ImageMagick](https://imagemagick.org/script/download.php) v6.x or higher is required.

Tool usage:

```bash
Usage: ./equalisation.sh [options] file1 [file2 ...]

Required:
  -t, --target <r,g,b>                Target background color (default: 255,255,255)

Optional:
  -v, --verbose                       Show more information about the process
  -c, --contrast <value>              Contrast enhancement value (default: 0). See ImageMagick's 'Sigmoidal Non-linearity Contrast'.
  -b, --background-reference <r,g,b>  Reference background color (default: 255,255,255). You typically don't want to set this parameter manually unless you really know what you are doing.
  -h, --help                          Display this help message
```

Note: equalised images will be written to `./output/`

## Example

```bash
./equalisation.sh -c 5,70% -t 255,255,255 -v example/*.jpg
```

## Additional thoughts

1. CIELAB aims to be perceptually uniform, so scaling a color should ideally result in a similar perceptual shift. However, CIELAB isn't perfectly uniform, and [some areas show distortions](https://web.archive.org/web/20080705054344/http://www.aim-dtp.net/aim/evaluation/cie_de/index.htm). Although there are complex colour difference formulas to compensate, we use Euclidean distance because it is easy to compute and because historical maps often have washed-out, brownish colours which we hope are less affected by CIELAB distortions (notably affecting very dark and very saturated colours).
2. ImageMagick clamps scaled colours that fall outside the sRGB gamut. This may result in colours being 'squashed' against the gamut edges when the scaling factor is large, especially if $c^∗$ and $t$ are far apart.
