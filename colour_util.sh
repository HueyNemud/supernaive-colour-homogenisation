#!/bin/bash
set -e

# The name of the convert command is passed as an argument. Default is "convert".
convert=${1:-convert}

function asvars() {
   # Replaces commas with spaces to facilitate Bash variable assignment
   echo "${1//,/ }"
}

function top_n_colors {
  local image=$1
  local n_colors=$2

  # 1. Force sRGB input to ensure correct math
  # 2. Apply blur and resize to speed up processing and average noise
  # 3. Reduce the palette to N colors
  # 4. Convert to CIELab and output the histogram
  $convert "$image" -set colorspace sRGB -blur 0x8 -resize 100x100 \
    +dither -colors "$n_colors" -colorspace CIELab \
    -format "%c" histogram:info:- | \
    sed -n 's/.*cielab(\([^)]*\)).*/\1/p' | tr ',' ' '
}

function colordifference {
  # Compute the color difference between two colors in LAB space
  # Formula: Euclidean distance = sqrt((L1 - L2)^2 + (a1 - a2)^2 + (b1 - b2)^2)
  # Arguments: L1 a1 b1 L2 a2 b2
  echo "scale=4; sqrt(($1 - $4)^2 + ($2 - $5)^2 + ($3 - $6)^2)" | bc
}

function n_closest() {
  local ref_l=$1
  local ref_a=$2
  local ref_b=$3
  local image=$4
  local n=$5

  local best_dist=999999
  local best_color=""

  # Read the space-separated output from top_n_colors
  while read -r l a b; do
    [[ -z "$l" ]] && continue
    
    # Calculate distance
    distance=$(colordifference "$l" "$a" "$b" "$ref_l" "$ref_a" "$ref_b")
    
    # Check if this is the closest color so far
    if (( $(echo "$distance < $best_dist" | bc -l) )); then
      best_dist=$distance
      best_color="$l $a $b"
    fi
  done < <(top_n_colors "$image" "$n")

  # Output ONLY the best color values for the calling script to capture
  echo "$best_color"
}

function rgb_norm {
  # Normalize RGB values (0-255) to 0-1 range
  for param in "$@"; do
    echo "scale=10; $param / 255" | bc
  done
}

function rgb2lab {
  # Convert sRGB (0-255) to CIELab using ImageMagick
  # L: 0 to 100 | a & b: -128 to 127
  local srgb="srgb($1,$2,$3)"
  local conversion_formula='%[fx:100*u.r],%[fx:255*(u.g-0.5)],%[fx:255*(u.b-0.5)]'
  
  local lab=$($convert xc:"$srgb" -colorspace CIELab -format "$conversion_formula" info:)
  echo "${lab//,/ }"
}

function lab2rgb {
  # Convert CIELab back to sRGB (0-255)
  local lab="cielab($1,$2,$3)"
  local conversion_formula='%[fx:round(255*u.r)],%[fx:round(255*u.g)],%[fx:round(255*u.b)]'
  
  local rgb=$($convert xc:"$lab" -colorspace sRGB -format "$conversion_formula" info:)
  echo "${rgb//,/ }"
}

function lab_denorm {
  # Denormalize a LAB color from 0-1 range to standard Lab scales
  l=$(echo "scale=10; $1 * 100" | bc)
  a=$(echo "scale=10; 255 * ($2 - 0.5)" | bc)
  b=$(echo "scale=10; 255 * ($3 - 0.5)" | bc)
  echo "$l $a $b"
}

function lab_norm {
  # Normalize a LAB color from standard scales to 0-1 range
  l=$(echo "scale=10; $1 / 100" | bc)
  a=$(echo "scale=10; $2 / 255 + 0.5" | bc)
  b=$(echo "scale=10; $3 / 255 + 0.5" | bc)
  echo "$l $a $b"
}