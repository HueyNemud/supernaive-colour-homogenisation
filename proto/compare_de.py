# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "scipy", "pillow"]
# ///
"""ΔE76 statistics between two renderings of the same image (e.g. Rust vs prototype).

Usage: uv run proto/compare_de.py a.png b.png [c.png d.png ...]
"""

import sys

import numpy as np
from PIL import Image

from homog_proto import srgb_to_lab


def load(path):
    im = Image.open(path)
    a = np.asarray(im)
    scale = 65535.0 if a.dtype == np.uint16 else 255.0
    return a[..., :3].astype(np.float64) / scale


def main():
    paths = sys.argv[1:]
    if len(paths) < 2 or len(paths) % 2:
        sys.exit(__doc__)
    for a, b in zip(paths[::2], paths[1::2]):
        de = np.linalg.norm(srgb_to_lab(load(a)) - srgb_to_lab(load(b)), axis=-1)
        print(f"{a} vs {b}: ΔE mean {de.mean():.3f}, p99 {np.percentile(de, 99):.3f}, max {de.max():.3f}")


if __name__ == "__main__":
    main()
