#!/usr/bin/env python3
"""Thumbnails and crops of the paired H&E images, to look at what the metrics of he_quality.py measure.

For each paired section (same images and pixel sizes as he_quality.py):
  - a thumbnail of the whole image at THUMB_UM microns per pixel
  - N_CROPS crops of CROP_UM x CROP_UM, read from the full-resolution image and rescaled to TARGET_UM (the pixel
    size the metrics use). Crop centres are random positions (fixed SEED) where the whole crop is tissue on the
    thumbnail (saturation above Otsu); set CROP_CENTRES_UM to place them by hand after looking at the thumbnails.

Writes to results/he_quality, for plot_he_overview.py:
  - thumbnails/<sample_id>.png
  - crops/<sample_id>_<n>.png
  - crops.csv: one row per crop, centre and size in microns from the image's top-left corner

  python code/figures/he_quality/he_crops.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from he_quality import MICRONS_PER_PIXEL, RESULTS_DIR, SAMPLES, TARGET_UM, open_he, rescaled_bands

THUMB_FACTOR = 16  # thumbnail pixel = 16 x 16 metric pixels
THUMB_UM = TARGET_UM * THUMB_FACTOR
CROP_UM = 300
N_CROPS = 2
SEED = 0
# sample_id: [(x_um, y_um), ...] crop centres from the image's top-left corner, replacing the random ones
CROP_CENTRES_UM: dict[str, list[tuple[float, float]]] = {}


def thumbnail(arr, um: float) -> np.ndarray:
    bands = []
    for band in rescaled_bands(arr, um / TARGET_UM):
        size = (round(band.shape[1] / THUMB_FACTOR), max(1, round(band.shape[0] / THUMB_FACTOR)))
        bands.append(cv2.resize(band, size, interpolation=cv2.INTER_AREA))
    return np.vstack(bands)


def random_centres(thumb: np.ndarray, rng: np.random.Generator) -> list[tuple[float, float]]:
    sat = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)[..., 1]
    tissue = (sat > threshold_otsu(sat)).astype(np.float32)
    window = round(CROP_UM / THUMB_UM)
    inside = ndi.uniform_filter(tissue, size=window, mode="constant") > 0.999
    ys, xs = np.nonzero(inside)
    if not len(ys):
        raise SystemExit(f"no {CROP_UM} um window is fully tissue; set CROP_CENTRES_UM by hand")
    pick = rng.choice(len(ys), size=min(N_CROPS, len(ys)), replace=False)
    return [((xs[i] + 0.5) * THUMB_UM, (ys[i] + 0.5) * THUMB_UM) for i in pick]


def crop(arr, um: float, x_um: float, y_um: float) -> np.ndarray:
    size = round(CROP_UM / um)
    x0 = min(max(round(x_um / um - size / 2), 0), arr.shape[1] - size)
    y0 = min(max(round(y_um / um - size / 2), 0), arr.shape[0] - size)
    region = np.asarray(arr[y0 : y0 + size, x0 : x0 + size, :3])
    out = round(CROP_UM / TARGET_UM)
    return region if size == out else cv2.resize(region, (out, out), interpolation=cv2.INTER_AREA)


def save_png(path: Path, rgb: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(path), cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR))


def main() -> None:
    rows = []
    for _, row in SAMPLES.iterrows():
        sample_id = str(row["sample_id"])
        arr, name = open_he(sample_id)
        um = MICRONS_PER_PIXEL[sample_id]
        print(f"{sample_id}: {name}", flush=True)
        thumb = thumbnail(arr, um)
        save_png(RESULTS_DIR / "thumbnails" / f"{sample_id}.png", thumb)
        centres = CROP_CENTRES_UM.get(sample_id) or random_centres(thumb, np.random.default_rng(SEED))
        for n, (x_um, y_um) in enumerate(centres, start=1):
            save_png(RESULTS_DIR / "crops" / f"{sample_id}_{n}.png", crop(arr, um, x_um, y_um))
            rows.append({"sample_id": sample_id, "patient_id": str(row["patient_id"]),
                         "method": sample_id.rsplit("_", 1)[1], "crop": n, "x_um": x_um, "y_um": y_um,
                         "size_um": CROP_UM, "crop_um_per_px": TARGET_UM, "thumb_um_per_px": THUMB_UM})
    out = RESULTS_DIR / "crops.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
