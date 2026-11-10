#!/usr/bin/env python3
"""Basic image metrics of the paired STAY vs standard H&E images, over the whole picture, at one pixel size.

For each paired section, reads data/raw_data/paired/<sample_id>/images/he_image/*.tif, rescales it to TARGET_UM
microns per pixel (the coarsest of the images, so images are only ever shrunk, by pixel averaging), and computes
over all pixels:
  - intensity: brightness, RMS contrast, dynamic range (1st-99th percentile), entropy of the grey histogram
  - sharpness: variance of the Laplacian, Tenengrad (mean squared Sobel gradient)
  - texture: GLCM contrast, homogeneity, energy and correlation (grey levels binned to GLCM_LEVELS, distance 1,
    averaged over 4 directions)
  - colour: mean red, green, blue and saturation

The standard images are too large to load whole (up to ~70 GB), so every image is read, rescaled and measured in
horizontal bands and the metrics are summed over bands (GLCM pixel pairs across band borders are not counted).

Writes results/he_quality/he_quality_sections.csv (one row per section); plot_he_quality.py compares and plots it.

  python code/figures/he_quality/he_quality.py
"""

from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import tifffile
import zarr
from skimage.feature import graycomatrix, graycoprops

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RAW_DIR = PROJECT_ROOT / "data" / "raw_data" / "paired"
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"
SAMPLES = pd.read_csv(PROJECT_ROOT / "resources" / "sample_list.csv", dtype=str).query("cohort == 'paired'")

# pixel size of each H&E image, from its TIFF tags
METAFER_200300 = 3450 / 200300000 * 1e4  # XResolution 200300000/3450 px per cm
METAFER_200400 = 3450 / 200400000 * 1e4  # XResolution 200400000/3450 px per cm
QUPATH_EXPORT = 1e6 / 2259600  # XResolution 2259600/1000000 px per unit, ImageJ unit=um
MICRONS_PER_PIXEL = {
    "Patient_1_STAY": METAFER_200300,
    "Patient_1_standard": METAFER_200300,
    "Patient_2_STAY": QUPATH_EXPORT,
    "Patient_2_standard": METAFER_200400,
    "Patient_3_STAY": QUPATH_EXPORT,
    "Patient_3_standard": METAFER_200400,
    "Patient_4_STAY": QUPATH_EXPORT,
    "Patient_4_standard": METAFER_200300,
}
TARGET_UM = max(MICRONS_PER_PIXEL.values())
BAND_ROWS = 1024  # output rows per band
GLCM_LEVELS = 32


def open_he(sample_id: str) -> tuple[zarr.Array, str]:
    """Full-resolution H&E image as a lazily read (rows, cols, RGB) array."""
    path = sorted((RAW_DIR / sample_id / "images" / "he_image").glob("*.tif*"))[0]
    arr = zarr.open(tifffile.imread(path, aszarr=True), mode="r")
    if isinstance(arr, zarr.Group):  # pyramidal file: level 0 is full resolution
        arr = arr["0"]
    if arr.ndim != 3 or arr.shape[-1] < 3 or arr.dtype != np.uint8:
        raise SystemExit(f"{sample_id}: expected an 8-bit RGB image, got {arr.shape} {arr.dtype}")
    return arr, path.name


def rescaled_bands(arr: zarr.Array, scale: float):
    """Yield the image rescaled by scale (<= 1), BAND_ROWS output rows at a time."""
    height, width = arr.shape[:2]
    if arr.chunks[0] >= height:  # stored as one block (the uncompressed QuPath exports): read it once
        arr = arr[:]
    out_h, out_w = round(height * scale), round(width * scale)
    for o0 in range(0, out_h, BAND_ROWS):
        o1 = min(o0 + BAND_ROWS, out_h)
        src = np.asarray(arr[round(o0 / scale) : min(height, round(o1 / scale)), :, :3])
        yield src if scale == 1 else cv2.resize(src, (out_w, o1 - o0), interpolation=cv2.INTER_AREA)


def image_metrics(bands) -> dict:
    hist = np.zeros(256)
    glcm = np.zeros((GLCM_LEVELS, GLCM_LEVELS, 1, 4))
    rgb = np.zeros(3)
    lap_sum = lap_sq = grad_sq = sat = 0.0
    n = n_filtered = height = width = 0
    carry = None
    for band in bands:
        gray = cv2.cvtColor(band, cv2.COLOR_RGB2GRAY)
        hist += cv2.calcHist([gray], [0], None, [256], [0, 256]).ravel()
        # 3x3 filters need the rows above and below: prepend the previous band's last two rows and keep only
        # rows with both neighbours, so every row except the image's first and last is filtered exactly once
        ext = gray if carry is None else np.vstack([carry, gray])
        carry = gray[-2:]
        lap = cv2.Laplacian(ext, cv2.CV_64F)[1:-1]
        lap_sum += lap.sum()
        lap_sq += np.square(lap).sum()
        grad_sq += np.square(cv2.Sobel(ext, cv2.CV_64F, 1, 0)[1:-1]).sum() + np.square(cv2.Sobel(ext, cv2.CV_64F, 0, 1)[1:-1]).sum()
        n_filtered += lap.size
        glcm += graycomatrix(
            (gray // (256 // GLCM_LEVELS)).astype(np.uint8),
            distances=[1],
            angles=[0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
            levels=GLCM_LEVELS,
            symmetric=True,
        )
        rgb += band.reshape(-1, 3).sum(axis=0, dtype=np.float64)
        sat += cv2.cvtColor(band, cv2.COLOR_RGB2HSV)[..., 1].sum(dtype=np.float64)
        n += gray.size
        height += band.shape[0]
        width = band.shape[1]

    p = hist / n
    levels = np.arange(256)
    mean = (levels * p).sum()
    cdf = np.cumsum(p)
    glcm /= glcm.sum(axis=(0, 1), keepdims=True)
    return {
        "height_px": height,
        "width_px": width,
        "brightness": float(mean),
        "rms_contrast": float(np.sqrt((levels**2 * p).sum() - mean**2)),
        "dynamic_range": float(np.searchsorted(cdf, 0.99) - np.searchsorted(cdf, 0.01)),
        "entropy": float(-(p[p > 0] * np.log2(p[p > 0])).sum()),
        "laplacian_var": float(lap_sq / n_filtered - (lap_sum / n_filtered) ** 2),
        "tenengrad": float(grad_sq / n_filtered),
        **{f"glcm_{prop}": float(graycoprops(glcm, prop).mean()) for prop in ["contrast", "homogeneity", "energy", "correlation"]},
        "mean_red": float(rgb[0] / n),
        "mean_green": float(rgb[1] / n),
        "mean_blue": float(rgb[2] / n),
        "saturation": float(sat / n),
    }


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for _, row in SAMPLES.iterrows():
        sample_id = str(row["sample_id"])
        arr, name = open_he(sample_id)
        um = MICRONS_PER_PIXEL[sample_id]
        print(f"{sample_id}: {name} {arr.shape[1]}x{arr.shape[0]} px at {um:.4f} um/px -> {TARGET_UM:.4f} um/px", flush=True)
        metrics = image_metrics(rescaled_bands(arr, um / TARGET_UM))
        rows.append({"sample_id": sample_id, "patient_id": str(row["patient_id"]), "method": sample_id.rsplit("_", 1)[1],
                     "file": name, "um_per_px_original": um, "um_per_px": TARGET_UM, **metrics})
    out = RESULTS_DIR / "he_quality_sections.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
