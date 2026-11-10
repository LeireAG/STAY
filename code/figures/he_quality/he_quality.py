#!/usr/bin/env python3
"""Basic image metrics of the paired STAY vs standard H&E images, over the whole picture.

For each paired section, reads data/raw_data/paired/<sample_id>/images/he_image/*.tif at its native resolution
and computes, over all pixels (no tissue detection, no regions, no TIFF metadata):
  - intensity: brightness, RMS contrast, dynamic range (1st-99th percentile), entropy of the grey histogram
  - sharpness: variance of the Laplacian, Tenengrad (mean squared Sobel gradient)
  - texture: GLCM contrast, homogeneity, energy and correlation (grey levels binned to GLCM_LEVELS, distance 1,
    averaged over 4 directions)
  - colour: mean red, green, blue and saturation

The whole picture includes the glass around the tissue, and images are compared at whatever resolution they were
scanned at, so check height_px/width_px in he_quality_sections.csv before reading sharpness and texture.
There are four patients with one section per protocol, so no p-values: the result is the STAY - standard
difference per patient, and for metrics with a better direction, which protocol is better in how many patients.

Writes to results/he_quality:
  - he_quality_sections.csv: one row per section
  - he_quality_paired.csv: per metric and patient, STAY - standard and which protocol is better
  - he_quality.pdf/png: per-patient paired values of every metric

  python code/figures/he_quality/he_quality.py
"""

from __future__ import annotations

from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from skimage.feature import graycomatrix, graycoprops

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RAW_DIR = PROJECT_ROOT / "data" / "raw_data" / "paired"
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"
SAMPLES = pd.read_csv(PROJECT_ROOT / "resources" / "sample_list.csv", dtype=str).query("cohort == 'paired'")

GLCM_LEVELS = 32
# per patient, a difference smaller than this fraction of the pair's mean counts as a tie
TIE = 0.05

METHODS = ["standard", "STAY"]
# metric: (label, better direction: 1 higher is better, -1 lower is better, 0 no better direction)
METRICS = {
    "brightness": ("Brightness", 0),
    "rms_contrast": ("RMS contrast", 1),
    "dynamic_range": ("Dynamic range (p1-p99)", 1),
    "entropy": ("Entropy (bits)", 1),
    "laplacian_var": ("Sharpness: Laplacian variance", 1),
    "tenengrad": ("Sharpness: Tenengrad", 1),
    "glcm_contrast": ("GLCM contrast", 0),
    "glcm_homogeneity": ("GLCM homogeneity", 0),
    "glcm_energy": ("GLCM energy", 0),
    "glcm_correlation": ("GLCM correlation", 0),
    "mean_red": ("Mean red", 0),
    "mean_green": ("Mean green", 0),
    "mean_blue": ("Mean blue", 0),
    "saturation": ("Mean saturation", 0),
}


def read_he(sample_id: str) -> tuple[np.ndarray, str]:
    path = sorted((RAW_DIR / sample_id / "images" / "he_image").glob("*.tif*"))[0]
    img = np.squeeze(tifffile.imread(path))
    if img.ndim == 3 and img.shape[0] in (3, 4) and img.shape[-1] not in (3, 4):
        img = np.moveaxis(img, 0, -1)
    if img.ndim != 3:
        raise SystemExit(f"{sample_id}: expected an RGB image, got shape {img.shape}")
    img = np.ascontiguousarray(img[..., :3])
    if img.dtype != np.uint8:
        img = (img.astype(np.float32) / img.max() * 255).astype(np.uint8)
    return img, path.name


def image_metrics(img: np.ndarray) -> dict:
    gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    hist = np.bincount(gray.ravel(), minlength=256) / gray.size
    p1, p99 = np.percentile(gray, [1, 99])
    gx = cv2.Sobel(gray, cv2.CV_32F, 1, 0)
    gy = cv2.Sobel(gray, cv2.CV_32F, 0, 1)
    tenengrad = float(np.mean(np.square(gx, out=gx), dtype=np.float64) + np.mean(np.square(gy, out=gy), dtype=np.float64))
    del gx, gy
    glcm = graycomatrix(
        (gray // (256 // GLCM_LEVELS)).astype(np.uint8),
        distances=[1],
        angles=[0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
        levels=GLCM_LEVELS,
        symmetric=True,
        normed=True,
    )
    return {
        "height_px": img.shape[0],
        "width_px": img.shape[1],
        "brightness": float(gray.mean()),
        "rms_contrast": float(gray.std()),
        "dynamic_range": float(p99 - p1),
        "entropy": float(-(hist[hist > 0] * np.log2(hist[hist > 0])).sum()),
        "laplacian_var": float(cv2.Laplacian(gray, cv2.CV_32F).var(dtype=np.float64)),
        "tenengrad": tenengrad,
        **{f"glcm_{p}": float(graycoprops(glcm, p).mean()) for p in ["contrast", "homogeneity", "energy", "correlation"]},
        "mean_red": float(img[..., 0].mean()),
        "mean_green": float(img[..., 1].mean()),
        "mean_blue": float(img[..., 2].mean()),
        "saturation": float(cv2.cvtColor(img, cv2.COLOR_RGB2HSV)[..., 1].mean()),
    }


def compare_pairs(sections: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric, (_, better) in METRICS.items():
        wide = sections.pivot(index="patient_id", columns="method", values=metric)
        for patient, row in wide.iterrows():
            diff = row["STAY"] - row["standard"]
            if better == 0 or not np.isfinite(diff):
                winner = "-"
            elif abs(diff) <= TIE * abs(row[METHODS].mean()):
                winner = "tie"
            else:
                winner = "STAY" if np.sign(diff) == better else "standard"
            rows.append({"metric": metric, "patient_id": patient, "standard": row["standard"], "STAY": row["STAY"],
                         "STAY_minus_standard": diff, "better": winner})
    return pd.DataFrame(rows)


def verdict(better: pd.Series) -> str:
    """How many patients each protocol is better in, e.g. "STAY 3, standard 0, tie 1"."""
    if (better == "-").all():
        return "no better direction"
    return ", ".join(f"{k} {(better == k).sum()}" for k in ["STAY", "standard", "tie"])


def plot_sections(sections: pd.DataFrame, pairs: pd.DataFrame) -> None:
    ncol = 5
    nrow = int(np.ceil(len(METRICS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.2 * ncol, 2.3 * nrow), squeeze=False)
    patients = sorted(sections["patient_id"].unique())
    colors = dict(zip(patients, plt.get_cmap("tab10").colors))
    for ax, (metric, (label, _)) in zip(axes.ravel(), METRICS.items()):
        wide = sections.pivot(index="patient_id", columns="method", values=metric)
        for patient, row in wide.iterrows():
            ax.plot([0, 1], row[METHODS], "-o", color=colors[patient], ms=4, lw=1, label=patient)
        ax.set_title(f"{label}\n{verdict(pairs.loc[pairs['metric'] == metric, 'better'])}")
        ax.set_xticks([0, 1], METHODS)
        ax.set_xlim(-0.3, 1.3)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.ravel()[len(METRICS) :]:
        ax.axis("off")
    axes[0, 0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_quality.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "pdf.fonttype": 42})
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for _, row in SAMPLES.iterrows():
        sample_id = str(row["sample_id"])
        img, name = read_he(sample_id)
        print(f"{sample_id}: {name} {img.shape[1]}x{img.shape[0]} px", flush=True)
        rows.append({"sample_id": sample_id, "patient_id": str(row["patient_id"]),
                     "method": sample_id.rsplit("_", 1)[1], "file": name, **image_metrics(img)})
        del img
    sections = pd.DataFrame(rows)
    pairs = compare_pairs(sections)

    sections.to_csv(RESULTS_DIR / "he_quality_sections.csv", index=False)
    pairs.to_csv(RESULTS_DIR / "he_quality_paired.csv", index=False)
    plot_sections(sections, pairs)

    summary = pairs.groupby("metric", sort=False).agg(
        better_in=("better", verdict),
        median_STAY_minus_standard=("STAY_minus_standard", "median"),
    )
    print(summary.to_string(), flush=True)
    print(f"wrote {RESULTS_DIR}", flush=True)


if __name__ == "__main__":
    main()
