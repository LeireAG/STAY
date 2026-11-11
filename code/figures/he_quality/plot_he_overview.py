#!/usr/bin/env python3
"""Arrange the thumbnails and crops written by he_crops.py into figures, without touching the raw images.

Writes to results/he_quality:
  - he_overview.pdf/png: whole sections, one row per patient, standard and STAY side by side; boxes mark the crops.
    Scale bars, OVERVIEW_BAR_UM.
  - he_crops_<patient_id>.pdf/png: per crop, standard and STAY at the metrics' pixel size, each next to a zoom
    on its centre (box in the crop). Scale bars, CROP_BAR_UM on crops and ZOOM_BAR_UM on zooms.

  python code/figures/he_quality/plot_he_overview.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"

MM = 1 / 25.4
METHODS = ["standard", "STAY"]
METHOD_LABELS = {"standard": "Standard", "STAY": "STAY"}
ZOOM_UM = 75
OVERVIEW_BAR_UM = 1000
CROP_BAR_UM = 50
ZOOM_BAR_UM = 10
BOX_COLOR = "#0072B2"


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def show(ax: plt.Axes, img, um_per_px: float, bar_um: float) -> None:
    """Image with a black scale bar in the bottom-right corner."""
    ax.imshow(img, interpolation="nearest")
    h, w = img.shape[:2]
    bar = bar_um / um_per_px
    ax.add_patch(Rectangle((w * 0.95 - bar, h * 0.93), bar, h * 0.015, color="black", lw=0))
    ax.set_axis_off()


def plot_overview(crops: pd.DataFrame) -> None:
    patients = sorted(crops["patient_id"].unique())
    fig, axes = plt.subplots(len(patients), 2, figsize=(120 * MM, 45 * MM * len(patients)), squeeze=False)
    for i, patient in enumerate(patients):
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = crops.loc[(crops["patient_id"] == patient) & (crops["method"] == method)]
            um = sub["thumb_um_per_px"].iloc[0]
            show(ax, plt.imread(RESULTS_DIR / "thumbnails" / f"{sub['sample_id'].iloc[0]}.png"), um, OVERVIEW_BAR_UM)
            for _, c in sub.iterrows():
                size = c["size_um"] / um
                ax.add_patch(Rectangle((c["x_um"] / um - size / 2, c["y_um"] / um - size / 2), size, size,
                                       fill=False, ec=BOX_COLOR, lw=0.8))
                ax.text(c["x_um"] / um + size / 2, c["y_um"] / um - size / 2, str(c["crop"]), color=BOX_COLOR,
                        ha="left", va="bottom", fontsize=6, fontweight="bold")
            if i == 0:
                ax.set_title(METHOD_LABELS[method], fontsize=7)
        axes[i, 0].text(-0.04, 0.5, patient.replace("_", " "), transform=axes[i, 0].transAxes, rotation=90,
                        ha="right", va="center", fontsize=7)
    fig.tight_layout()
    save(fig, "he_overview")


def plot_crops(crops: pd.DataFrame, patient: str) -> None:
    sub = crops.loc[crops["patient_id"] == patient]
    n_crops = int(sub["crop"].max())
    fig, axes = plt.subplots(n_crops, 4, figsize=(183 * MM, 47 * MM * n_crops), squeeze=False,
                             gridspec_kw={"wspace": 0.04, "hspace": 0.08})
    for i in range(n_crops):
        for j, method in enumerate(METHODS):
            c = sub.loc[(sub["method"] == method) & (sub["crop"] == i + 1)].iloc[0]
            um = c["crop_um_per_px"]
            img = plt.imread(RESULTS_DIR / "crops" / f"{c['sample_id']}_{i + 1}.png")
            h, w = img.shape[:2]
            z = ZOOM_UM / um
            x0, y0 = round((w - z) / 2), round((h - z) / 2)
            ax_crop, ax_zoom = axes[i, 2 * j], axes[i, 2 * j + 1]
            show(ax_crop, img, um, CROP_BAR_UM)
            ax_crop.add_patch(Rectangle((x0, y0), z, z, fill=False, ec=BOX_COLOR, lw=0.8))
            show(ax_zoom, img[y0 : y0 + round(z), x0 : x0 + round(z)], um, ZOOM_BAR_UM)
            if i == 0:
                ax_crop.set_title(METHOD_LABELS[method], fontsize=7, loc="left")
        axes[i, 0].text(-0.04, 0.5, f"Crop {i + 1}", transform=axes[i, 0].transAxes, rotation=90, ha="right",
                        va="center", fontsize=7)
    save(fig, f"he_crops_{patient}")


def save(fig: plt.Figure, name: str) -> None:
    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"wrote {RESULTS_DIR / name}.png", flush=True)


def main() -> None:
    apply_style()
    crops = pd.read_csv(RESULTS_DIR / "crops.csv")
    plot_overview(crops)
    for patient in sorted(crops["patient_id"].unique()):
        plot_crops(crops, patient)


if __name__ == "__main__":
    main()
