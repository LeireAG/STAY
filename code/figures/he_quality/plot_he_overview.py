#!/usr/bin/env python3
"""Arrange the thumbnails and crops written by he_crops.py into one figure, without touching the raw images.

Writes results/he_quality/he_images.pdf/png: one row per section (patients in order, standard above STAY) with the
whole section (boxes mark the crops), then each crop at the metrics' pixel size followed by a zoom on its centre
(box in the crop). Scale bars, OVERVIEW_BAR_UM on sections, CROP_BAR_UM on crops and ZOOM_BAR_UM on zooms.

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


def plot_all(crops: pd.DataFrame) -> None:
    sections = crops.drop_duplicates("sample_id").copy()
    sections["order"] = sections["method"].map({m: i for i, m in enumerate(METHODS)})
    sections = sections.sort_values(["patient_id", "order"])
    n_crops = int(crops["crop"].max())
    ncol = 1 + 2 * n_crops
    fig, axes = plt.subplots(len(sections), ncol, figsize=(183 * MM, 33 * MM * len(sections)), squeeze=False,
                             gridspec_kw={"wspace": 0.05, "hspace": 0.06})

    for i, (_, section) in enumerate(sections.iterrows()):
        sample_id = section["sample_id"]
        sub = crops.loc[crops["sample_id"] == sample_id]
        um = section["thumb_um_per_px"]
        ax = axes[i, 0]
        show(ax, plt.imread(RESULTS_DIR / "thumbnails" / f"{sample_id}.png"), um, OVERVIEW_BAR_UM)
        for _, c in sub.iterrows():
            size = c["size_um"] / um
            ax.add_patch(Rectangle((c["x_um"] / um - size / 2, c["y_um"] / um - size / 2), size, size,
                                   fill=False, ec=BOX_COLOR, lw=0.6))
            ax.text(c["x_um"] / um + size / 2, c["y_um"] / um - size / 2, str(c["crop"]), color=BOX_COLOR,
                    ha="left", va="bottom", fontsize=5, fontweight="bold")
        ax.text(-0.06, 0.5, f"{section['patient_id'].replace('_', ' ')}\n{METHOD_LABELS[section['method']]}",
                transform=ax.transAxes, ha="right", va="center", fontsize=6)

        for _, c in sub.iterrows():
            img = plt.imread(RESULTS_DIR / "crops" / f"{sample_id}_{c['crop']}.png")
            cum = c["crop_um_per_px"]
            h, w = img.shape[:2]
            z = ZOOM_UM / cum
            x0, y0 = round((w - z) / 2), round((h - z) / 2)
            ax_crop, ax_zoom = axes[i, 2 * c["crop"] - 1], axes[i, 2 * c["crop"]]
            show(ax_crop, img, cum, CROP_BAR_UM)
            ax_crop.add_patch(Rectangle((x0, y0), z, z, fill=False, ec=BOX_COLOR, lw=0.6))
            show(ax_zoom, img[y0 : y0 + round(z), x0 : x0 + round(z)], cum, ZOOM_BAR_UM)

    headers = ["Section"] + [label for n in range(1, n_crops + 1) for label in (f"Crop {n}", "Zoom")]
    for ax, header in zip(axes[0], headers):
        ax.set_title(header, fontsize=7)

    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_images.{ext}", bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"wrote {RESULTS_DIR / 'he_images.pdf'}", flush=True)


def main() -> None:
    apply_style()
    plot_all(pd.read_csv(RESULTS_DIR / "crops.csv"))


if __name__ == "__main__":
    main()
