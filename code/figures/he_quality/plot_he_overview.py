#!/usr/bin/env python3
"""Arrange the thumbnails and crops written by he_crops.py into one figure, without touching the raw images.

Writes results/he_quality/he_images.pdf/png: one row per patient, standard on the left and STAY on the right, each
with the whole section (boxes mark the crops), then each crop at the metrics' pixel size followed by a zoom on its
centre (box in the crop). Crops and zooms are embedded in the PDF at their own pixels, without resampling.

  python code/figures/he_quality/plot_he_overview.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"

MM = 1 / 25.4
PANEL_MM = 34
METHODS = ["standard", "STAY"]
METHOD_LABELS = {"standard": "Standard", "STAY": "STAY"}
ZOOM_UM = 75
SECTION_BAR_UM = 1000
CROP_BAR_UM = 50
ZOOM_BAR_UM = 10
BOX_COLOR = "#0072B2"
# resolution the section thumbnails are resampled to in the PDF and PNG (crops and zooms keep their own pixels)
DPI = 900


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


def show(ax: plt.Axes, img, um_per_px: float, bar_um: float, interpolation: str) -> None:
    """Image with a labelled scale bar in the bottom-right corner."""
    ax.imshow(img, interpolation=interpolation)
    h, w = img.shape[:2]
    bar = bar_um / um_per_px
    x1, y = w * 0.95, h * 0.93
    ax.add_patch(Rectangle((x1 - bar, y), bar, h * 0.015, color="black", lw=0))
    label = f"{bar_um / 1000:g} mm" if bar_um >= 1000 else f"{bar_um:g} µm"
    ax.text(x1 - bar / 2, y - h * 0.01, label, ha="center", va="bottom", fontsize=5,
            path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])
    ax.set_axis_off()


def plot_section(ax: plt.Axes, sample_id: str, sub: pd.DataFrame) -> None:
    um = sub["thumb_um_per_px"].iloc[0]
    # thumbnails are much larger than their panels: let matplotlib downsample them smoothly
    show(ax, plt.imread(RESULTS_DIR / "thumbnails" / f"{sample_id}.png"), um, SECTION_BAR_UM, "antialiased")
    for _, c in sub.iterrows():
        size = c["size_um"] / um
        ax.add_patch(Rectangle((c["x_um"] / um - size / 2, c["y_um"] / um - size / 2), size, size,
                               fill=False, ec=BOX_COLOR, lw=0.6))
        ax.text(c["x_um"] / um + size / 2, c["y_um"] / um - size / 2, str(c["crop"]), color=BOX_COLOR,
                ha="left", va="bottom", fontsize=5, fontweight="bold")


def plot_crop(ax_crop: plt.Axes, ax_zoom: plt.Axes, c: pd.Series) -> None:
    img = plt.imread(RESULTS_DIR / "crops" / f"{c['sample_id']}_{c['crop']}.png")
    um = c["crop_um_per_px"]
    h, w = img.shape[:2]
    z = round(ZOOM_UM / um)
    x0, y0 = round((w - z) / 2), round((h - z) / 2)
    # "none": the PDF gets the crop's own pixels, so zooming in shows exactly what the metrics measured
    show(ax_crop, img, um, CROP_BAR_UM, "none")
    ax_crop.add_patch(Rectangle((x0, y0), z, z, fill=False, ec=BOX_COLOR, lw=0.6))
    show(ax_zoom, img[y0 : y0 + z, x0 : x0 + z], um, ZOOM_BAR_UM, "none")


def plot_all(crops: pd.DataFrame) -> None:
    patients = sorted(crops["patient_id"].unique())
    n_crops = int(crops["crop"].max())
    per_method = 1 + 2 * n_crops
    # a narrow empty column separates standard (left) from STAY (right)
    widths = [1] * per_method + [0.25] + [1] * per_method
    fig = plt.figure(figsize=(PANEL_MM * MM * sum(widths) * 1.04, PANEL_MM * MM * len(patients) * 1.04 + 12 * MM))
    grid = GridSpec(len(patients), len(widths), figure=fig, width_ratios=widths, wspace=0.04, hspace=0.04,
                    left=0.035, right=0.995, top=1 - 12 * MM / fig.get_figheight(), bottom=0.005)
    headers = ["Section"] + [label for n in range(1, n_crops + 1) for label in (f"Crop {n}", "Zoom")]

    for i, patient in enumerate(patients):
        for m, method in enumerate(METHODS):
            first = m * (per_method + 1)
            sub = crops.loc[(crops["patient_id"] == patient) & (crops["method"] == method)].sort_values("crop")
            axes = [fig.add_subplot(grid[i, first + k]) for k in range(per_method)]
            plot_section(axes[0], sub["sample_id"].iloc[0], sub)
            for _, c in sub.iterrows():
                plot_crop(axes[2 * c["crop"] - 1], axes[2 * c["crop"]], c)
            if i == 0:
                for ax, header in zip(axes, headers):
                    ax.set_title(header, fontsize=6, pad=3)
                axes[per_method // 2].annotate(METHOD_LABELS[method], xy=(0.5, 1), xycoords="axes fraction",
                                               xytext=(0, 13), textcoords="offset points", ha="center",
                                               va="bottom", fontsize=7, fontweight="bold")
            if m == 0:
                axes[0].text(-0.05, 0.5, patient.replace("_", " "), transform=axes[0].transAxes, rotation=90,
                             ha="right", va="center", fontsize=7)

    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_images.{ext}", dpi=DPI)
    plt.close(fig)
    print(f"wrote {RESULTS_DIR / 'he_images.pdf'}", flush=True)


def main() -> None:
    apply_style()
    plot_all(pd.read_csv(RESULTS_DIR / "crops.csv"))


if __name__ == "__main__":
    main()
