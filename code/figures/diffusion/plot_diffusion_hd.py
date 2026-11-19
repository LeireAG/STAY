#!/usr/bin/env python3
"""Lateral diffusion in the external Visium HD STAY sections (SIMPlex, data/external/SIMPlex/HD), at 8 um bins.

Same logic as plot_diffusion.py, on the square 8 um grid instead of the Visium spots: tissue bins are Space Ranger's
in_tissue bins (no detachment annotation here), the other bins form contiguous empty areas (bins connected through
direct neighbours), and empty areas smaller than plot_diffusion.py's MIN_AREA_SPOTS Visium spots (as area) are left
out. Sections in TOP_LEFT_HALF are sliced along a diagonal and only their top-left half is used. Each bin gets its
distance to the tissue edge on the grid: empty bins to the nearest tissue bin (positive), tissue bins to
the nearest empty bin (negative). UMIs and detected genes come from the raw 8 um matrix; per section and per
8 um distance step, the mean is taken (most empty bins far from the tissue have 0 UMIs, so their median is 0), absolute
and as % of the section's mean tissue bin.

Writes results/diffusion/hd/diffusion_hd_per_section.csv, results/diffusion/hd/diffusion_hd.{pdf,png} and, per section,
results/diffusion/hd/sections/<section>.png: the bins on the hires H&E coloured by distance, next to the CytAssist and
H&E images.

  python code/figures/diffusion/plot_diffusion_hd.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_diffusion import (MAX_DISTANCE_UM, METRICS, MIN_AREA_SPOTS, MIN_DISTANCE_UM, MIN_OUTSIDE_SPOTS, PROJECT_ROOT,
                            SPOT_PITCH_UM, apply_style, distance_ticks, read_image, style_ax)

HD_DIR = PROJECT_ROOT / "data/external/SIMPlex/HD"
RESULTS = PROJECT_ROOT / "results/diffusion/hd"
SECTIONS = ["patient4_HD", "patient5_HD"]
BIN_UM = 8
# the Visium thresholds as area: one Visium spot covers 100 x 100 x sqrt(3) / 2 um2 of the hexagonal grid
BINS_PER_VISIUM_SPOT = SPOT_PITCH_UM**2 * np.sqrt(3) / 2 / BIN_UM**2
MIN_OUTSIDE_BINS = int(MIN_OUTSIDE_SPOTS * BINS_PER_VISIUM_SPOT)
MIN_AREA_BINS = int(MIN_AREA_SPOTS * BINS_PER_VISIUM_SPOT)
COLORS = {"patient4_HD": "#E41A1C", "patient5_HD": "#377EB8"}
# sections sliced along the diagonal from the bottom-left to the top-right corner of the capture area (as in the
# section PNG): only the top-left half is used, where the tissue edges border real glass
TOP_LEFT_HALF = {"patient5_HD"}


def bin_dir(section: str) -> Path:
    return HD_DIR / section / "binned_outputs/square_008um"


def read_bins(section: str) -> pd.DataFrame:
    """All 8 um bins: grid position, tissue call, position on the hires H&E, UMIs and detected genes (raw matrix)."""
    bins = pd.read_parquet(bin_dir(section) / "spatial/tissue_positions.parquet")
    scale = json.loads((bin_dir(section) / "spatial/scalefactors_json.json").read_text())["tissue_hires_scalef"]
    bins["hires_x"] = bins["pxl_col_in_fullres"] * scale
    bins["hires_y"] = bins["pxl_row_in_fullres"] * scale
    with h5py.File(bin_dir(section) / "raw_feature_bc_matrix.h5", "r") as f:
        barcodes = f["matrix/barcodes"][:].astype(str)
        data, indptr = f["matrix/data"][:], f["matrix/indptr"][:]
    cumulative = np.concatenate([[0], np.cumsum(data, dtype=np.int64)])
    metrics = pd.DataFrame({"barcode": barcodes, "total_counts": cumulative[indptr[1:]] - cumulative[indptr[:-1]],
                            "n_genes": np.diff(indptr)})
    return bins.merge(metrics, on="barcode")


def section_profile(section: str) -> pd.DataFrame | None:
    """Mean UMIs and genes per distance from the tissue edge, absolute and as % of the tissue mean."""
    bins = read_bins(section)
    shape = (bins["array_row"].max() + 1, bins["array_col"].max() + 1)
    rows, cols = bins["array_row"].to_numpy(), bins["array_col"].to_numpy()
    tissue = np.zeros(shape, bool)
    tissue[rows, cols] = bins["in_tissue"].astype(bool)
    keep = np.ones(shape, bool)
    if section in TOP_LEFT_HALF:
        u = (bins["hires_x"] - bins["hires_x"].min()) / (bins["hires_x"].max() - bins["hires_x"].min())
        v = (bins["hires_y"] - bins["hires_y"].min()) / (bins["hires_y"].max() - bins["hires_y"].min())
        keep[rows, cols] = (u + v < 1).to_numpy()
    outside = ~tissue & keep
    if outside.sum() < MIN_OUTSIDE_BINS:
        print(f"{section}: {outside.sum()} bins outside tissue, left out", flush=True)
        return None

    # contiguous empty areas (direct neighbours); small areas are left out
    areas, _ = ndimage.label(outside)
    sizes = np.bincount(areas.ravel())
    empty = outside & (sizes[areas] >= MIN_AREA_BINS)
    if not empty.any():
        print(f"{section}: no empty area of >= {MIN_AREA_BINS} bins, left out", flush=True)
        return None

    # signed distance to the tissue edge, in um: empty bins + (to the nearest tissue bin), tissue bins - (to the
    # nearest empty bin), over the whole section; only bins kept by the slice get one
    distance = np.full(shape, np.nan)
    distance[empty] = ndimage.distance_transform_edt(~tissue)[empty]
    distance[tissue & keep] = -ndimage.distance_transform_edt(~empty)[tissue & keep]
    bins["distance_um"] = np.round(distance[rows, cols]) * BIN_UM
    shown = bins[bins["distance_um"].between(MIN_DISTANCE_UM, MAX_DISTANCE_UM)]
    plot_section_map(section, bins, shown)

    out = shown.groupby("distance_um")[list(METRICS)].mean()
    out.insert(0, "n_bins", shown.groupby("distance_um").size())
    for metric in METRICS:
        out[f"{metric}_pct_tissue"] = 100 * out[metric] / bins.loc[(tissue & keep)[rows, cols], metric].mean()
    return out.reset_index().assign(section=section)


def plot_section_map(section: str, bins: pd.DataFrame, shown: pd.DataFrame) -> None:
    """Bins on the hires H&E coloured by distance to the tissue edge (negative: inside), next to CytAssist and H&E."""
    image = read_image(bin_dir(section) / "spatial/tissue_hires_image.png")
    fig, (ax, ax_cyt, ax_he) = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    ax.imshow(image, alpha=0.6)
    dots = ax.scatter(shown["hires_x"], shown["hires_y"], s=1, c=shown["distance_um"], cmap="viridis",
                      vmin=MIN_DISTANCE_UM, vmax=MAX_DISTANCE_UM, linewidths=0, rasterized=True)
    fig.colorbar(dots, ax=ax, shrink=0.6, label="Distance from tissue edge (µm)")
    ax.set_title(section)
    ax_cyt.imshow(read_image(bin_dir(section) / "spatial/cytassist_image.tiff"))
    ax_cyt.set_title("CytAssist")
    ax_he.imshow(image)
    ax_he.set_title("H&E")
    for a in (ax, ax_cyt, ax_he):
        a.set_axis_off()
    (RESULTS / "sections").mkdir(parents=True, exist_ok=True)
    fig.savefig(RESULTS / "sections" / f"{section}.png", dpi=200)
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    profiles = [p for s in SECTIONS if (p := section_profile(s)) is not None]
    table = pd.concat(profiles, ignore_index=True)
    table.to_csv(RESULTS / "diffusion_hd_per_section.csv", index=False)

    apply_style()
    fig, axs = plt.subplots(2, 2, figsize=(5.0, 4.2), layout="constrained")
    for col, (metric, label) in enumerate(METRICS.items()):
        label = label.replace("spot", "8 µm bin")
        rows = [(metric, label), (f"{metric}_pct_tissue", f"{label}\n(% of tissue mean)")]
        for row, (column, ylabel) in enumerate(rows):
            ax = axs[row, col]
            for section, profile in table.groupby("section"):
                ax.plot(profile["distance_um"], profile[column], color=COLORS[section], linewidth=1, label=section)
            ax.axvline(0, color="grey", linewidth=0.6, linestyle="--")  # tissue edge
            ax.set(yscale="log", ylabel=ylabel,
                   xlim=(MIN_DISTANCE_UM - BIN_UM, MAX_DISTANCE_UM + BIN_UM))
            if row == 1:
                ax.set_xlabel("Distance from tissue edge (µm)")
                ax.yaxis.set_major_formatter("{x:g}%")
            distance_ticks(ax)
            style_ax(ax)
    axs[0, 0].legend(frameon=False)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"diffusion_hd.{ext}")
    print(f"wrote {RESULTS / 'diffusion_hd.pdf'}")


if __name__ == "__main__":
    main()
