#!/usr/bin/env python3
"""Lateral diffusion: UMIs and genes in spots outside the tissue, by distance from the tissue edge, STAY vs standard.

For every STAY and standard section with at least MIN_OUTSIDE_SPOTS spots outside the tissue: tissue spots are Space
Ranger's in_tissue spots plus the annotated detached spots (run_detachment_pipeline.py), so detached spots are never
outside. The other spots form contiguous empty areas (spots connected through direct neighbours); areas of fewer than
MIN_AREA_SPOTS spots, such as small holes in the tissue, are left out. Each spot gets its distance to the tissue edge on
the Visium array grid (100 um between neighbouring spots): empty spots to the nearest tissue spot (positive), tissue
spots to the nearest empty spot (negative, down to MIN_DISTANCE_UM, for scale). UMIs and detected genes come from the
raw matrix. Per section, the median per grid distance is taken, absolute and as % of the section's median
tissue spot; per cohort, the median and interquartile range over sections are plotted. Sections in
EXCLUDED_SECTIONS (tissue missed by Space Ranger's mask) are left out, and their maps from earlier runs deleted.

Writes results/diffusion/standard/diffusion_per_section.csv, results/diffusion/standard/diffusion.{pdf,png} and, per
section, results/diffusion/standard/sections/<sample_id>.png: the spots on the hires H&E coloured by their distance
to the tissue edge, next to the CytAssist and H&E images.

  python code/figures/diffusion/plot_diffusion.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figure_2"))
from plot_detachment_cohort import COHORT_COLORS, apply_style, read_spot_metrics, style_ax
from run_detachment_pipeline import (COHORT, PROJECT_ROOT, SAMPLES, cohort_dir, read_cytassist_image, read_image,
                                     read_spots)

RESULTS = PROJECT_ROOT / "results/diffusion/standard"
SPOT_PITCH_UM = 100  # Visium centre-to-centre distance of neighbouring spots
MIN_DISTANCE_UM = -500  # negative: inside the tissue, for scale
MAX_DISTANCE_UM = 1000
MIN_OUTSIDE_SPOTS = 200  # sections with fewer spots outside the tissue are left out
MIN_AREA_SPOTS = 200  # contiguous empty areas with fewer spots (e.g. small holes in the tissue) are left out
# sections where Space Ranger's tissue mask misses tissue, judged on results/diffusion/standard/sections/<sample_id>.png
EXCLUDED_SECTIONS: list[str] = ["P03", "P06", "P19", "P20", "P28"]
METRICS = {"total_counts": "UMIs per spot", "n_genes": "Genes per spot"}


def section_profile(sample_id: str) -> pd.DataFrame | None:
    """Median UMIs and genes per distance from the tissue edge, absolute and as % of the tissue median."""
    if sample_id in EXCLUDED_SECTIONS:
        (RESULTS / "sections" / f"{sample_id}.png").unlink(missing_ok=True)  # map from an earlier run
        print(f"{sample_id}: excluded (tissue mask)", flush=True)
        return None
    spots = read_spots(cohort_dir("spaceranger", sample_id) / sample_id)
    calls = pd.read_csv(cohort_dir("detachment", sample_id) / f"{sample_id}.csv")
    spots = spots.merge(calls[["barcode", "tissue"]], on="barcode")
    spots = spots.merge(read_spot_metrics(sample_id, spots["barcode"]), on="barcode")
    tissue = spots["tissue"].to_numpy()
    if (~tissue).sum() < MIN_OUTSIDE_SPOTS:
        print(f"{sample_id}: {(~tissue).sum()} spots outside tissue, left out", flush=True)
        return None

    # array grid in um: neighbours in a row are 2 columns (100 um) apart, rows are 100 * sqrt(3) / 2 um apart
    xy = np.c_[spots["array_col"] * SPOT_PITCH_UM / 2, spots["array_row"] * SPOT_PITCH_UM * np.sqrt(3) / 2]
    # contiguous empty areas: outside spots connected through direct neighbours; small areas are left out
    outside = np.flatnonzero(~tissue)
    pairs = cKDTree(xy[outside]).query_pairs(SPOT_PITCH_UM * 1.01, output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(outside), len(outside)))
    area = connected_components(graph, directed=False)[1]
    empty = np.zeros(len(spots), bool)
    empty[outside[np.bincount(area)[area] >= MIN_AREA_SPOTS]] = True
    if not empty.any():
        print(f"{sample_id}: no empty area of >= {MIN_AREA_SPOTS} spots, left out", flush=True)
        return None

    # signed distance to the tissue edge: empty spots + (to the nearest tissue spot), tissue spots - (to the nearest
    # empty spot); spots of the small areas left out get none
    spots["distance_um"] = np.nan
    spots.loc[empty, "distance_um"] = cKDTree(xy[tissue]).query(xy[empty])[0]
    spots.loc[tissue, "distance_um"] = -cKDTree(xy[empty]).query(xy[tissue])[0]
    spots["distance_um"] = spots["distance_um"].round(1)
    shown = spots[spots["distance_um"].between(MIN_DISTANCE_UM, MAX_DISTANCE_UM)]
    plot_section_map(sample_id, spots, shown)

    out = shown.groupby("distance_um")[list(METRICS)].median()
    out.insert(0, "n_spots", shown.groupby("distance_um").size())
    for metric in METRICS:
        out[f"{metric}_pct_tissue"] = 100 * out[metric] / spots.loc[tissue, metric].median()
    return out.reset_index().assign(sample_id=sample_id, cohort=COHORT[sample_id])


def plot_section_map(sample_id: str, spots: pd.DataFrame, shown: pd.DataFrame) -> None:
    """Spots on the hires H&E coloured by distance to the tissue edge (negative: inside; the others grey), next to
    the CytAssist and H&E images, to check the tissue mask."""
    image = read_image(cohort_dir("spaceranger", sample_id) / sample_id / "spatial/tissue_hires_image.png")
    fig, (ax, ax_cyt, ax_he) = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    ax.imshow(image, alpha=0.6)
    beyond = spots.drop(shown.index)
    ax.scatter(beyond["hires_x"], beyond["hires_y"], s=1, color="lightgrey", linewidths=0)
    dots = ax.scatter(shown["hires_x"], shown["hires_y"], s=1.5, c=shown["distance_um"],
                      cmap="viridis", vmin=MIN_DISTANCE_UM, vmax=MAX_DISTANCE_UM, linewidths=0)
    fig.colorbar(dots, ax=ax, shrink=0.6, label="Distance from tissue edge (µm)")
    ax.set_title(f"{sample_id} ({COHORT[sample_id]})")
    ax_cyt.imshow(read_cytassist_image(sample_id))
    ax_cyt.set_title("CytAssist")
    ax_he.imshow(image)
    ax_he.set_title("H&E")
    for a in (ax, ax_cyt, ax_he):
        a.set_axis_off()
    (RESULTS / "sections").mkdir(parents=True, exist_ok=True)
    fig.savefig(RESULTS / "sections" / f"{sample_id}.png", dpi=200)
    plt.close(fig)


def distance_ticks(ax) -> None:
    """Distance axis labelled every 100 um (vertical labels so they fit), minor ticks every 50 um."""
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.tick_params(axis="x", labelrotation=90)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    profiles = [p for s in SAMPLES["sample_id"] if (p := section_profile(s)) is not None]
    table = pd.concat(profiles, ignore_index=True)
    table.to_csv(RESULTS / "diffusion_per_section.csv", index=False)
    n = table.groupby("cohort")["sample_id"].nunique()
    print(f"sections used: {n.to_dict()}")

    apply_style()
    fig, axs = plt.subplots(2, 2, figsize=(5.0, 4.2), layout="constrained")
    for col, (metric, label) in enumerate(METRICS.items()):
        rows = [(metric, label), (f"{metric}_pct_tissue", f"{label}\n(% of tissue median)")]
        for row, (column, ylabel) in enumerate(rows):
            ax = axs[row, col]
            for cohort in ["standard", "STAY"]:
                stats = table[table["cohort"] == cohort].groupby("distance_um")[column].quantile([0.25, 0.5, 0.75])
                stats = stats.unstack()
                color = COHORT_COLORS[cohort]
                ax.fill_between(stats.index, stats[0.25], stats[0.75], color=color, alpha=0.35, linewidth=0)
                ax.plot(stats.index, stats[0.5], "o-", color=color, markersize=2, linewidth=1,
                        markeredgecolor="black", markeredgewidth=0.3, label=f"{cohort} (n = {n.get(cohort, 0)})")
            ax.axvline(0, color="grey", linewidth=0.6, linestyle="--")  # tissue edge
            ax.set(yscale="log", ylabel=ylabel,
                   xlim=(MIN_DISTANCE_UM - SPOT_PITCH_UM / 2, MAX_DISTANCE_UM + SPOT_PITCH_UM / 2))
            if row == 1:
                ax.set_xlabel("Distance from tissue edge (µm)")
            distance_ticks(ax)
            style_ax(ax)
    axs[0, 0].legend(frameon=False)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"diffusion.{ext}")
    print(f"wrote {RESULTS / 'diffusion.pdf'}")


if __name__ == "__main__":
    main()
