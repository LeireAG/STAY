#!/usr/bin/env python3
"""Lateral diffusion: UMIs and genes in spots outside the tissue, by distance from the tissue edge, STAY vs standard.

For every STAY and standard section with spots outside the tissue: tissue spots are Space Ranger's in_tissue spots plus
the annotated detached spots (run_detachment_pipeline.py), all other spots of the capture area are outside. Each spot
gets its distance to the tissue edge from the Visium array grid (100 um between neighbouring spots): outside spots to
the nearest tissue spot (positive), tissue spots to the nearest outside spot (negative, down to MIN_DISTANCE_UM, for
scale), and its UMIs and detected genes from the raw matrix. Per section, the median per DISTANCE_BIN_UM ring (with
at least MIN_SPOTS spots) is taken, absolute and as % of the section's median tissue spot; per cohort, the median and
interquartile range over sections are plotted.

Writes results/diffusion/diffusion_per_section.csv, results/diffusion/diffusion.{pdf,png} and, per section,
results/diffusion/sections/<sample_id>.png: the spots on the hires H&E coloured by their distance to the tissue edge.

  python code/figures/diffusion/plot_diffusion.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figure_2"))
from plot_detachment_cohort import COHORT_COLORS, apply_style, read_spot_metrics, style_ax
from run_detachment_pipeline import COHORT, PROJECT_ROOT, SAMPLES, cohort_dir, read_image, read_spots

RESULTS = PROJECT_ROOT / "results/diffusion"
DISTANCE_BIN_UM = 100
MIN_DISTANCE_UM = -500  # negative: inside the tissue, for scale
MAX_DISTANCE_UM = 1000
MIN_SPOTS = 10  # spots needed in a ring for a section's median
METRICS = {"total_counts": "UMIs per spot", "n_genes": "Genes per spot"}


def section_profile(sample_id: str) -> pd.DataFrame | None:
    """Median UMIs and genes of the outside spots per distance ring, absolute and as % of the tissue median."""
    spots = read_spots(cohort_dir("spaceranger", sample_id) / sample_id)
    calls = pd.read_csv(cohort_dir("detachment", sample_id) / f"{sample_id}.csv")
    spots = spots.merge(calls[["barcode", "tissue"]], on="barcode")
    spots = spots.merge(read_spot_metrics(sample_id, spots["barcode"]), on="barcode")
    if spots["tissue"].all():
        return None

    # array grid in um: neighbours in a row are 2 columns (100 um) apart, rows are 100 * sqrt(3) / 2 um apart
    xy = np.c_[spots["array_col"] * DISTANCE_BIN_UM / 2, spots["array_row"] * DISTANCE_BIN_UM * np.sqrt(3) / 2]
    # signed distance to the tissue edge: outside spots + (to the nearest tissue spot), tissue spots - (to the nearest
    # outside spot)
    tissue = spots["tissue"].to_numpy()
    spots["distance_um"] = 0.0
    spots.loc[~tissue, "distance_um"] = cKDTree(xy[tissue]).query(xy[~tissue])[0]
    spots.loc[tissue, "distance_um"] = -cKDTree(xy[~tissue]).query(xy[tissue])[0]
    d = spots["distance_um"].round(1)
    spots["ring_um"] = np.sign(d) * np.ceil(d.abs() / DISTANCE_BIN_UM) * DISTANCE_BIN_UM
    spots_in_range = spots[spots["ring_um"].between(MIN_DISTANCE_UM, MAX_DISTANCE_UM)]
    plot_section_map(sample_id, spots, spots_in_range)

    rings = spots_in_range.groupby("ring_um")[list(METRICS)].agg(["median", "size"])
    rings = rings[rings[(list(METRICS)[0], "size")] >= MIN_SPOTS]
    out = pd.DataFrame({"n_spots": rings[(list(METRICS)[0], "size")]})
    for metric in METRICS:
        out[metric] = rings[(metric, "median")]
        out[f"{metric}_pct_tissue"] = 100 * out[metric] / spots.loc[tissue, metric].median()
    return out.reset_index().assign(sample_id=sample_id, cohort=COHORT[sample_id])


def plot_section_map(sample_id: str, spots: pd.DataFrame, spots_in_range: pd.DataFrame) -> None:
    """Spots on the hires H&E: coloured by signed distance ring (negative inside the tissue); spots further away grey."""
    image = read_image(cohort_dir("spaceranger", sample_id) / sample_id / "spatial/tissue_hires_image.png")
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(image, alpha=0.6)
    beyond = spots.drop(spots_in_range.index)
    ax.scatter(beyond["hires_x"], beyond["hires_y"], s=1, color="lightgrey", linewidths=0)
    dots = ax.scatter(spots_in_range["hires_x"], spots_in_range["hires_y"], s=1.5, c=spots_in_range["ring_um"],
                      cmap="viridis", vmin=MIN_DISTANCE_UM, vmax=MAX_DISTANCE_UM, linewidths=0)
    fig.colorbar(dots, ax=ax, shrink=0.6, label="Distance from tissue edge (µm)")
    ax.set_title(f"{sample_id} ({COHORT[sample_id]})")
    ax.set_axis_off()
    (RESULTS / "sections").mkdir(parents=True, exist_ok=True)
    fig.savefig(RESULTS / "sections" / f"{sample_id}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    profiles = [p for s in SAMPLES["sample_id"] if (p := section_profile(s)) is not None]
    table = pd.concat(profiles, ignore_index=True)
    table.to_csv(RESULTS / "diffusion_per_section.csv", index=False)
    n = table.groupby("cohort")["sample_id"].nunique()
    print(f"sections with spots outside tissue: {n.to_dict()}")

    apply_style()
    fig, axs = plt.subplots(2, 2, figsize=(5.0, 4.2), layout="constrained")
    for col, (metric, label) in enumerate(METRICS.items()):
        for row, (column, ylabel) in enumerate([(metric, label), (f"{metric}_pct_tissue", f"{label}\n(% of tissue)")]):
            ax = axs[row, col]
            for cohort in ["standard", "STAY"]:
                stats = table[table["cohort"] == cohort].groupby("ring_um")[column].quantile([0.25, 0.5, 0.75])
                stats = stats.unstack()
                color = COHORT_COLORS[cohort]
                ax.fill_between(stats.index, stats[0.25], stats[0.75], color=color, alpha=0.35, linewidth=0)
                ax.plot(stats.index, stats[0.5], "o-", color=color, markersize=2.5, linewidth=1,
                        markeredgecolor="black", markeredgewidth=0.3, label=f"{cohort} (n = {n.get(cohort, 0)})")
            ax.axvline(0, color="grey", linewidth=0.6, linestyle="--")  # tissue edge
            ax.set(yscale="log", ylabel=ylabel,
                   xlim=(MIN_DISTANCE_UM - DISTANCE_BIN_UM / 2, MAX_DISTANCE_UM + DISTANCE_BIN_UM / 2))
            if row == 1:
                ax.set_xlabel("Distance from tissue edge (µm)")
            style_ax(ax)
    axs[0, 0].legend(frameon=False)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"diffusion.{ext}")
    print(f"wrote {RESULTS / 'diffusion.pdf'}")


if __name__ == "__main__":
    main()
