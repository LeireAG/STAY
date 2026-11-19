#!/usr/bin/env python3
"""Lateral diffusion: UMIs and genes in spots outside the tissue, by distance from the tissue edge, STAY vs standard.

For every STAY and standard section with spots outside the tissue: tissue spots are Space Ranger's in_tissue spots plus
the annotated detached spots (run_detachment_pipeline.py), all other spots of the capture area are outside. Each outside
spot gets the distance to its nearest tissue spot, from the Visium array grid (100 um between neighbouring spots), and
its UMIs and detected genes from the raw matrix. Per section, the median per DISTANCE_BIN_UM ring (with at least
MIN_SPOTS spots) is taken, absolute and as % of the section's median tissue spot; per cohort, the median and
interquartile range over sections are plotted.

Writes results/diffusion/diffusion_per_section.csv and results/diffusion/diffusion.{pdf,png}.

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
from run_detachment_pipeline import COHORT, PROJECT_ROOT, SAMPLES, cohort_dir, read_spots

RESULTS = PROJECT_ROOT / "results/diffusion"
DISTANCE_BIN_UM = 100
MAX_DISTANCE_UM = 1000
MIN_SPOTS = 10  # outside spots needed in a ring for a section's median
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
    tissue = spots["tissue"].to_numpy()
    outside = spots.loc[~tissue].copy()
    outside["distance_um"] = cKDTree(xy[tissue]).query(xy[~tissue])[0]
    outside["ring_um"] = np.ceil(outside["distance_um"].round(1) / DISTANCE_BIN_UM) * DISTANCE_BIN_UM
    outside = outside[outside["ring_um"] <= MAX_DISTANCE_UM]

    rings = outside.groupby("ring_um")[list(METRICS)].agg(["median", "size"])
    rings = rings[rings[(list(METRICS)[0], "size")] >= MIN_SPOTS]
    out = pd.DataFrame({"n_spots": rings[(list(METRICS)[0], "size")]})
    for metric in METRICS:
        out[metric] = rings[(metric, "median")]
        out[f"{metric}_pct_tissue"] = 100 * out[metric] / spots.loc[tissue, metric].median()
    return out.reset_index().assign(sample_id=sample_id, cohort=COHORT[sample_id])


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
            ax.set(yscale="log", ylabel=ylabel, xlim=(0, MAX_DISTANCE_UM + DISTANCE_BIN_UM / 2))
            if row == 1:
                ax.set_xlabel("Distance from tissue (µm)")
            style_ax(ax)
    axs[0, 0].legend(frameon=False)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"diffusion.{ext}")
    print(f"wrote {RESULTS / 'diffusion.pdf'}")


if __name__ == "__main__":
    main()
