#!/usr/bin/env python3
"""Compare and plot the H&E image metrics written by he_quality.py, without touching the images.

Reads results/he_quality/he_quality_sections.csv. The unit is the patient (four pairs of one STAY and one
standard section). Per metric: a two-sided paired t-test of the raw differences STAY - standard (3 degrees of
freedom), the mean difference with its 95% confidence interval, and Benjamini-Hochberg q-values over all metrics.

Writes to results/he_quality:
  - he_quality_paired.csv: per metric and patient, both values and the difference
  - he_quality_stats.csv: per metric, mean difference, 95% CI, t, P and q
  - he_quality.pdf/png: one panel per metric with the paired raw values, double Nature column (183 mm)

  python code/figures/he_quality/plot_he_quality.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from scipy.stats import false_discovery_control, t, ttest_rel

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"

MM = 1 / 25.4
# group: metric: y-axis label
METRICS = {
    "Intensity": {
        "brightness": "Brightness (grey level)",
        "rms_contrast": "RMS contrast (grey level)",
        "dynamic_range": "Dynamic range (grey level)",
        "entropy": "Entropy (bits)",
    },
    "Sharpness": {
        "laplacian_var": "Laplacian variance (grey levels²)",
        "tenengrad": "Tenengrad (grey levels²)",
    },
    "Texture": {
        "glcm_contrast": "GLCM contrast",
        "glcm_homogeneity": "GLCM homogeneity",
        "glcm_energy": "GLCM energy",
        "glcm_correlation": "GLCM correlation",
    },
    "Colour": {
        "mean_red": "Red channel (0–255)",
        "mean_green": "Green channel (0–255)",
        "mean_blue": "Blue channel (0–255)",
        "saturation": "Saturation (0–255)",
    },
}
METHODS = ["standard", "STAY"]
# Okabe-Ito
PATIENT_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]


def paired_values(sections: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metrics in METRICS.values():
        for metric in metrics:
            wide = sections.pivot(index="patient_id", columns="method", values=metric)
            for patient, row in wide.iterrows():
                rows.append({"metric": metric, "patient_id": patient, "standard": row["standard"], "STAY": row["STAY"],
                             "STAY_minus_standard": row["STAY"] - row["standard"]})
    return pd.DataFrame(rows)


def paired_stats(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric, sub in pairs.groupby("metric", sort=False):
        d = sub["STAY_minus_standard"].to_numpy()
        res = ttest_rel(sub["STAY"], sub["standard"])
        half = t.ppf(0.975, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
        rows.append({"metric": metric, "n_patients": len(d), "mean_standard": sub["standard"].mean(),
                     "mean_STAY": sub["STAY"].mean(), "mean_difference": d.mean(), "ci95_low": d.mean() - half,
                     "ci95_high": d.mean() + half, "t": res.statistic, "p": res.pvalue})
    stats = pd.DataFrame(rows)
    stats["q"] = false_discovery_control(stats["p"])
    return stats


def format_q(q: float) -> str:
    return f"{q:.2f}" if q >= 0.01 else f"{q:.0e}".replace("e-0", "e-")


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "font.size": 6,
            "axes.labelsize": 6,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.5,
            "legend.fontsize": 5.5,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "ytick.major.size": 2,
            "xtick.major.size": 2,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def plot_panels(pairs: pd.DataFrame, stats: pd.DataFrame) -> None:
    patients = sorted(pairs["patient_id"].unique())
    colors = dict(zip(patients, PATIENT_COLORS))
    # patient 1 is the only pair scanned on the same scanner on the same day
    markers = {p: ("D" if i == 0 else "o") for i, p in enumerate(patients)}
    stats = stats.set_index("metric")
    ncol = max(len(m) for m in METRICS.values())

    fig = plt.figure(figsize=(183 * MM, 150 * MM))
    grid = GridSpec(len(METRICS), ncol, figure=fig, wspace=0.75, hspace=0.45, left=0.09, right=0.99, top=0.97,
                    bottom=0.04)
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    for i, (group, metrics) in enumerate(METRICS.items()):
        for j, (metric, label) in enumerate(metrics.items()):
            ax = fig.add_subplot(grid[i, j])
            ax.annotate(next(letters), xy=(0, 1), xycoords="axes fraction", xytext=(-24, 6), textcoords="offset points",
                        ha="left", va="bottom", fontweight="bold", fontsize=8)
            for _, r in pairs.loc[pairs["metric"] == metric].iterrows():
                p = r["patient_id"]
                ax.plot([0, 1], [r["standard"], r["STAY"]], color=colors[p], lw=0.8, marker=markers[p],
                        ms=3 if markers[p] == "D" else 3.5, mec="none", clip_on=False)
            ax.set_xticks([0, 1], ["Standard", "STAY"])
            ax.set_xlim(-0.25, 1.25)
            ax.set_ylabel(label)
            ax.text(0.5, 1.0, f"$q$ = {format_q(stats.loc[metric, 'q'])}", transform=ax.transAxes, ha="center",
                    va="bottom", fontsize=5.5)
            ax.spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax.annotate(group, xy=(0, 0.5), xycoords="axes fraction", xytext=(-38, 0), textcoords="offset points",
                            ha="center", va="center", rotation=90, fontweight="bold", fontsize=6.5)
        if len(metrics) < ncol:  # use the free space of a short row for the legend
            lax = fig.add_subplot(grid[i, len(metrics) :])
            lax.axis("off")
            handles = [Line2D([], [], color=colors[p], lw=0.8, marker=markers[p], ms=3 if markers[p] == "D" else 3.5,
                              mec="none", label=p.replace("_", " ")) for p in patients]
            lax.legend(handles=handles, loc="center", frameon=False, handlelength=2)

    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_quality.{ext}", dpi=600)
    plt.close(fig)


def main() -> None:
    apply_style()
    sections = pd.read_csv(RESULTS_DIR / "he_quality_sections.csv")
    pairs = paired_values(sections)
    stats = paired_stats(pairs)
    pairs.to_csv(RESULTS_DIR / "he_quality_paired.csv", index=False)
    stats.to_csv(RESULTS_DIR / "he_quality_stats.csv", index=False)
    plot_panels(pairs, stats)
    print(stats.round(4).to_string(index=False), flush=True)
    print(f"wrote {RESULTS_DIR}", flush=True)


if __name__ == "__main__":
    main()
