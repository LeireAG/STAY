#!/usr/bin/env python3
"""Compare and plot the H&E image metrics written by he_quality.py, without touching the images.

Reads results/he_quality/he_quality_sections.csv. The unit is the patient (four pairs of one STAY and one
standard section). Per metric and patient the effect is log2(STAY / standard); per metric, a two-sided paired
t-test of those four values against 0 (one-sample t-test on the log2 fold changes, 3 degrees of freedom), its 95%
confidence interval, and Benjamini-Hochberg q-values over all metrics.

Writes to results/he_quality:
  - he_quality_paired.csv: per metric and patient, both values and the log2 fold change
  - he_quality_stats.csv: per metric, mean log2 fold change, 95% CI, t, P and q
  - he_quality.pdf/png: forest plot, single Nature column (89 mm)

  python code/figures/he_quality/plot_he_quality.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import false_discovery_control, t, ttest_1samp

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"

MM = 1 / 25.4
METRICS = {
    "Intensity": {
        "brightness": "Brightness",
        "rms_contrast": "RMS contrast",
        "dynamic_range": "Dynamic range",
        "entropy": "Entropy",
    },
    "Sharpness": {
        "laplacian_var": "Laplacian variance",
        "tenengrad": "Tenengrad",
    },
    "Texture (GLCM)": {
        "glcm_contrast": "Contrast",
        "glcm_homogeneity": "Homogeneity",
        "glcm_energy": "Energy",
        "glcm_correlation": "Correlation",
    },
    "Colour": {
        "mean_red": "Red",
        "mean_green": "Green",
        "mean_blue": "Blue",
        "saturation": "Saturation",
    },
}
# Okabe-Ito
PATIENT_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]


def paired_log2fc(sections: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metrics in METRICS.values():
        for metric in metrics:
            wide = sections.pivot(index="patient_id", columns="method", values=metric)
            for patient, row in wide.iterrows():
                rows.append({"metric": metric, "patient_id": patient, "standard": row["standard"], "STAY": row["STAY"],
                             "log2_fc": np.log2(row["STAY"] / row["standard"])})
    return pd.DataFrame(rows)


def paired_stats(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric, sub in pairs.groupby("metric", sort=False):
        x = sub["log2_fc"].to_numpy()
        res = ttest_1samp(x, 0.0)
        half = t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
        rows.append({"metric": metric, "n_patients": len(x), "mean_log2_fc": x.mean(), "ci95_low": x.mean() - half,
                     "ci95_high": x.mean() + half, "t": res.statistic, "p": res.pvalue})
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
            "ytick.labelsize": 6,
            "legend.fontsize": 5.5,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0,
            "xtick.major.size": 2,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def plot_forest(pairs: pd.DataFrame, stats: pd.DataFrame) -> None:
    patients = sorted(pairs["patient_id"].unique())
    colors = dict(zip(patients, PATIENT_COLORS))
    # patient 1 is the only pair scanned on the same scanner on the same day
    markers = {p: ("D" if i == 0 else "o") for i, p in enumerate(patients)}
    offsets = dict(zip(patients, np.linspace(-0.2, 0.2, len(patients))))
    stats = stats.set_index("metric")

    # top to bottom: a header per group, then its metrics
    ypos, headers, yy = {}, [], 0.0
    for i, (group, metrics) in enumerate(METRICS.items()):
        yy -= 0.5 if i else 0
        headers.append((group, yy))
        yy -= 0.9
        for metric in metrics:
            ypos[metric] = yy
            yy -= 1

    fig, ax = plt.subplots(figsize=(89 * MM, 92 * MM))
    label_pad = plt.rcParams["ytick.major.pad"]
    ax.axvline(0, color="0.65", lw=0.5, zorder=0)
    for metric, y in ypos.items():
        s = stats.loc[metric]
        ax.plot([s["ci95_low"], s["ci95_high"]], [y, y], color="black", lw=0.8, solid_capstyle="butt", zorder=2)
        ax.plot(s["mean_log2_fc"], y, marker="|", ms=7, mew=1.2, color="black", zorder=5)
        for _, r in pairs.loc[pairs["metric"] == metric].iterrows():
            p = r["patient_id"]
            ax.scatter(r["log2_fc"], y + offsets[p], s=8 if markers[p] == "D" else 10, marker=markers[p],
                       color=colors[p], edgecolors="none", zorder=4)
        ax.text(1.03, y, format_q(s["q"]), transform=ax.get_yaxis_transform(), ha="left", va="center")
    for group, y in headers:
        ax.annotate(group, xy=(0, y), xycoords=ax.get_yaxis_transform(), xytext=(-label_pad, 0),
                    textcoords="offset points", ha="right", va="center", fontweight="bold")
    ax.text(1.03, headers[0][1], "q", transform=ax.get_yaxis_transform(), ha="left", va="center", fontstyle="italic")

    ax.set_yticks(list(ypos.values()), [METRICS_FLAT[m] for m in ypos])
    ax.set_ylim(yy + 0.4, 0.5)
    ax.set_xlabel("log$_2$ fold change (STAY / standard)")
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

    handles = [Line2D([], [], ls="", marker=markers[p], ms=3.2 if markers[p] == "o" else 2.8, mec="none",
                      color=colors[p], label=p.replace("_", " ")) for p in patients]
    handles.append(Line2D([], [], ls="-", marker="|", ms=5, mew=1.2, lw=0.8, color="black", label="Mean ± 95% CI"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.1), ncol=5, frameon=False,
              handletextpad=0.3, columnspacing=0.9, handlelength=1.4)

    fig.tight_layout()
    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_quality.{ext}", bbox_inches="tight", dpi=600)
    plt.close(fig)


METRICS_FLAT = {m: label for metrics in METRICS.values() for m, label in metrics.items()}


def main() -> None:
    apply_style()
    sections = pd.read_csv(RESULTS_DIR / "he_quality_sections.csv")
    pairs = paired_log2fc(sections)
    stats = paired_stats(pairs)
    pairs.to_csv(RESULTS_DIR / "he_quality_paired.csv", index=False)
    stats.to_csv(RESULTS_DIR / "he_quality_stats.csv", index=False)
    plot_forest(pairs, stats)
    print(stats.round(4).to_string(index=False), flush=True)
    print(f"wrote {RESULTS_DIR}", flush=True)


if __name__ == "__main__":
    main()
