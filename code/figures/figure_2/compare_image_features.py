#!/usr/bin/env python3
"""Compare H&E image features between the STAY and standard cohorts.

Reads the per-spot H&E features written by run_detachment_pipeline.py (he_* columns) and the epithelial
annotation. Sections come from different tissues and there are no replicates, so spots are never the unit of
comparison; patients are (sections of the same patient are pooled). Per compartment (all, epithelial,
non-epithelial):
  - per patient: median (feature level) and IQR (intra-tissue spread) of each feature
  - Mann-Whitney U between cohorts on those per-patient values, with Cliff's delta and BH-FDR over all tests
  - spot-level variance split into between-cohort, between-patient (within cohort) and within-patient parts,
    to show how large any cohort difference is next to the variation inside each tissue

Writes to results/figure_2:
  - image_features_per_patient.csv, image_features_stats.csv, image_features_variance.csv
  - image_features.pdf/png

  python code/figures/figure_2/compare_image_features.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_detachment_cohort import BOXPLOT_STYLE, COHORT_COLORS, COHORT_ORDER, RESULTS_DIR, apply_style, cohort_legend, style_ax
from run_detachment_pipeline import SAMPLES, cohort_dir

FEATURES = {
    "he_r": "Red",
    "he_g": "Green",
    "he_b": "Blue",
    "he_brightness": "Brightness",
    "he_ASM": "Texture (ASM)",
}
COMPARTMENTS = ["all", "epithelial", "non-epithelial"]
# patients with fewer spots in a compartment are left out of that compartment
MIN_SPOTS = 20


def load_spots() -> pd.DataFrame:
    """One row per spot and compartment: every spot is in "all", annotated sections also in their compartment."""
    frames = []
    for _, row in SAMPLES.iterrows():
        sample_id = str(row["sample_id"])
        spots = pd.read_csv(cohort_dir("detachment_score", sample_id) / f"{sample_id}.csv", dtype={"barcode": str})
        epi_path = cohort_dir("histopathology_annotation", sample_id) / sample_id / "csv" / f"{sample_id}_epithelial.csv"
        if epi_path.exists():
            epi = pd.read_csv(epi_path, dtype={"barcode": str}).dropna(subset=["annotation"])
            spots["compartment"] = np.where(spots["barcode"].isin(epi["barcode"]), "epithelial", "non-epithelial")
        else:
            spots["compartment"] = np.nan
        frames.append(spots.assign(patient_id=str(row["patient_id"]), cohort=str(row["cohort"])))
    spots = pd.concat(frames, ignore_index=True)
    return pd.concat([spots.assign(compartment="all"), spots.dropna(subset=["compartment"])], ignore_index=True)


def per_patient(spots: pd.DataFrame) -> pd.DataFrame:
    grouped = spots.groupby(["compartment", "cohort", "patient_id"])[list(FEATURES)]
    out = pd.concat(
        [
            grouped.size().rename("n_spots"),
            grouped.median().add_suffix("_median"),
            (grouped.quantile(0.75) - grouped.quantile(0.25)).add_suffix("_iqr"),
        ],
        axis=1,
    ).reset_index()
    return out.loc[out["n_spots"] >= MIN_SPOTS]


def compare_cohorts(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for comp, sub in patients.groupby("compartment", sort=False):
        for feature in FEATURES:
            for stat in ["median", "iqr"]:
                a = sub.loc[sub["cohort"] == "standard", f"{feature}_{stat}"].dropna()
                b = sub.loc[sub["cohort"] == "STAY", f"{feature}_{stat}"].dropna()
                u, p = mannwhitneyu(a, b, alternative="two-sided")
                rows.append(
                    {
                        "compartment": comp,
                        "feature": feature,
                        "per_patient": stat,
                        "n_standard": len(a),
                        "n_STAY": len(b),
                        "median_standard": float(a.median()),
                        "median_STAY": float(b.median()),
                        "U": float(u),
                        "p": float(p),
                        # P(STAY > standard) - P(STAY < standard), from -1 to 1
                        "cliffs_delta": 1 - 2 * float(u) / (len(a) * len(b)),
                    }
                )
    stats = pd.DataFrame(rows)
    stats["q"] = false_discovery_control(stats["p"])
    return stats


def variance_split(spots: pd.DataFrame) -> pd.DataFrame:
    """Fractions of the spot-level sum of squares; patients are nested in cohorts so the three parts add up to 1."""
    rows = []
    for comp, sub in spots.groupby("compartment", sort=False):
        for feature in FEATURES:
            x = sub[[feature, "cohort", "patient_id"]].dropna()
            values = x[feature]
            grand = values.mean()
            cohort_mean = x.groupby("cohort")[feature].transform("mean")
            patient_mean = x.groupby("patient_id")[feature].transform("mean")
            total = ((values - grand) ** 2).sum()
            rows.append(
                {
                    "compartment": comp,
                    "feature": feature,
                    "between_cohort": ((cohort_mean - grand) ** 2).sum() / total,
                    "between_patient": ((patient_mean - cohort_mean) ** 2).sum() / total,
                    "within_patient": ((values - patient_mean) ** 2).sum() / total,
                }
            )
    return pd.DataFrame(rows)


def plot_feature(ax: plt.Axes, sub: pd.DataFrame, feature: str, q: float) -> None:
    rng = np.random.default_rng(0)
    data = [sub.loc[sub["cohort"] == c, f"{feature}_median"].dropna().to_numpy() for c in COHORT_ORDER]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.42, **BOXPLOT_STYLE)
    for patch, cohort in zip(bp["boxes"], COHORT_ORDER):
        patch.set_facecolor(COHORT_COLORS[cohort])
        patch.set_edgecolor("black")
    for pos, pts in enumerate(data):
        ax.scatter(pos + rng.uniform(-0.06, 0.06, len(pts)), pts, s=8, c="black", alpha=0.55, linewidths=0)
    ax.set_xticks([0, 1], [f"{c}\n(n={len(d)})" for c, d in zip(COHORT_ORDER, data)])
    ax.set_title(f"{FEATURES[feature]}\nq = {q:.2g}")
    style_ax(ax)


def plot_variance(ax: plt.Axes, var: pd.DataFrame) -> None:
    parts = {"between_cohort": "cohort", "between_patient": "patient", "within_patient": "within tissue"}
    colors = ["#E41A1C", "#999999", "#DDDDDD"]
    y = np.arange(len(var))
    left = np.zeros(len(var))
    for (col, label), color in zip(parts.items(), colors):
        ax.barh(y, var[col], left=left, color=color, edgecolor="black", linewidth=0.4, label=label)
        left += var[col].to_numpy()
    ax.set_yticks(y, [FEATURES[f] for f in var["feature"]])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of spot variance")
    ax.set_title("Variance split")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> None:
    apply_style()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    spots = load_spots()
    patients = per_patient(spots)
    stats = compare_cohorts(patients)
    var = variance_split(spots)
    patients.to_csv(RESULTS_DIR / "image_features_per_patient.csv", index=False)
    stats.to_csv(RESULTS_DIR / "image_features_stats.csv", index=False)
    var.to_csv(RESULTS_DIR / "image_features_variance.csv", index=False)

    comps = [c for c in COMPARTMENTS if c in set(patients["compartment"])]
    fig, axes = plt.subplots(len(comps), len(FEATURES) + 1, figsize=(11, 2.3 * len(comps)), squeeze=False)
    q = stats.loc[stats["per_patient"] == "median"].set_index(["compartment", "feature"])["q"]
    for i, comp in enumerate(comps):
        sub = patients.loc[patients["compartment"] == comp]
        for j, feature in enumerate(FEATURES):
            plot_feature(axes[i, j], sub, feature, q[(comp, feature)])
        axes[i, 0].set_ylabel(f"{comp}\nper-patient median")
        plot_variance(axes[i, -1], var.loc[var["compartment"] == comp])
    axes[0, 0].legend(handles=cohort_legend(), frameon=False, loc="best", fontsize=6)
    axes[-1, -1].legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=3)
    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = RESULTS_DIR / f"image_features.{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=300)
        print(f"wrote {out}", flush=True)
    plt.close(fig)

    print(stats.sort_values("q").to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
