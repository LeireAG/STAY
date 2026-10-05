#!/usr/bin/env python3
"""Single composite cohort figure for detachment and QC metrics.

Outputs one PDF/PNG with:
  - sample-level detachment rate by cohort
  - unique genes, UMI count, and % mitochondrial reads by cohort
  - detachment rate by epithelial vs non-epithelial compartment
  - sample counts shown on cohort x-axis labels

Also writes two summary CSVs:
  - detachment_per_sample_summary.csv: one row per sample
  - detachment_cohort_summary.csv: one row per (metric, cohort) with
    n/mean/median/std/min/max and a Mann-Whitney U cohort comparison

  python scripts/plot_detachment_cohort.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_detachment_pipeline import MATRIX_NAME, SPACERANGER_DIR

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CFG = yaml.safe_load((PROJECT_ROOT / "resources" / "detachment.yaml").read_text(encoding="utf-8"))

SAMPLE_LIST = PROJECT_ROOT / CFG["paths"]["sample_list"]
SCORE_DIR = PROJECT_ROOT / CFG["paths"]["detachment_score"]
EPITHELIAL_DIR = PROJECT_ROOT / CFG["paths"]["epithelial_annotation"]
RESULTS_DIR = PROJECT_ROOT / CFG["paths"]["results"]
COHORT_COLORS = CFG["colors"]["cohort"]
COHORT_ORDER = ["standard", "STAY"]
COHORT_LABELS = {"standard": "standard", "STAY": "STAY"}
SCORE_THRES = float(CFG["scoring"]["score_thres"])
SCORE_MAD_THRES = float(CFG["scoring"]["score_mad_thres"])
MIN_SPOTS = int(CFG["scoring"]["min_spots"])

BOXPLOT_STYLE = dict(
    patch_artist=True,
    showfliers=False,
    medianprops={"color": "black", "linewidth": 1.0},
    boxprops={"linewidth": 0.7},
    whiskerprops={"linewidth": 0.7},
    capprops={"linewidth": 0.7},
)


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 7,
            "axes.titlesize": 7,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "figure.dpi": 150,
        }
    )


def detachment_threshold(values: pd.Series) -> float:
    if not len(values):
        return np.nan
    median = float(values.median())
    mad = float(np.median(np.abs(values - median))) * 1.4826
    return max(SCORE_THRES, median + SCORE_MAD_THRES * mad)


def assign_detached(scores: pd.Series) -> pd.Series:
    values = pd.to_numeric(scores, errors="coerce")
    threshold = detachment_threshold(values.dropna())
    called = (values >= threshold).fillna(False)
    if int(called.sum()) < MIN_SPOTS:
        return pd.Series(False, index=scores.index)
    return called


def call_rate(scores: pd.Series) -> tuple[float, float, int]:
    values = pd.to_numeric(scores, errors="coerce").dropna()
    n_high = int(assign_detached(scores).sum())
    threshold = detachment_threshold(values)
    return 100.0 * n_high / len(values), threshold, n_high


def mannwhitney(a: pd.Series, b: pd.Series) -> dict:
    a = a.dropna()
    b = b.dropna()
    if len(a) < 2 or len(b) < 2:
        return {
            "n_standard": len(a),
            "n_STAY": len(b),
            "median_standard": float(a.median()) if len(a) else np.nan,
            "median_STAY": float(b.median()) if len(b) else np.nan,
            "U": np.nan,
            "p": np.nan,
        }
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_standard": int(len(a)),
        "n_STAY": int(len(b)),
        "median_standard": float(a.median()),
        "median_STAY": float(b.median()),
        "U": float(stat),
        "p": float(p),
    }


def p_to_stars(p: float) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def format_p(p: float) -> str:
    if not np.isfinite(p):
        return "p = NA"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3g}"


def style_ax(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.5)


def cohort_xtick_labels(summaries: pd.DataFrame) -> list[str]:
    return [
        f"{COHORT_LABELS[c]}\n(n={int((summaries['cohort'] == c).sum())})"
        for c in COHORT_ORDER
    ]


def cohort_legend() -> list[Patch]:
    return [Patch(facecolor=COHORT_COLORS[c], edgecolor="black", label=COHORT_LABELS[c]) for c in COHORT_ORDER]


def read_spot_metrics(sample_id: str, barcodes: pd.Series) -> pd.DataFrame:
    path = SPACERANGER_DIR / sample_id / MATRIX_NAME
    with h5py.File(path, "r") as f:
        mat = f["matrix"]
        all_bc = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["barcodes"][:]]
        gene_names = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["features/name"][:]]
        data = mat["data"][:]
        indices = mat["indices"][:]
        indptr = mat["indptr"][:]
    nnz = np.diff(indptr)
    x = csc_matrix((data, indices, indptr), shape=(len(gene_names), len(all_bc)))
    total_counts = np.asarray(x.sum(axis=0)).ravel()
    mt_mask = np.array([g.startswith("MT-") for g in gene_names])
    mt_counts = np.asarray(x[mt_mask, :].sum(axis=0)).ravel() if mt_mask.any() else np.zeros(len(all_bc))
    pct_mt = 100.0 * mt_counts / np.maximum(total_counts, 1)
    df = pd.DataFrame(
        {
            "barcode": all_bc,
            "n_genes": nnz,
            "total_counts": total_counts,
            "pct_mt": pct_mt,
        }
    )
    barcodes_str = barcodes.astype(str)
    return df.loc[df["barcode"].isin(barcodes_str)].copy()


def load_summaries() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-sample and per-(sample, compartment) detachment/QC summaries."""
    samples = pd.read_csv(SAMPLE_LIST)
    sample_rows = []
    comp_rows = []
    for _, row in samples.iterrows():
        sample_id = str(row["sample_id"])
        cohort = str(row["cohort"])
        scores = pd.read_csv(SCORE_DIR / f"{sample_id}.csv")
        scores["barcode"] = scores["barcode"].astype(str)
        scores["detachment_score"] = pd.to_numeric(scores["detachment_score"], errors="coerce")
        scores["detached"] = assign_detached(scores["detachment_score"])
        rate, threshold, n_high = call_rate(scores["detachment_score"])
        metrics = read_spot_metrics(sample_id, scores["barcode"])
        sample_rows.append(
            {
                "sample_id": sample_id,
                "cohort": cohort,
                "n_spots": int(scores["detachment_score"].notna().sum()),
                "n_detached": n_high,
                "detachment_rate_pct": rate,
                "threshold": threshold,
                "median_n_genes": float(metrics["n_genes"].median()),
                "median_total_counts": float(metrics["total_counts"].median()),
                "median_pct_mt": float(metrics["pct_mt"].median()),
            }
        )

        epi_path = EPITHELIAL_DIR / sample_id / "csv" / f"{sample_id}_epithelial.csv"
        if epi_path.exists():
            epi = pd.read_csv(epi_path)
            epi["barcode"] = epi["barcode"].astype(str)
            epi = epi[epi["barcode"].isin(scores["barcode"])]
        else:
            epi = pd.DataFrame({"barcode": [], "annotation": []})

        merged = scores[["barcode", "detached"]].merge(epi, on="barcode", how="left")
        merged["compartment"] = np.where(merged["annotation"].notna(), "epithelial", "non-epithelial")
        for comp, sub in merged.groupby("compartment", sort=False):
            comp_rows.append(
                {
                    "sample_id": sample_id,
                    "cohort": cohort,
                    "compartment": comp,
                    "n_spots": int(len(sub)),
                    "n_detached": int(sub["detached"].sum()),
                    "detachment_rate_pct": 100.0 * float(sub["detached"].sum()) / len(sub),
                }
            )
    return pd.DataFrame(sample_rows), pd.DataFrame(comp_rows)


def describe_stats(values: pd.Series) -> dict:
    values = pd.to_numeric(values, errors="coerce").dropna()
    return {
        "n": int(len(values)),
        "mean": float(values.mean()) if len(values) else np.nan,
        "median": float(values.median()) if len(values) else np.nan,
        "std": float(values.std()) if len(values) else np.nan,
        "min": float(values.min()) if len(values) else np.nan,
        "max": float(values.max()) if len(values) else np.nan,
    }


def build_cohort_summary(per_sample: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "detachment_rate_pct",
        "median_n_genes",
        "median_total_counts",
        "median_pct_mt",
        "epithelial_detachment_rate_pct",
        "nonepithelial_detachment_rate_pct",
    ]
    rows = []
    for metric in metrics:
        stats = mannwhitney(
            per_sample.loc[per_sample["cohort"] == COHORT_ORDER[0], metric],
            per_sample.loc[per_sample["cohort"] == COHORT_ORDER[1], metric],
        )
        for cohort in COHORT_ORDER:
            rows.append(
                {
                    "metric": metric,
                    "cohort": cohort,
                    **describe_stats(per_sample.loc[per_sample["cohort"] == cohort, metric]),
                    "p_value": stats["p"],
                    "p_stars": p_to_stars(stats["p"]),
                }
            )
    return pd.DataFrame(rows)


def add_bracket(ax: plt.Axes, x0: float, x1: float, y: float, h: float, p: float) -> None:
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color="black", lw=0.7, clip_on=False)
    ax.text((x0 + x1) / 2, y + h, p_to_stars(p), ha="center", va="bottom", fontsize=7, clip_on=False)


def plot_cohort_panel(
    ax: plt.Axes,
    summaries: pd.DataFrame,
    metric: str,
    ylabel: str,
    title: str,
    show_legend: bool = False,
) -> dict:
    rng = np.random.default_rng(0)
    positions = [0.0, 1.0]
    data = [summaries.loc[summaries["cohort"] == c, metric].dropna().to_numpy() for c in COHORT_ORDER]
    bp = ax.boxplot(data, positions=positions, widths=0.42, **BOXPLOT_STYLE)
    for patch, cohort in zip(bp["boxes"], COHORT_ORDER):
        patch.set_facecolor(COHORT_COLORS[cohort])
        patch.set_edgecolor("black")
        patch.set_alpha(0.92)
    for pos, pts in zip(positions, data):
        if len(pts):
            ax.scatter(pos + rng.uniform(-0.06, 0.06, len(pts)), pts, s=12, c="black", alpha=0.55, linewidths=0)

    stats = mannwhitney(
        summaries.loc[summaries["cohort"] == COHORT_ORDER[0], metric],
        summaries.loc[summaries["cohort"] == COHORT_ORDER[1], metric],
    )
    y0, y1 = ax.get_ylim()
    y = y1 - 0.08 * (y1 - y0)
    add_bracket(ax, positions[0], positions[1], y, 0.03 * (y1 - y0), stats["p"])
    ax.set_xticks(positions, cohort_xtick_labels(summaries))
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nMann–Whitney U, {format_p(stats['p'])}", pad=8)
    if show_legend:
        ax.legend(handles=cohort_legend(), frameon=False, loc="best", fontsize=6.5)
    style_ax(ax)
    return stats


def plot_compartment_panel(ax: plt.Axes, comp_df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    compartments = sorted([c for c in comp_df["compartment"].unique().tolist() if c != "non-epithelial"])
    if "non-epithelial" in comp_df["compartment"].unique():
        compartments.append("non-epithelial")
    base_x = np.arange(len(compartments), dtype=float) * 2.1
    offsets = {COHORT_ORDER[0]: -0.24, COHORT_ORDER[1]: 0.24}
    stats_rows = []

    for i, comp in enumerate(compartments):
        x_center = base_x[i]
        for cohort in COHORT_ORDER:
            pts = comp_df.loc[
                (comp_df["compartment"] == comp) & (comp_df["cohort"] == cohort),
                "detachment_rate_pct",
            ].dropna().to_numpy()
            if len(pts) == 0:
                continue
            pos = x_center + offsets[cohort]
            bp = ax.boxplot([pts], positions=[pos], widths=0.38, **BOXPLOT_STYLE)
            for patch in bp["boxes"]:
                patch.set_facecolor(COHORT_COLORS[cohort])
                patch.set_edgecolor("black")
                patch.set_alpha(0.92)
            ax.scatter(pos + rng.uniform(-0.05, 0.05, len(pts)), pts, s=10, c="black", alpha=0.55, linewidths=0)

        st = mannwhitney(
            comp_df.loc[(comp_df["compartment"] == comp) & (comp_df["cohort"] == COHORT_ORDER[0]), "detachment_rate_pct"],
            comp_df.loc[(comp_df["compartment"] == comp) & (comp_df["cohort"] == COHORT_ORDER[1]), "detachment_rate_pct"],
        )
        stats_rows.append({"compartment": comp, **st, "stars": p_to_stars(st["p"])})

    y0, y1 = ax.get_ylim()
    for i, st in enumerate(stats_rows):
        if np.isfinite(st["p"]):
            ymax = float(comp_df.loc[comp_df["compartment"] == st["compartment"], "detachment_rate_pct"].max())
            ax.text(base_x[i], ymax + 0.04 * (y1 - y0 + 1e-9), st["stars"], ha="center", va="bottom", fontsize=7)

    ax.set_xticks(base_x, compartments)
    ax.set_ylabel("Detached spots (%)")
    ax.set_title("Detachment by compartment", pad=8)
    ax.legend(handles=cohort_legend(), frameon=False, loc="upper left", fontsize=6.5)
    style_ax(ax)
    return pd.DataFrame(stats_rows)


def main() -> None:
    apply_style()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    summaries, comp_df = load_summaries()

    epi = comp_df.loc[comp_df["compartment"] == "epithelial"].set_index("sample_id")
    epi = epi[["n_spots", "n_detached", "detachment_rate_pct"]].add_prefix("epithelial_")
    non = comp_df.loc[comp_df["compartment"] == "non-epithelial"].set_index("sample_id")
    non = non[["n_spots", "n_detached", "detachment_rate_pct"]].add_prefix("nonepithelial_")
    per_sample = summaries.set_index("sample_id").join(epi).join(non).reset_index()
    per_sample.to_csv(RESULTS_DIR / "detachment_per_sample_summary.csv", index=False)
    build_cohort_summary(per_sample).to_csv(RESULTS_DIR / "detachment_cohort_summary.csv", index=False)

    fig = plt.figure(figsize=(8.6, 6.0))
    gs = GridSpec(
        2,
        8,
        figure=fig,
        height_ratios=[1.45, 1.0],
        hspace=0.58,
        wspace=0.40,
        top=0.96,
        bottom=0.11,
        left=0.07,
        right=0.99,
    )

    ax_rate = fig.add_subplot(gs[0, 0:2])
    ax_genes = fig.add_subplot(gs[0, 2:4])
    ax_umi = fig.add_subplot(gs[0, 4:6])
    ax_mt = fig.add_subplot(gs[0, 6:8])
    ax_comp = fig.add_subplot(gs[1, 0:3])
    comp_pos = ax_comp.get_position()
    ax_comp.set_position([comp_pos.x0, comp_pos.y0, comp_pos.width * 0.7, comp_pos.height])

    plot_cohort_panel(
        ax_rate,
        summaries,
        "detachment_rate_pct",
        "Detached spots (%)",
        "Detachment rate",
        show_legend=True,
    )
    plot_cohort_panel(
        ax_genes,
        summaries,
        "median_n_genes",
        "Unique genes",
        "Unique genes",
    )
    plot_cohort_panel(
        ax_umi,
        summaries,
        "median_total_counts",
        "UMI count",
        "UMI count",
    )
    plot_cohort_panel(
        ax_mt,
        summaries,
        "median_pct_mt",
        "Mitochondrial reads (%)",
        "% mitochondrial",
    )
    plot_compartment_panel(ax_comp, comp_df)

    out_pdf = RESULTS_DIR / "detachment_cohort.pdf"
    out_png = RESULTS_DIR / "detachment_cohort.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"wrote {out_pdf}", flush=True)
    print(f"wrote {out_png}", flush=True)


if __name__ == "__main__":
    main()
