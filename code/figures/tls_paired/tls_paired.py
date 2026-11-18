#!/usr/bin/env python3
"""TLS score of the paired STAY and standard sections, mapped on the H&E.

Scores every tissue spot with the TLS model of Andersson et al. (2021), https://doi.org/10.1038/s41467-021-26271-2
(adapted from https://github.com/almaan/her2st/blob/master/scripts/predict.py): spots with 0 counts dropped, counts
normalised to 10,000 per spot and log1p-transformed, then each spot divided by its total and each gene by its standard
deviation across spots, and the dot product taken with the model coefficients
(resources/gene_lists/andersson_2021_coef_full_visium.csv, gene names mapped to the Visium probe set).

Writes data/analysis/paired/tls/<sample_id>.csv (barcode, TLS_score) and results/tls_paired/tls_paired.{pdf,png}:
Patients 4 and 3, standard and STAY, with the spots above TLS_THRESHOLD (TLS-high cut-off) on the hires H&E.

  python code/figures/tls_paired/tls_paired.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.collections import EllipseCollection

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SPACERANGER = PROJECT_ROOT / "data/spaceranger/paired"
SCORE_DIR = PROJECT_ROOT / "data/analysis/paired/tls"
RESULTS = PROJECT_ROOT / "results/tls_paired"
COEFS = PROJECT_ROOT / "resources/gene_lists/andersson_2021_coef_full_visium.csv"

PANELS = [["Patient_4_standard", "Patient_4_STAY"], ["Patient_3_standard", "Patient_3_STAY"]]
TLS_THRESHOLD = 0.01  # TLS-high cut-off: spots above it are shown
VMAX = 0.05  # top of the colour scale
MARGIN_PX = 20  # hires pixels around the tissue spots


def read_section(sample_id: str):
    adata = sc.read_visium(SPACERANGER / sample_id, library_id=sample_id)  # tissue spots (filtered matrix)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_counts=1)
    return adata


def tls_score(adata, sample_id: str) -> pd.Series:
    """Andersson et al. (2021) TLS model."""
    norm = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)["X"]
    x = np.asarray(sc.pp.log1p(norm).todense(), dtype=float)
    x = x / x.sum(axis=1, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    x = np.divide(x, std, out=x.copy(), where=std > 0)  # genes with SD 0 left unchanged
    cnt = pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names.str.upper().str.replace(" ", ""))
    cnt["INTERCEPT"] = 1.0
    cnt = cnt.loc[:, ~cnt.columns.duplicated(keep=False)]

    coefs = pd.read_csv(COEFS, index_col=0, header=None).iloc[:, 0]
    coefs.index = coefs.index.str.upper().str.replace(" ", "")
    genes = coefs.index.intersection(cnt.columns)
    print(f"{sample_id}: {len(genes)} / {len(coefs)} model genes present")
    return pd.Series(cnt[genes].to_numpy() @ coefs[genes].to_numpy(), index=adata.obs_names, name="TLS_score")


def plot_section(ax, adata, score: pd.Series, sample_id: str):
    spatial = adata.uns["spatial"][sample_id]
    scale = spatial["scalefactors"]["tissue_hires_scalef"]
    xy = adata.obsm["spatial"] * scale  # hires pixels
    x0, y0 = np.floor(xy.min(axis=0) - MARGIN_PX).clip(0).astype(int)
    x1, y1 = np.ceil(xy.max(axis=0) + MARGIN_PX).astype(int)
    ax.imshow(spatial["images"]["hires"][y0:y1, x0:x1], extent=(x0, x1, y1, y0))

    high = (score > TLS_THRESHOLD).to_numpy()
    order = np.argsort(score.to_numpy()[high])  # highest scores drawn on top
    d = spatial["scalefactors"]["spot_diameter_fullres"] * scale
    spots = EllipseCollection(d, d, 0, units="xy", offsets=xy[high][order], offset_transform=ax.transData,
                              cmap="viridis", norm=plt.Normalize(0, VMAX), linewidths=0)
    spots.set_array(score.to_numpy()[high][order])
    ax.add_collection(spots)
    ax.set(xlim=(x0, x1), ylim=(y1, y0), xticks=[], yticks=[])
    for side in ax.spines.values():
        side.set_visible(False)
    patient, protocol = sample_id.rsplit("_", 1)
    ax.set_title(f"{patient.replace('_', ' ')} - {protocol}", fontsize=9)
    return spots


def main() -> None:
    SCORE_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(2, 2, figsize=(6, 6.4), layout="constrained")
    for row, samples in zip(axs, PANELS):
        for ax, sample_id in zip(row, samples):
            adata = read_section(sample_id)
            score = tls_score(adata, sample_id)
            score.rename_axis("barcode").to_csv(SCORE_DIR / f"{sample_id}.csv")
            spots = plot_section(ax, adata, score, sample_id)

    cbar = fig.colorbar(spots, ax=axs, location="top", shrink=0.4, aspect=12, ticks=[0, VMAX])
    cbar.ax.set_xticklabels(["Low", "High"])
    cbar.ax.invert_xaxis()  # High on the left, as in the original figure
    cbar.set_label("TLS score", fontsize=9)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"tls_paired.{ext}", dpi=300)
    print(f"wrote {RESULTS / 'tls_paired.pdf'}")


if __name__ == "__main__":
    main()
