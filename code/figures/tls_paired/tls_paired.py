#!/usr/bin/env python3
"""TLS score of the paired STAY and standard sections, mapped on the H&E.

TLS model of Andersson et al. (2021), https://doi.org/10.1038/s41467-021-26271-2 (adapted from
https://github.com/almaan/her2st/blob/master/scripts/predict.py): counts of the tissue spots normalised to 10,000 per
spot and log1p-transformed, then each spot divided by its total and each gene by its standard deviation across spots,
and the dot product taken with the model coefficients (resources/gene_lists/andersson_2021_coef_full_visium.csv).

The H&E and spot positions are the cropped ones in data/analysis/paired/<sample_id>/ (01_calculate_crop_coords.py),
as for the semla object. Writes data/analysis/paired/tls/<sample_id>.csv (barcode, TLS_score) and
results/tls_paired/tls_paired.{pdf,png}: Patients 4 and 3, standard and STAY, spots with TLS score > TLS_THRESHOLD.

  python code/figures/tls_paired/tls_paired.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.collections import EllipseCollection

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SPACERANGER = PROJECT_ROOT / "data/spaceranger/paired"
CROPPED = PROJECT_ROOT / "data/analysis/paired"
RESULTS = PROJECT_ROOT / "results/tls_paired"
COEFS = PROJECT_ROOT / "resources/gene_lists/andersson_2021_coef_full_visium.csv"

PANELS = [["Patient_4_standard", "Patient_4_STAY"], ["Patient_3_standard", "Patient_3_STAY"]]
TLS_THRESHOLD = 0.01  # TLS-high cut-off: spots above it are shown
VMAX = 0.05  # top of the colour scale


def tls_score(sample_id: str) -> pd.Series:
    """Andersson et al. (2021) TLS model, per tissue spot."""
    adata = sc.read_10x_h5(SPACERANGER / sample_id / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    x = adata.X.toarray().astype(float)
    x = x / x.sum(axis=1, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    x = np.divide(x, std, out=x.copy(), where=std > 0)  # genes with SD 0 left unchanged
    cnt = pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names.str.upper().str.replace(" ", ""))
    cnt["INTERCEPT"] = 1.0
    cnt = cnt.loc[:, ~cnt.columns.duplicated(keep=False)]

    coefs = pd.read_csv(COEFS, index_col=0, header=None).iloc[:, 0]
    coefs.index = coefs.index.str.upper().str.replace(" ", "")
    genes = coefs.index.intersection(cnt.columns)
    print(f"{sample_id}: {len(genes)} / {len(coefs)} model genes present", flush=True)
    return pd.Series(cnt[genes].to_numpy() @ coefs[genes].to_numpy(), index=adata.obs_names, name="TLS_score")


def plot_section(ax, sample_id: str, score: pd.Series):
    scale = json.loads((SPACERANGER / sample_id / "spatial/scalefactors_json.json").read_text())
    spots = pd.read_csv(CROPPED / sample_id / "tissue_positions_shifted.csv", index_col="barcode")
    spots = spots.join(score[score > TLS_THRESHOLD], how="inner").sort_values("TLS_score")  # highest on top
    xy = spots[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy() * scale["tissue_hires_scalef"]
    d = scale["spot_diameter_fullres"] * scale["tissue_hires_scalef"]

    ax.imshow(plt.imread(CROPPED / sample_id / "tissue_hires_image.png"))
    dots = EllipseCollection(d, d, 0, units="xy", offsets=xy, offset_transform=ax.transData,
                             cmap="viridis", norm=plt.Normalize(0, VMAX), linewidths=0)
    dots.set_array(spots["TLS_score"].to_numpy())
    ax.add_collection(dots)
    ax.set_axis_off()
    patient, protocol = sample_id.rsplit("_", 1)
    ax.set_title(f"{patient.replace('_', ' ')} - {protocol}", fontsize=9)
    return dots


def main() -> None:
    (CROPPED / "tls").mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    fig, axs = plt.subplots(2, 2, figsize=(6, 6.4), layout="constrained")
    for row, samples in zip(axs, PANELS):
        for ax, sample_id in zip(row, samples):
            score = tls_score(sample_id)
            score.rename_axis("barcode").to_csv(CROPPED / "tls" / f"{sample_id}.csv")
            dots = plot_section(ax, sample_id, score)

    cbar = fig.colorbar(dots, ax=axs, location="top", shrink=0.4, aspect=12, ticks=[0, VMAX])
    cbar.ax.set_xticklabels(["Low", "High"])
    cbar.ax.invert_xaxis()  # High on the left, as in the original figure
    cbar.set_label("TLS score", fontsize=9)
    for ext in ("pdf", "png"):
        fig.savefig(RESULTS / f"tls_paired.{ext}", dpi=300)
    print(f"wrote {RESULTS / 'tls_paired.pdf'}")


if __name__ == "__main__":
    main()
