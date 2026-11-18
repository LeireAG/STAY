#!/usr/bin/env python3
"""TLS score of the paired STAY and standard sections, mapped on the H&E.

TLS model of Andersson et al. (2021), https://doi.org/10.1038/s41467-021-26271-2 (adapted from
https://github.com/almaan/her2st/blob/master/scripts/predict.py): counts of the tissue spots normalised to 10,000 per
spot and log1p-transformed, then each spot divided by its total and each gene by its standard deviation across spots,
and the dot product taken with the model coefficients (resources/gene_lists/andersson_2021_coef_full_visium.csv).

For every paired section in resources/sample_list.csv, writes data/analysis/paired/tls/<sample_id>.csv (barcode,
TLS_score) and results/tls_paired/<sample_id>.pdf: spots with TLS score > TLS_THRESHOLD on the H&E (sc.pl.spatial).
The H&E and spot positions are the cropped ones in data/analysis/paired/<sample_id>/ (01_calculate_crop_coords.py),
as for the semla object.

  python code/figures/tls_paired/tls_paired.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SPACERANGER = PROJECT_ROOT / "data/spaceranger/paired"
CROPPED = PROJECT_ROOT / "data/analysis/paired"
RESULTS = PROJECT_ROOT / "results/tls_paired"
COEFS = PROJECT_ROOT / "resources/gene_lists/andersson_2021_coef_full_visium.csv"
SAMPLES = pd.read_csv(PROJECT_ROOT / "resources/sample_list.csv").query("cohort == 'paired'")["sample_id"]

TLS_THRESHOLD = 0.01  # TLS-high cut-off: spots above it are shown
VMAX = 0.05  # top of the colour scale, the same for all sections
SPOT_SCALE = 1.6  # spots drawn 1.6x their 55 um diameter (88 um, below the 100 um spacing) so they read at panel size


def read_section(sample_id: str):
    """Tissue spots with the cropped hires H&E and the matching (shifted) spot positions."""
    adata = sc.read_10x_h5(SPACERANGER / sample_id / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_counts=1)
    spots = pd.read_csv(CROPPED / sample_id / "tissue_positions_shifted.csv", index_col="barcode")
    adata.obsm["spatial"] = spots.loc[adata.obs_names, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    adata.uns["spatial"] = {sample_id: {
        "images": {"hires": plt.imread(CROPPED / sample_id / "tissue_hires_image.png")},
        "scalefactors": json.loads((SPACERANGER / sample_id / "spatial/scalefactors_json.json").read_text()),
    }}
    return adata


def tls_score(adata) -> np.ndarray:
    """Andersson et al. (2021) TLS model, per spot."""
    norm = sc.pp.log1p(sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)["X"])
    x = norm.toarray().astype(float)
    x = x / x.sum(axis=1, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    x = np.divide(x, std, out=x.copy(), where=std > 0)  # genes with SD 0 left unchanged
    cnt = pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names.str.upper().str.replace(" ", ""))
    cnt["INTERCEPT"] = 1.0
    cnt = cnt.loc[:, ~cnt.columns.duplicated(keep=False)]

    coefs = pd.read_csv(COEFS, index_col=0, header=None).iloc[:, 0]
    coefs.index = coefs.index.str.upper().str.replace(" ", "")
    genes = coefs.index.intersection(cnt.columns)
    print(f"{len(genes)} / {len(coefs)} model genes present", flush=True)
    return cnt[genes].to_numpy() @ coefs[genes].to_numpy()


def main() -> None:
    (CROPPED / "tls").mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    for sample_id in SAMPLES:
        print(sample_id, end=": ", flush=True)
        adata = read_section(sample_id)
        adata.obs["TLS_score"] = tls_score(adata)
        adata.obs["TLS_score"].rename_axis("barcode").to_csv(CROPPED / "tls" / f"{sample_id}.csv")

        adata.obs["TLS score"] = adata.obs["TLS_score"].where(adata.obs["TLS_score"] > TLS_THRESHOLD)  # others hidden
        with plt.rc_context({"font.size": 7, "axes.titlesize": 8}):
            fig, ax = plt.subplots(figsize=(2.4, 2.2))
            sc.pl.spatial(adata, color="TLS score", vmin=0, vmax=VMAX, na_color=(0, 0, 0, 0), ax=ax, show=False,
                          title=sample_id.replace("_", " "), alpha_img=0.6, size=SPOT_SCALE, frameon=False,
                          colorbar_loc=None)
            fig.colorbar(plt.cm.ScalarMappable(plt.Normalize(0, VMAX), "viridis"), ax=ax, shrink=0.5, aspect=10,
                         ticks=[0, VMAX], label="TLS score")
            fig.savefig(RESULTS / f"{sample_id}.pdf", bbox_inches="tight", dpi=300)
            plt.close(fig)
    print(f"wrote {RESULTS}")


if __name__ == "__main__":
    main()
