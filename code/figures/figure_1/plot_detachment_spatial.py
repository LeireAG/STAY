#!/usr/bin/env python3
"""Per-sample spatial panels (PowerPoint-friendly 2x3).

Row 1: H&E; CytAssist; pathology/epithelial annotation
Row 2: detachment score; detachment assignment; unique genes

  python scripts/plot_detachment_spatial.py
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import anndata as ad
import cv2
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

cv2.utils.logging.setLogLevel(cv2.utils.logging.LOG_LEVEL_ERROR)
warnings.filterwarnings("ignore", message="Use `squidpy.pl.spatial_scatter`")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_detachment_cohort import assign_detached
from run_detachment_pipeline import MATRIX_NAME, SPACERANGER_DIR, read_cytassist_coords, read_image, read_spots

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CFG = yaml.safe_load((PROJECT_ROOT / "resources" / "detachment.yaml").read_text(encoding="utf-8"))

SAMPLE_LIST = PROJECT_ROOT / CFG["paths"]["sample_list"]
SCORE_DIR = PROJECT_ROOT / CFG["paths"]["detachment_score"]
EPITHELIAL_DIR = PROJECT_ROOT / CFG["paths"]["epithelial_annotation"]
OUT_DIR = PROJECT_ROOT / CFG["paths"]["results"] / "detachment_spatial"

ASSIGN_PALETTE = {"not detached": "#D0D0D0", "detached": "#C0392B"}


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 8,
            "axes.titlesize": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "figure.dpi": 150,
        }
    )


def score_limits(sample_ids: list[str]) -> tuple[float, float]:
    vals = []
    for sample_id in sample_ids:
        s = pd.to_numeric(pd.read_csv(SCORE_DIR / f"{sample_id}.csv")["detachment_score"], errors="coerce")
        vals.append(s.dropna())
    all_scores = pd.concat(vals, ignore_index=True)
    return float(all_scores.quantile(0.01)), float(all_scores.quantile(0.99))


def align_cytassist_to_hires(
    cyt_img: np.ndarray, hires_xy: np.ndarray, cyt_xy: np.ndarray, out_shape: tuple[int, ...]
) -> np.ndarray:
    valid = np.isfinite(hires_xy).all(axis=1) & np.isfinite(cyt_xy).all(axis=1)
    h, w = int(out_shape[0]), int(out_shape[1])
    if int(valid.sum()) < 3:
        return cyt_img
    src = np.ascontiguousarray(cyt_xy[valid], dtype=np.float32)
    dst = np.ascontiguousarray(hires_xy[valid], dtype=np.float32)
    M, _ = cv2.estimateAffine2D(src, dst, method=cv2.RANSAC, ransacReprojThreshold=5.0)
    if M is None:
        return cyt_img
    return cv2.warpAffine(cyt_img, M, (w, h), flags=cv2.INTER_LINEAR, borderValue=(0, 0, 0))


def read_unique_genes(sample_id: str, barcodes: pd.Series) -> pd.Series:
    path = SPACERANGER_DIR / sample_id / MATRIX_NAME
    with h5py.File(path, "r") as f:
        mat = f["matrix"]
        all_bc = [x.decode() if isinstance(x, bytes) else str(x) for x in mat["barcodes"][:]]
        nnz_per_cell = np.diff(mat["indptr"][:])
    idx = pd.Index(all_bc).get_indexer(barcodes.astype(str))
    keep = idx >= 0
    vals = nnz_per_cell[idx[keep]]
    return pd.Series(vals, index=barcodes.astype(str).to_numpy()[keep], name="n_genes")


def load_path_annotation(sample_id: str, barcodes: pd.Series) -> pd.Series:
    epi_path = EPITHELIAL_DIR / sample_id / "csv" / f"{sample_id}_epithelial.csv"
    barcodes_str = barcodes.astype(str)
    if not epi_path.exists():
        return pd.Series(["non-epithelial"] * len(barcodes_str), index=barcodes_str, name="path_annotation")
    epi = pd.read_csv(epi_path)
    epi["barcode"] = epi["barcode"].astype(str)
    out = pd.Series(["non-epithelial"] * len(barcodes_str), index=barcodes_str, name="path_annotation")
    epi = epi.loc[epi["barcode"].isin(barcodes_str), ["barcode", "annotation"]]
    out.loc[epi["barcode"].to_numpy()] = epi["annotation"].astype(str).to_numpy()
    return out


def build_path_palette(sample_ids: list[str]) -> dict[str, str]:
    labels: set[str] = set()
    for sid in sample_ids:
        epi_path = EPITHELIAL_DIR / sid / "csv" / f"{sid}_epithelial.csv"
        if epi_path.exists():
            epi = pd.read_csv(epi_path, usecols=["annotation"])
            labels.update(epi["annotation"].astype(str).unique().tolist())
    labels = sorted(labels)
    palette: dict[str, str] = {"non-epithelial": "#D0D0D0"}
    cmap = plt.get_cmap("tab20")
    for i, lab in enumerate(labels):
        rgb = np.array(cmap(i % 20)[:3]) * 255.0
        palette[lab] = "#{:02X}{:02X}{:02X}".format(int(rgb[0]), int(rgb[1]), int(rgb[2]))
    return palette


def load_adata(sample_id: str) -> tuple[ad.AnnData, np.ndarray, np.ndarray]:
    sample_dir = SPACERANGER_DIR / sample_id
    scores = pd.read_csv(SCORE_DIR / f"{sample_id}.csv")
    scores["barcode"] = scores["barcode"].astype(str)

    spots = read_cytassist_coords(sample_id, read_spots(sample_dir))
    spots["barcode"] = spots["barcode"].astype(str)
    spots = spots.merge(scores, on="barcode", how="inner").set_index("barcode")

    spots["detached"] = assign_detached(spots["detachment_score"])
    spots["detachment"] = pd.Categorical(
        np.where(spots["detached"], "detached", "not detached"),
        categories=["not detached", "detached"],
    )
    spots["n_genes"] = read_unique_genes(sample_id, scores["barcode"]).reindex(spots.index).to_numpy(dtype=float)
    spots["path_annotation"] = load_path_annotation(sample_id, scores["barcode"]).reindex(spots.index).to_numpy(dtype=str)
    spots["path_annotation"] = pd.Categorical(spots["path_annotation"], categories=sorted(spots["path_annotation"].unique().tolist()))

    spatial = sample_dir / "spatial"
    cyt_path = spatial / "cytassist_image.tiff"
    if not cyt_path.exists():
        cyt_path = spatial / "cytassist_image.tif"
    hires = read_image(spatial / "tissue_hires_image.png")
    cyt = read_image(cyt_path)

    with (spatial / "scalefactors_json.json").open() as handle:
        scalefactors = json.load(handle)
    scalefactors["tissue_cytassist_scalef"] = 1.0

    hires_xy = np.column_stack(
        [
            spots["pxl_col_in_fullres"].to_numpy(dtype=float) * float(scalefactors["tissue_hires_scalef"]),
            spots["pxl_row_in_fullres"].to_numpy(dtype=float) * float(scalefactors["tissue_hires_scalef"]),
        ]
    )
    cyt_xy = spots[["cyt_x", "cyt_y"]].to_numpy(dtype=float)
    cyt_aligned = align_cytassist_to_hires(cyt, hires_xy, cyt_xy, hires.shape)

    adata = ad.AnnData(obs=spots[["detachment_score", "detachment", "path_annotation", "n_genes"]].copy())
    adata.obsm["spatial"] = spots[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    adata.uns["spatial"] = {
        sample_id: {
            "images": {"hires": hires, "cytassist": cyt_aligned},
            "scalefactors": scalefactors,
        }
    }
    return adata, hires, cyt_aligned


def _clear_embedding_axes(ax) -> None:
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_sample(sample_id: str, vmin: float, vmax: float, path_palette: dict[str, str]) -> None:
    sc.settings.verbosity = 0
    adata, hires, cyt = load_adata(sample_id)
    ng = pd.to_numeric(adata.obs["n_genes"], errors="coerce").dropna()
    ng_vmin = float(ng.quantile(0.01)) if len(ng) else 0.0
    ng_vmax = float(ng.quantile(0.99)) if len(ng) else 1.0

    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.7), constrained_layout=True)
    axes[0, 0].imshow(hires)
    axes[0, 0].set_title("H&E")
    axes[0, 0].axis("off")

    axes[0, 1].imshow(cyt)
    axes[0, 1].set_title("CytAssist")
    axes[0, 1].axis("off")

    sc.pl.spatial(
        adata,
        color="path_annotation",
        img_key="hires",
        basis="spatial",
        ax=axes[0, 2],
        show=False,
        size=1.5,
        title="Path annotation",
        palette=path_palette,
        legend_loc="upper right",
        legend_fontsize=7,
    )
    _clear_embedding_axes(axes[0, 2])

    sc.pl.spatial(
        adata,
        color="detachment_score",
        img_key="hires",
        basis="spatial",
        ax=axes[1, 0],
        show=False,
        vmin=vmin,
        vmax=vmax,
        cmap="YlOrRd",
        size=1.5,
        title="Detachment score",
        colorbar_loc="right",
    )
    _clear_embedding_axes(axes[1, 0])

    sc.pl.spatial(
        adata,
        color="detachment",
        img_key="hires",
        basis="spatial",
        ax=axes[1, 1],
        show=False,
        palette=ASSIGN_PALETTE,
        size=1.5,
        title="Detachment",
        legend_loc="upper right",
        legend_fontsize=7,
    )
    _clear_embedding_axes(axes[1, 1])

    sc.pl.spatial(
        adata,
        color="n_genes",
        img_key="hires",
        basis="spatial",
        ax=axes[1, 2],
        show=False,
        vmin=ng_vmin,
        vmax=ng_vmax,
        cmap="viridis",
        size=1.5,
        title="Unique genes",
        colorbar_loc="right",
    )
    _clear_embedding_axes(axes[1, 2])

    fig.suptitle(sample_id, fontsize=10, fontweight="bold")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_pdf = OUT_DIR / f"{sample_id}.pdf"
    out_png = OUT_DIR / f"{sample_id}.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"wrote {out_pdf}", flush=True)


def main() -> None:
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sample_ids = pd.read_csv(SAMPLE_LIST)["sample_id"].astype(str).tolist()
    vmin, vmax = score_limits(sample_ids)
    path_palette = build_path_palette(sample_ids)
    for sample_id in sample_ids:
        plot_sample(sample_id, vmin, vmax, path_palette)


if __name__ == "__main__":
    main()
