#!/usr/bin/env python3
"""H&E image quality of the paired STAY vs standard sections, from the raw H&E images only (no spots).

For each paired section, reads data/raw_data/paired/<sample_id>/images/he_image/*.tif and:
  1. rescales it to TILE_UM microns per pixel (tile metrics) and MASK_UM (section metrics), so sections scanned
     at different resolutions are measured on the same physical scale
  2. finds the tissue (Otsu on saturation, fragments smaller than MIN_FRAGMENT_MM2 dropped)
  3. per TILE_PX x TILE_PX tile that is at least MIN_TISSUE tissue: sharpness (variance of the Laplacian),
     contrast, haematoxylin and eosin optical density (skimage rgb2hed)
  4. per section: tile medians, worst-decile sharpness, stain unevenness across tiles (CV), and artefacts:
     folds (tissue darker than FOLD_OD_FACTOR x the section's median optical density), white space inside the
     tissue outline (tears, but also fat and lumina), fragments, debris and brightness of the background

There are four patients with one section per protocol, so no p-values: the result is the STAY - standard
difference per patient, and for metrics with a better direction, in how many patients STAY is better.
Colour metrics also depend on the scanner; check the scanner columns of he_quality_sections.csv.

Writes to results/he_quality:
  - he_quality_tiles.csv: one row per tile
  - he_quality_sections.csv: one row per section, with pixel size and scanner metadata from the TIFF tags
  - he_quality_paired.csv: per metric and patient, STAY - standard and which protocol is better
  - he_quality.pdf/png: per-patient paired values of every section metric
  - he_quality_tiles.pdf/png: tile distributions per section (intra-section variation)
  - masks/<sample_id>.png: the section next to its masks (red folds, blue white space, green debris), to check

  python code/figures/he_quality/he_quality.py
"""

from __future__ import annotations

from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.color import rgb2hed, rgb2hsv
from skimage.filters import threshold_otsu
from skimage.morphology import disk

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RAW_DIR = PROJECT_ROOT / "data" / "raw_data" / "paired"
RESULTS_DIR = PROJECT_ROOT / "results" / "he_quality"
SAMPLES = pd.read_csv(PROJECT_ROOT / "resources" / "sample_list.csv", dtype=str).query("cohort == 'paired'")

# pixel size per sample_id, only for images whose TIFF tags don't hold it
MICRONS_PER_PIXEL: dict[str, float] = {}
TILE_UM = 1.0
TILE_PX = 256
MASK_UM = 8.0
MIN_TISSUE = 0.8
MIN_FRAGMENT_MM2 = 0.02
FOLD_OD_FACTOR = 2.0
# per patient, a difference smaller than this fraction of the pair's mean counts as a tie
TIE = 0.05

METHODS = ["standard", "STAY"]
METHOD_COLORS = {"STAY": "#FBB4AE", "standard": "#B3CDE3"}
# section metric: (label, better direction: 1 higher is better, -1 lower is better, 0 no better direction)
METRICS = {
    "sharpness_median": ("Sharpness (median tile)", 1),
    "sharpness_p10": ("Sharpness (worst 10% tiles)", 1),
    "contrast_median": ("Contrast", 1),
    "haematoxylin_median": ("Haematoxylin OD", 0),
    "eosin_median": ("Eosin OD", 0),
    "he_ratio": ("Haematoxylin / eosin", 0),
    "haematoxylin_cv": ("Haematoxylin unevenness (CV)", -1),
    "eosin_cv": ("Eosin unevenness (CV)", -1),
    "fold_fraction": ("Folds (fraction of tissue)", -1),
    "white_space_fraction": ("White space in tissue", -1),
    "n_fragments": ("Tissue fragments", -1),
    "debris_fraction": ("Background debris", -1),
    "background_brightness": ("Background brightness", 1),
    "tissue_area_mm2": ("Tissue area (mm²)", 0),
}
TILE_METRICS = {"sharpness": "Sharpness", "contrast": "Contrast", "haematoxylin": "Haematoxylin OD", "eosin": "Eosin OD"}


def read_he(sample_id: str) -> tuple[np.ndarray, float, dict]:
    path = sorted((RAW_DIR / sample_id / "images" / "he_image").glob("*.tif*"))[0]
    with tifffile.TiffFile(path) as tif:
        tags = {tag.name: tag.value for tag in tif.pages[0].tags}
        img = np.squeeze(tif.series[0].asarray())
    if img.ndim == 3 and img.shape[0] in (3, 4) and img.shape[-1] not in (3, 4):
        img = np.moveaxis(img, 0, -1)
    if img.ndim != 3:
        raise SystemExit(f"{sample_id}: expected an RGB image, got shape {img.shape}")
    img = img[..., :3]
    if img.dtype != np.uint8:
        img = (img.astype(np.float32) / img.max() * 255).astype(np.uint8)

    um = MICRONS_PER_PIXEL.get(sample_id)
    if um is None:
        unit = {2: 25400.0, 3: 10000.0}.get(int(tags.get("ResolutionUnit", 1)))
        res = tags.get("XResolution")
        um = unit * res[1] / res[0] if unit and res and res[0] else None
    # tags often hold a placeholder like 72 dpi; anything outside microscope range is not a real pixel size
    if um is None or not 0.05 < um < 5:
        raise SystemExit(f"{sample_id}: no usable pixel size in the TIFF tags ({um}); add it to MICRONS_PER_PIXEL")
    if um > TILE_UM:
        raise SystemExit(f"{sample_id}: {um:.2f} um/px is coarser than TILE_UM = {TILE_UM}")

    meta = {
        "file": path.name,
        "height_px": img.shape[0],
        "width_px": img.shape[1],
        "um_per_px": um,
        "make": str(tags.get("Make", "")),
        "model": str(tags.get("Model", "")),
        "software": str(tags.get("Software", "")),
        "description": str(tags.get("ImageDescription", ""))[:200].replace("\n", " "),
    }
    return img, um, meta


def rescale(img: np.ndarray, factor: float) -> np.ndarray:
    return cv2.resize(img, None, fx=factor, fy=factor, interpolation=cv2.INTER_AREA)


def section_masks(low: np.ndarray) -> dict[str, np.ndarray]:
    """Masks at MASK_UM resolution."""
    sat = rgb2hsv(low)[..., 1]
    raw = sat > threshold_otsu(sat)
    labels, _ = ndi.label(raw)
    sizes = np.bincount(labels.ravel())
    big = sizes >= MIN_FRAGMENT_MM2 / (MASK_UM / 1000) ** 2
    big[0] = False
    tissue = big[labels]
    outline = ndi.binary_fill_holes(ndi.binary_closing(tissue, structure=disk(3)))
    od = -np.log10((low.astype(np.float32) + 1) / 256).mean(axis=-1)
    return {
        "tissue": tissue,
        "outline": outline,
        "white": outline & ~tissue,
        "fold": tissue & (od > FOLD_OD_FACTOR * np.median(od[tissue])),
        "debris": raw & ~tissue,
        "background": ~ndi.binary_dilation(outline, structure=disk(5)),
    }


def section_metrics(low: np.ndarray, masks: dict[str, np.ndarray]) -> dict:
    gray = cv2.cvtColor(low, cv2.COLOR_RGB2GRAY)
    background = masks["background"]
    return {
        "tissue_area_mm2": float(masks["outline"].sum() * (MASK_UM / 1000) ** 2),
        "n_fragments": int(ndi.label(masks["outline"])[1]),
        "white_space_fraction": float(masks["white"].sum() / masks["outline"].sum()),
        "fold_fraction": float(masks["fold"].sum() / masks["tissue"].sum()),
        "debris_fraction": float(masks["debris"][background].mean()) if background.any() else np.nan,
        "background_brightness": float(np.median(gray[background])) if background.any() else np.nan,
    }


def tile_metrics(img: np.ndarray, tissue: np.ndarray) -> pd.DataFrame:
    rows = []
    for y in range(0, img.shape[0] - TILE_PX + 1, TILE_PX):
        for x in range(0, img.shape[1] - TILE_PX + 1, TILE_PX):
            t = tissue[y : y + TILE_PX, x : x + TILE_PX]
            if t.mean() < MIN_TISSUE:
                continue
            tile = img[y : y + TILE_PX, x : x + TILE_PX]
            gray = cv2.cvtColor(tile, cv2.COLOR_RGB2GRAY)
            hed = rgb2hed(tile)
            rows.append(
                {
                    "x": x,
                    "y": y,
                    "sharpness": float(cv2.Laplacian(gray, cv2.CV_64F).var()),
                    "contrast": float(gray[t].std()),
                    "haematoxylin": float(hed[..., 0][t].mean()),
                    "eosin": float(hed[..., 1][t].mean()),
                }
            )
    return pd.DataFrame(rows)


def tile_summary(tiles: pd.DataFrame) -> dict:
    return {
        "n_tiles": len(tiles),
        "sharpness_median": float(tiles["sharpness"].median()),
        "sharpness_p10": float(tiles["sharpness"].quantile(0.1)),
        "contrast_median": float(tiles["contrast"].median()),
        "haematoxylin_median": float(tiles["haematoxylin"].median()),
        "eosin_median": float(tiles["eosin"].median()),
        "he_ratio": float(tiles["haematoxylin"].median() / tiles["eosin"].median()),
        "haematoxylin_cv": float(tiles["haematoxylin"].std() / tiles["haematoxylin"].mean()),
        "eosin_cv": float(tiles["eosin"].std() / tiles["eosin"].mean()),
    }


def save_mask_check(sample_id: str, low: np.ndarray, masks: dict[str, np.ndarray]) -> None:
    over = low.copy()
    dim = ~masks["outline"]
    over[dim] = over[dim] // 2
    over[masks["white"]] = (0, 120, 255)
    over[masks["fold"]] = (255, 0, 0)
    over[masks["debris"]] = (0, 200, 0)
    (RESULTS_DIR / "masks").mkdir(parents=True, exist_ok=True)
    plt.imsave(RESULTS_DIR / "masks" / f"{sample_id}.png", np.hstack([low, over]))


def analyse(sample_id: str) -> tuple[dict, pd.DataFrame]:
    full, um, meta = read_he(sample_id)
    print(f"{sample_id}: {meta['file']} {meta['width_px']}x{meta['height_px']} px, {um:.3f} um/px, "
          f"{meta['make']} {meta['model']} {meta['software']}".rstrip(), flush=True)
    img = rescale(full, um / TILE_UM)
    del full
    low = rescale(img, TILE_UM / MASK_UM)
    masks = section_masks(low)
    tissue = cv2.resize(masks["tissue"].astype(np.uint8), (img.shape[1], img.shape[0]), interpolation=cv2.INTER_NEAREST).astype(bool)
    tiles = tile_metrics(img, tissue)
    save_mask_check(sample_id, low, masks)
    return {**meta, **section_metrics(low, masks), **tile_summary(tiles)}, tiles


def compare_pairs(sections: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric, (_, better) in METRICS.items():
        wide = sections.pivot(index="patient_id", columns="method", values=metric)
        for patient, row in wide.iterrows():
            diff = row["STAY"] - row["standard"]
            if better == 0 or not np.isfinite(diff):
                winner = "-"
            elif abs(diff) <= TIE * abs(row[METHODS].mean()):
                winner = "tie"
            else:
                winner = "STAY" if np.sign(diff) == better else "standard"
            rows.append({"metric": metric, "patient_id": patient, "standard": row["standard"], "STAY": row["STAY"],
                         "STAY_minus_standard": diff, "better": winner})
    return pd.DataFrame(rows)


def verdict(better: pd.Series) -> str:
    """How many patients each protocol is better in, e.g. "STAY 3, standard 0, tie 1"."""
    if (better == "-").all():
        return "no better direction"
    return ", ".join(f"{k} {(better == k).sum()}" for k in ["STAY", "standard", "tie"])


def plot_sections(sections: pd.DataFrame, pairs: pd.DataFrame) -> None:
    ncol = 5
    nrow = int(np.ceil(len(METRICS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.2 * ncol, 2.3 * nrow), squeeze=False)
    patients = sorted(sections["patient_id"].unique())
    colors = dict(zip(patients, plt.get_cmap("tab10").colors))
    for ax, (metric, (label, _)) in zip(axes.ravel(), METRICS.items()):
        wide = sections.pivot(index="patient_id", columns="method", values=metric)
        for patient, row in wide.iterrows():
            ax.plot([0, 1], row[METHODS], "-o", color=colors[patient], ms=4, lw=1, label=patient)
        ax.set_title(f"{label}\n{verdict(pairs.loc[pairs['metric'] == metric, 'better'])}")
        ax.set_xticks([0, 1], METHODS)
        ax.set_xlim(-0.3, 1.3)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.ravel()[len(METRICS) :]:
        ax.axis("off")
    axes[0, 0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_quality.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_tiles(tiles: pd.DataFrame, sections: pd.DataFrame) -> None:
    order = sections.assign(rank=sections["method"].map({m: i for i, m in enumerate(METHODS)})).sort_values(["patient_id", "rank"])
    fig, axes = plt.subplots(len(TILE_METRICS), 1, figsize=(1.0 * len(order) + 1, 2.0 * len(TILE_METRICS)), sharex=True)
    for ax, (metric, label) in zip(axes, TILE_METRICS.items()):
        data = [tiles.loc[tiles["sample_id"] == s, metric].to_numpy() for s in order["sample_id"]]
        parts = ax.violinplot(data, showmedians=True, showextrema=False)
        for body, method in zip(parts["bodies"], order["method"]):
            body.set_facecolor(METHOD_COLORS[method])
            body.set_edgecolor("black")
            body.set_alpha(0.9)
        parts["cmedians"].set_color("black")
        if metric == "sharpness":
            ax.set_yscale("log")
        ax.set_ylabel(label)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks(np.arange(1, len(order) + 1), order["sample_id"], rotation=45, ha="right")
    fig.tight_layout()
    for ext in ["pdf", "png"]:
        fig.savefig(RESULTS_DIR / f"he_quality_tiles.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "pdf.fonttype": 42})
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    section_rows, tile_frames = [], []
    for _, row in SAMPLES.iterrows():
        sample_id = str(row["sample_id"])
        info = {"sample_id": sample_id, "patient_id": str(row["patient_id"]), "method": sample_id.rsplit("_", 1)[1]}
        metrics, tiles = analyse(sample_id)
        section_rows.append({**info, **metrics})
        tile_frames.append(tiles.assign(**info))
    sections = pd.DataFrame(section_rows)
    tiles = pd.concat(tile_frames, ignore_index=True)
    pairs = compare_pairs(sections)

    tiles.to_csv(RESULTS_DIR / "he_quality_tiles.csv", index=False)
    sections.to_csv(RESULTS_DIR / "he_quality_sections.csv", index=False)
    pairs.to_csv(RESULTS_DIR / "he_quality_paired.csv", index=False)
    plot_sections(sections, pairs)
    plot_tiles(tiles, sections)

    summary = pairs.groupby("metric", sort=False).agg(
        better_in=("better", verdict),
        median_STAY_minus_standard=("STAY_minus_standard", "median"),
    )
    print(summary.to_string(), flush=True)
    print(f"wrote {RESULTS_DIR}", flush=True)


if __name__ == "__main__":
    main()
