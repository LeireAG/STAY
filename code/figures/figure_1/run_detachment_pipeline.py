#!/usr/bin/env python3
"""Score per-spot detachment from images and write CSVs.

Reads config.yaml. Uses the raw barcode set, drops empty-tissue spots,
writes data/detachment_score/{sample_id}.csv with barcode and detachment_score.

  python scripts/run_detachment_pipeline.py
"""

from __future__ import annotations

import json
from pathlib import Path

import cv2
import h5py
import numpy as np
import pandas as pd
import yaml
from skimage.feature import graycomatrix, graycoprops

cv2.utils.logging.setLogLevel(cv2.utils.logging.LOG_LEVEL_ERROR)

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CFG = yaml.safe_load((PROJECT_ROOT / "resources" / "detachment.yaml").read_text(encoding="utf-8"))

SAMPLE_LIST = PROJECT_ROOT / CFG["paths"]["sample_list"]
SPACERANGER_DIR = PROJECT_ROOT / CFG["paths"]["spaceranger"]
RAW_DIR = PROJECT_ROOT / CFG["paths"]["raw"]
OUT_DIR = PROJECT_ROOT / CFG["paths"]["detachment_score"]

MATRIX_NAME = "filtered_feature_bc_matrix.h5" if CFG["scoring"]["matrix"] == "filtered" else "raw_feature_bc_matrix.h5"
PATCH_SIZE_CYTASSIST = int(CFG["scoring"]["patch_size_cytassist"])
PATCH_SIZE_HIRES = int(CFG["scoring"]["patch_size_hires"])
EMPTY_TISSUE_EXPR = str(CFG["scoring"]["empty_tissue_expr"])
DETACHMENT_SCORE_EXPR = str(CFG["scoring"]["detachment_score_expr"])


def read_barcodes(sample_dir: Path) -> list[str]:
    with h5py.File(sample_dir / MATRIX_NAME, "r") as f:
        raw = f["matrix/barcodes"][:]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]


def read_spots(sample_dir: Path) -> pd.DataFrame:
    tp = pd.read_csv(sample_dir / "spatial/tissue_positions.csv")
    tp["barcode"] = tp["barcode"].astype(str)
    keep = set(read_barcodes(sample_dir))
    tp = tp.loc[tp["barcode"].isin(keep)].copy()
    with (sample_dir / "spatial/scalefactors_json.json").open() as handle:
        scale = float(json.load(handle)["tissue_hires_scalef"])
    tp["hires_x"] = pd.to_numeric(tp["pxl_col_in_fullres"]) * scale
    tp["hires_y"] = pd.to_numeric(tp["pxl_row_in_fullres"]) * scale
    return tp


def read_cytassist_coords(sample_id: str, spots: pd.DataFrame) -> pd.DataFrame:
    json_path = sorted((RAW_DIR / sample_id / "manual_alignment").glob("*.json"))[0]
    with json_path.open() as handle:
        oligo = json.load(handle)["oligo"]
    cyt = pd.DataFrame(
        {"array_row": d["row"], "array_col": d["col"], "cyt_x": d["imageX"], "cyt_y": d["imageY"]}
        for d in oligo
    )
    return spots.merge(cyt, on=["array_row", "array_col"], how="left")


def read_image(path: Path) -> np.ndarray:
    img = cv2.imread(str(path), cv2.IMREAD_UNCHANGED)
    if img.ndim == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    if img.dtype != np.uint8:
        img = (img * 255).astype(np.uint8) if img.max() <= 1.0 else img.astype(np.uint8)
    return img


def patch_asm(gray: np.ndarray) -> float:
    patch_q = (gray.astype(np.float32) / 16.0).astype(np.uint8)
    glcm = graycomatrix(
        patch_q,
        distances=[1],
        angles=[0, np.pi / 4, np.pi / 2, 3 * np.pi / 4],
        levels=16,
        symmetric=True,
        normed=True,
    )
    return float(np.mean(graycoprops(glcm, "ASM")))


def patch_stats(img: np.ndarray, xs, ys, patch_size: int, scale: float = 1.0, asm: bool = False) -> pd.DataFrame:
    half_lo = patch_size // 2
    half_hi = patch_size - half_lo
    rows = []
    for x, y in zip(xs, ys):
        if not (np.isfinite(x) and np.isfinite(y)):
            rows.append((0.0, 0.0, 0.0, 0.0, 0.0) if asm else (0.0, 0.0, 0.0, 0.0))
            continue
        ix, iy = int(x * scale), int(y * scale)
        patch = img[iy - half_lo : iy + half_hi, ix - half_lo : ix + half_hi]
        if patch.size == 0:
            rows.append((0.0, 0.0, 0.0, 0.0, 0.0) if asm else (0.0, 0.0, 0.0, 0.0))
            continue
        med = np.median(patch, axis=(0, 1))
        gray = patch if patch.ndim == 2 else cv2.cvtColor(patch, cv2.COLOR_RGB2GRAY)
        brightness = float(np.mean(gray))
        if asm:
            rows.append((float(med[0]), float(med[1]), float(med[2]), brightness, patch_asm(gray)))
        else:
            rows.append((float(med[0]), float(med[1]), float(med[2]), brightness))
    cols = ["r", "g", "b", "brightness"] + (["ASM"] if asm else [])
    return pd.DataFrame(rows, columns=cols)


def align_dynamic_range(source: pd.Series, target: pd.Series) -> pd.Series:
    s_min, s_max = float(source.quantile(0.01)), float(source.quantile(0.99))
    t_min, t_max = float(target.quantile(0.01)), float(target.quantile(0.99))
    if s_max == s_min:
        return source.copy()
    return (source - s_min) / (s_max - s_min) * (t_max - t_min) + t_min


def score_sample(sample_id: str) -> pd.DataFrame:
    sample_dir = SPACERANGER_DIR / sample_id
    spots = read_cytassist_coords(sample_id, read_spots(sample_dir))

    spatial = sample_dir / "spatial"
    cyt_path = spatial / "cytassist_image.tiff"
    if not cyt_path.exists():
        cyt_path = spatial / "cytassist_image.tif"
    cyt_img = read_image(cyt_path)
    hires_img = read_image(spatial / "tissue_hires_image.png")

    cyt = patch_stats(cyt_img, spots["cyt_x"], spots["cyt_y"], PATCH_SIZE_CYTASSIST)
    hires = patch_stats(hires_img, spots["hires_x"], spots["hires_y"], PATCH_SIZE_HIRES, asm=True)

    valid = (cyt["brightness"] != 0) & (hires["brightness"] != 0)
    hires = hires.copy()
    hires.loc[~valid, :] = np.nan
    cyt.loc[~valid, :] = np.nan

    obs = pd.DataFrame(
        {
            "dr_cyt_r": align_dynamic_range(cyt["r"], hires["r"]),
            "dr_cyt_g": align_dynamic_range(cyt["g"], hires["g"]),
            "dr_cyt_b": align_dynamic_range(cyt["b"], hires["b"]),
            "dr_hires_g": hires["g"],
            "dr_hires_ASM": hires["ASM"],
        }
    )
    obs["delta_g"] = obs["dr_cyt_g"] - obs["dr_hires_g"]

    empty = obs.eval(EMPTY_TISSUE_EXPR).fillna(False).astype(bool)
    keep = valid.to_numpy() & ~empty.to_numpy()

    ns = {col: obs.loc[keep, col] for col in obs.columns}
    ns["np"] = np
    scores = pd.to_numeric(eval(DETACHMENT_SCORE_EXPR, ns), errors="coerce").replace([np.inf, -np.inf], np.nan)

    return pd.DataFrame(
        {
            "barcode": spots.loc[keep, "barcode"].to_numpy(),
            "detachment_score": scores.to_numpy(dtype=np.float32),
        }
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(SAMPLE_LIST)
    for sample_id in samples["sample_id"].astype(str):
        out = score_sample(sample_id)
        path = OUT_DIR / f"{sample_id}.csv"
        out.to_csv(path, index=False)
        print(f"wrote {path} ({len(out)} spots)", flush=True)


if __name__ == "__main__":
    main()
