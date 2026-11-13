#!/usr/bin/env python3
"""Per-spot detachment calls from the manual annotations made with annotate_detachment.py.

For each STAY and standard section in resources/sample_list.csv, fills the annotated polygons (CytAssist image
pixels) into a mask and writes data/analysis/{cohort}/detachment/{sample_id}.csv with, per spot of the capture area:
  - in_tissue: Space Ranger's tissue call
  - detached: the spot centre (CytAssist coordinates from the Loupe alignment) is inside an annotated polygon
  - tissue: in_tissue or detached, because Space Ranger's tissue mask can leave out detached tissue
and data/analysis/{cohort}/detachment/detached_area.csv with one row per section: tissue and detached spot counts and
the annotated area in mm2. Microns per CytAssist pixel come from the Visium spot pitch (SPOT_PITCH_UM between
neighbouring spots of a row, which are two array columns apart).

  python code/figures/figure_2/run_detachment_pipeline.py
"""

from __future__ import annotations

import json
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CFG = yaml.safe_load((PROJECT_ROOT / "resources" / "detachment.yaml").read_text(encoding="utf-8"))

SAMPLE_LIST = PROJECT_ROOT / CFG["paths"]["sample_list"]
# Figure 2 compares the STAY and standard cohorts; sample_list.csv also lists the other cohorts
SAMPLES = pd.read_csv(SAMPLE_LIST, dtype=str).query("cohort in ['STAY', 'standard']")
COHORT = SAMPLES.set_index("sample_id")["cohort"].to_dict()

# raw matrix: has every spot of the capture area, including detached ones outside Space Ranger's tissue mask
MATRIX_NAME = "raw_feature_bc_matrix.h5"
SPOT_PITCH_UM = 100.0  # Visium centre-to-centre distance of neighbouring spots


def cohort_dir(key: str, sample_id: str) -> Path:
    """Folder of the sample's cohort, e.g. data/spaceranger/STAY, from the {cohort} path template."""
    return PROJECT_ROOT / CFG["paths"][key].format(cohort=COHORT[sample_id])


def annotation_path(sample_id: str) -> Path:
    return cohort_dir("detachment_annotation", sample_id) / f"{sample_id}.json"


def read_spots(sample_dir: Path) -> pd.DataFrame:
    """All spots of the capture area with Space Ranger's tissue call and their position on the hires H&E."""
    tp = pd.read_csv(sample_dir / "spatial/tissue_positions.csv", dtype={"barcode": str})
    with (sample_dir / "spatial/scalefactors_json.json").open() as handle:
        scale = float(json.load(handle)["tissue_hires_scalef"])
    tp["hires_x"] = pd.to_numeric(tp["pxl_col_in_fullres"]) * scale
    tp["hires_y"] = pd.to_numeric(tp["pxl_row_in_fullres"]) * scale
    return tp


def read_cytassist_coords(sample_id: str, spots: pd.DataFrame) -> pd.DataFrame:
    """Adds each spot's position on the CytAssist image (cyt_x, cyt_y) from the Loupe alignment."""
    json_path = sorted((cohort_dir("raw", sample_id) / sample_id / "manual_alignment").glob("*.json"))[0]
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


def read_cytassist_image(sample_id: str) -> np.ndarray:
    spatial = cohort_dir("spaceranger", sample_id) / sample_id / "spatial"
    path = spatial / "cytassist_image.tiff"
    return read_image(path if path.exists() else spatial / "cytassist_image.tif")


def warp_image(img: np.ndarray, src_xy: np.ndarray, dst_xy: np.ndarray, out_shape: tuple[int, ...]) -> np.ndarray:
    """Warp img onto another image through matching spot positions (affine fit, RANSAC)."""
    valid = np.isfinite(src_xy).all(axis=1) & np.isfinite(dst_xy).all(axis=1)
    M, _ = cv2.estimateAffine2D(np.float32(src_xy[valid]), np.float32(dst_xy[valid]), method=cv2.RANSAC,
                                ransacReprojThreshold=5.0)
    return cv2.warpAffine(img, M, (int(out_shape[1]), int(out_shape[0])), flags=cv2.INTER_LINEAR,
                          borderValue=(255, 255, 255))


def cytassist_um_per_px(spots: pd.DataFrame) -> float:
    """Microns per CytAssist pixel: SPOT_PITCH_UM over the median distance between row neighbours."""
    a = spots[["array_row", "array_col", "cyt_x", "cyt_y"]].dropna()
    pairs = a.merge(a.assign(array_col=a["array_col"] - 2), on=["array_row", "array_col"], suffixes=("", "_next"))
    distance = np.hypot(pairs["cyt_x_next"] - pairs["cyt_x"], pairs["cyt_y_next"] - pairs["cyt_y"])
    return SPOT_PITCH_UM / float(distance.median())


def outline_mask(polygons: list, height: int, width: int) -> np.ndarray:
    """Union of the filled outlines. Each outline is filled on its own and from its outer edge, so overlapping
    outlines count once and a self-crossing outline has no holes (filling all at once would leave overlaps empty)."""
    mask = np.zeros((height, width), np.uint8)
    for polygon in polygons:
        one = np.zeros_like(mask)
        cv2.fillPoly(one, [np.round(np.asarray(polygon)).astype(np.int32)], 1)
        outer, _ = cv2.findContours(one, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        cv2.drawContours(mask, outer, -1, 1, thickness=cv2.FILLED)
    return mask


def detachment_calls(sample_id: str) -> tuple[pd.DataFrame, dict]:
    sample_dir = cohort_dir("spaceranger", sample_id) / sample_id
    spots = read_cytassist_coords(sample_id, read_spots(sample_dir))
    height, width = read_cytassist_image(sample_id).shape[:2]
    polygons = json.loads(annotation_path(sample_id).read_text())["polygons"]

    mask = outline_mask(polygons, height, width)
    x = spots["cyt_x"].round().to_numpy()
    y = spots["cyt_y"].round().to_numpy()
    on_image = np.isfinite(x) & np.isfinite(y) & (x >= 0) & (x < width) & (y >= 0) & (y < height)
    detached = np.zeros(len(spots), bool)
    detached[on_image] = mask[y[on_image].astype(int), x[on_image].astype(int)] == 1

    calls = pd.DataFrame({"barcode": spots["barcode"], "in_tissue": spots["in_tissue"].astype(bool),
                          "detached": detached})
    calls["tissue"] = calls["in_tissue"] | calls["detached"]
    um = cytassist_um_per_px(spots)
    area = {
        "sample_id": sample_id,
        "n_tissue_spots": int(calls["tissue"].sum()),
        "n_detached_spots": int(calls["detached"].sum()),
        "n_detached_spots_outside_tissue_mask": int((calls["detached"] & ~calls["in_tissue"]).sum()),
        "detached_area_mm2": float(mask.sum()) * um**2 / 1e6,
        "cytassist_um_per_px": um,
    }
    return calls, area


def main() -> None:
    missing = [s for s in SAMPLES["sample_id"] if not annotation_path(s).exists()]
    if missing:
        # counts only: listing sample IDs would unblind the annotation
        raise SystemExit(f"{len(missing)} of {len(SAMPLES)} sections are not annotated yet; run annotate_detachment.py")

    areas = []
    for sample_id in SAMPLES["sample_id"]:
        calls, area = detachment_calls(sample_id)
        out_dir = cohort_dir("detachment", sample_id)
        out_dir.mkdir(parents=True, exist_ok=True)
        calls.to_csv(out_dir / f"{sample_id}.csv", index=False)
        areas.append({**area, "cohort": COHORT[sample_id]})
        print(f"{sample_id}: {area['n_detached_spots']} of {area['n_tissue_spots']} tissue spots detached, "
              f"{area['detached_area_mm2']:.3f} mm2", flush=True)
    for cohort, rows in pd.DataFrame(areas).groupby("cohort"):
        path = PROJECT_ROOT / CFG["paths"]["detachment"].format(cohort=cohort) / "detached_area.csv"
        rows.drop(columns="cohort").to_csv(path, index=False)
        print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
