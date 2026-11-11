#!/usr/bin/env python3
"""Key metadata of the paired H&E and CytAssist images, from their TIFF tags only (no pixels are read).

For every TIFF in data/raw_data/paired/<sample_id>/images/*/: format (size, compression, tiling, pyramid), pixel
size, scanner and software tags, dates, and the source of exported images. Pixel size comes from XResolution with
ResolutionUnit (cm or inch), or with the ImageJ unit when ResolutionUnit is none. The Metafer (45xxx) and CytAssist
(65xxx) tags are private and undocumented; their names here are inferred from their values (e.g. 45232 = 20030 on a
20x scan, 45233 = 3450 for a 3.45 um camera pixel). For H&E images, um_per_px_he_quality is the value
he_quality.py uses, to check against um_per_px.

Writes results/he_quality/image_metadata.csv, one row per image.

  python code/figures/he_quality/he_metadata.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from he_quality import MICRONS_PER_PIXEL, RAW_DIR, RESULTS_DIR, SAMPLES

UNIT_UM = {2: 25400.0, 3: 10000.0}  # TIFF ResolutionUnit: inch, cm
# private tag: column name, inferred from the values
PRIVATE_TAGS = {
    45232: "metafer_magnification_x1000",
    45233: "metafer_camera_pixel_nm",
    45241: "metafer_scanner_id",
    45243: "metafer_owner",
    45278: "metafer_scan_start",
    45279: "metafer_scan_end",
    65000: "cytassist_name",
    65001: "cytassist_run_uuid",
    65002: "cytassist_slide_id",
}


def tag(page: tifffile.TiffPage, key) -> str:
    t = page.tags.get(key)
    return "" if t is None else str(t.value).strip()


def pixel_size_um(page: tifffile.TiffPage, ij: dict | None) -> tuple[float, str]:
    res = page.tags.get("XResolution")
    if res is None or not res.value[0]:
        return np.nan, ""
    px_per_unit = res.value[0] / res.value[1]
    unit = int(page.tags["ResolutionUnit"].value) if "ResolutionUnit" in page.tags else 1
    if unit in UNIT_UM:
        return UNIT_UM[unit] / px_per_unit, f"XResolution per {'inch' if unit == 2 else 'cm'}"
    ij_unit = str((ij or {}).get("unit", ""))
    if re.search(r"µm|u00B5m|micron|^um$", ij_unit):
        return 1 / px_per_unit, "XResolution per ImageJ unit (um)"
    return np.nan, "no physical unit"


def image_metadata(path: Path) -> dict:
    with tifffile.TiffFile(path) as tif:
        page = tif.pages[0]
        series = tif.series[0]
        ij = tif.imagej_metadata
        um, um_source = pixel_size_um(page, ij)
        description = tag(page, "ImageDescription")
        info = str((ij or {}).get("Info", ""))
        source = re.search(r"file:(\S+)", info)
        return {
            "file": path.name,
            "file_size_gb": round(path.stat().st_size / 1e9, 3),
            "height_px": series.shape[0],
            "width_px": series.shape[1],
            "dtype": str(series.dtype),
            "compression": page.compression.name,
            "photometric": page.photometric.name,
            "tile_px": f"{page.tilewidth}x{page.tilelength}" if page.is_tiled else "",
            "pyramid_levels": len(series.levels),
            "bigtiff": tif.is_bigtiff,
            "um_per_px": um,
            "um_per_px_source": um_source,
            "make": tag(page, "Make"),
            "model": tag(page, "Model"),
            "software": tag(page, "Software"),
            "datetime": tag(page, "DateTime"),
            "document_name": tag(page, "DocumentName"),
            "image_description": description[:150].replace("\n", " "),
            "imagej_version": (re.match(r"ImageJ=(\S+)", description) or [None, ""])[1],
            "imagej_source": source.group(1) if source else "",
            **{name: tag(page, code) for code, name in PRIVATE_TAGS.items()},
        }


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for _, row in SAMPLES.iterrows():
        sample_id = str(row["sample_id"])
        for path in sorted((RAW_DIR / sample_id / "images").glob("*/*.tif*")):
            image = "he" if path.parent.name == "he_image" else "cytassist"
            rows.append({"sample_id": sample_id, "patient_id": str(row["patient_id"]),
                         "method": sample_id.rsplit("_", 1)[1], "image": image, "folder": path.parent.name,
                         **image_metadata(path),
                         "um_per_px_he_quality": MICRONS_PER_PIXEL[sample_id] if image == "he" else np.nan})
            print(f"{sample_id}: {path.parent.name}/{path.name}", flush=True)
    out = RESULTS_DIR / "image_metadata.csv"
    # drop columns that are empty for every image
    table = pd.DataFrame(rows).replace("", np.nan).dropna(axis=1, how="all")
    table.to_csv(out, index=False)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
