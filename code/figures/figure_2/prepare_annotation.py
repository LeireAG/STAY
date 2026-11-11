#!/usr/bin/env python3
"""Light images for annotate_detachment.py, so the annotation can run on another machine than the data.

For each STAY and standard section, writes to data/analysis/{cohort}/detachment_annotation_input/ two JPEGs, both
downscaled by ANNOTATION_SCALE:
  - <sample_id>_cytassist.jpg: the CytAssist image (Space Ranger's spatial/cytassist_image.tiff)
  - <sample_id>_he.jpg: Space Ranger's hires H&E, warped onto the CytAssist image through the spot positions both
    share (Loupe alignment)
These files are tracked by git (see .gitignore): run this where the data is, push, and pull where you annotate.

  python code/figures/figure_2/prepare_annotation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import cv2
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_detachment_pipeline import (ANNOTATION_SCALE, SAMPLES, annotation_input_path, cohort_dir,
                                     read_cytassist_coords, read_cytassist_image, read_image, read_spots, warp_image)

JPEG_QUALITY = 90


def save_jpeg(path: Path, rgb: np.ndarray) -> None:
    small = cv2.resize(rgb, None, fx=ANNOTATION_SCALE, fy=ANNOTATION_SCALE, interpolation=cv2.INTER_AREA)
    path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(path), cv2.cvtColor(small, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, JPEG_QUALITY])


def main() -> None:
    for sample_id in SAMPLES["sample_id"]:
        sample_dir = cohort_dir("spaceranger", sample_id) / sample_id
        cyt = read_cytassist_image(sample_id)
        spots = read_cytassist_coords(sample_id, read_spots(sample_dir))
        hires = read_image(sample_dir / "spatial" / "tissue_hires_image.png")
        he = warp_image(hires, spots[["hires_x", "hires_y"]].to_numpy(float), spots[["cyt_x", "cyt_y"]].to_numpy(float),
                        cyt.shape)
        save_jpeg(annotation_input_path(sample_id, "cytassist"), cyt)
        save_jpeg(annotation_input_path(sample_id, "he"), he)
        print(f"{sample_id}: wrote annotation images", flush=True)


if __name__ == "__main__":
    main()
