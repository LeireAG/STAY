#!/usr/bin/env python3
"""Blinded manual annotation of detached tissue on the CytAssist images, in napari.

Needs only the light images of prepare_annotation.py (data/analysis/{cohort}/detachment_annotation_input/, tracked
by git), so it can run on a laptop without the Space Ranger data. Opens the STAY and standard sections one at a time,
in a new random order on every run. Neither the window nor the terminal ever shows the sample ID or cohort, only
"Section k of n". Layers:
  - CytAssist: the CytAssist image
  - H&E: the H&E warped onto the CytAssist image. Change its opacity in the layer controls to fade between the H&E
    (before) and the CytAssist image.
  - Detached: draw a polygon around each detached area (polygon or rectangle tools; ellipses and lines are not
    saved). Earlier polygons of the section are loaded for editing.
Close the window when done; the terminal then asks whether to save. Polygons are saved, in full-size CytAssist image
pixels (x right, y down; the images shown are downscaled by ANNOTATION_SCALE), to
data/analysis/{cohort}/detachment_annotation/<sample_id>.json, which run_detachment_pipeline.py reads. A section saved
with no polygons counts as annotated without detachment.

  python code/figures/figure_2/annotate_detachment.py            # sections not annotated yet
  python code/figures/figure_2/annotate_detachment.py --review   # every section, to check or correct
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import napari
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_detachment_pipeline import ANNOTATION_SCALE, SAMPLES, annotation_input_path, annotation_path, read_image

SAVED_SHAPES = {"polygon", "rectangle"}


def ask(question: str, answers: str) -> str:
    while (answer := input(question).strip().lower()) not in answers:
        pass
    return answer


def annotate(sample_id: str, k: int, n: int) -> bool:
    """Annotate one section; False when the annotator wants to stop."""
    path = annotation_path(sample_id)
    earlier = json.loads(path.read_text())["polygons"] if path.exists() else []

    viewer = napari.Viewer(title=f"Section {k} of {n}")
    viewer.add_image(read_image(annotation_input_path(sample_id, "cytassist")), name="CytAssist")
    viewer.add_image(read_image(annotation_input_path(sample_id, "he")), name="H&E", opacity=0.0)
    # napari shapes are (row, col) = (y, x), on the downscaled image
    shown = [np.asarray(p)[:, ::-1] * ANNOTATION_SCALE for p in earlier]
    layer = viewer.add_shapes(shown or None, shape_type="polygon", name="Detached",
                              edge_color="yellow", edge_width=3, face_color=[1, 1, 0, 0.15])
    layer.mode = "add_polygon"
    print(f"Section {k} of {n}: draw the detached areas, then close the window", flush=True)
    napari.run()

    polygons = [(np.asarray(d)[:, ::-1] / ANNOTATION_SCALE).round(1).tolist() for d, t in zip(layer.data, layer.shape_type)
                if t in SAVED_SHAPES and len(d) >= 3]
    skipped = sum(t not in SAVED_SHAPES for t in layer.shape_type)
    if skipped:
        print(f"  {skipped} shape(s) that are not polygons or rectangles will not be saved", flush=True)
    answer = ask(f"  {len(polygons)} polygon(s): [s]ave, [d]iscard, or [q]uit without saving? ", "sdq")
    if answer == "s":
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"polygons": polygons}))
    return answer != "q"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--review", action="store_true", help="go through every section, not only unannotated ones")
    args = parser.parse_args()

    missing = sum(not annotation_input_path(s, "cytassist").exists() for s in SAMPLES["sample_id"])
    if missing:
        raise SystemExit(f"{missing} sections have no annotation images; run prepare_annotation.py where the data "
                         "is, push, and pull here")
    samples = [s for s in SAMPLES["sample_id"] if args.review or not annotation_path(s).exists()]
    order = np.random.default_rng().permutation(samples)  # unseeded: a new order on every run
    for k, sample_id in enumerate(order, start=1):
        if not annotate(sample_id, k, len(order)):
            break
    done = sum(annotation_path(s).exists() for s in SAMPLES["sample_id"])
    print(f"{done} of {len(SAMPLES)} sections annotated", flush=True)


if __name__ == "__main__":
    main()
