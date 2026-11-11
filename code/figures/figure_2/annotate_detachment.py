#!/usr/bin/env python3
"""Blinded manual annotation of detached tissue on the CytAssist images, in the browser (Panel/HoloViews).

Shows the STAY and standard sections one at a time in a new random order per browser session, labelled only
"Section k of n" (never the sample ID or cohort): the CytAssist image with the H&E warped onto it (H&E opacity slider
to fade between them) and a polygon layer. Polygons are saved in CytAssist image pixels (x right, y down) to
data/analysis/{cohort}/detachment_annotation/<sample_id>.json, which run_detachment_pipeline.py reads; a section saved
with no polygons counts as annotated without detachment.

  panel serve code/figures/figure_2/annotate_detachment.py --port 5006   # on the server
  ssh -N -L 5006:localhost:5006 <user>@<server>                           # on the laptop, then open
  http://localhost:5006/annotate_detachment                               # (add ?review=1 to revisit all)
"""

import json
import sys
from pathlib import Path

import cv2
import holoviews as hv
import numpy as np
import panel as pn
from holoviews import streams

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_detachment_pipeline import (SAMPLES, annotation_path, cohort_dir, read_cytassist_coords,
                                     read_cytassist_image, read_image, read_spots, warp_image)

hv.extension("bokeh")
DISPLAY_PX = 1200  # images are sent to the browser at this size; polygons stay in full-size pixels


def load(sample_id: str) -> tuple[int, int, np.ndarray, np.ndarray]:
    """Height, width and display-size CytAssist and H&E (warped onto the CytAssist image through the spots)."""
    sample_dir = cohort_dir("spaceranger", sample_id) / sample_id
    cyt = read_cytassist_image(sample_id)
    spots = read_cytassist_coords(sample_id, read_spots(sample_dir))
    he = warp_image(read_image(sample_dir / "spatial" / "tissue_hires_image.png"),
                    spots[["hires_x", "hires_y"]].to_numpy(float), spots[["cyt_x", "cyt_y"]].to_numpy(float), cyt.shape)
    f = DISPLAY_PX / max(cyt.shape[:2])
    small = [cv2.resize(img, None, fx=f, fy=f, interpolation=cv2.INTER_AREA) for img in (cyt, he)]
    return *cyt.shape[:2], *small


review = "review" in pn.state.session_args
order = np.random.default_rng().permutation(  # unseeded: a new order per session
    [s for s in SAMPLES["sample_id"] if review or not annotation_path(s).exists()]).tolist()
state = {"k": 0}
opacity = pn.widgets.FloatSlider(name="H&E opacity", start=0, end=1, step=0.05, value=0)
save, skip = pn.widgets.Button(name="Save and next", button_type="primary"), pn.widgets.Button(name="Skip")
status, plot = pn.pane.Markdown(), pn.pane.HoloViews()


def show() -> None:
    k = state["k"]
    if k >= len(order):
        status.object, plot.object, save.disabled, skip.disabled = "**All sections done.**", None, True, True
        return
    sample_id = order[k]
    h, w, cyt, he = load(sample_id)
    path = annotation_path(sample_id)
    earlier = json.loads(path.read_text())["polygons"] if path.exists() else []
    # the plot's y axis points up, image pixels point down: y_plot = h - y_pixel
    polys = hv.Polygons([{"x": [x for x, _ in p], "y": [h - y for _, y in p]} for p in earlier]).opts(
        fill_alpha=0.2, fill_color="yellow", line_color="yellow", line_width=2)
    state.update(sample_id=sample_id, h=h, draw=streams.PolyDraw(source=polys, drag=True, show_vertices=True))
    streams.PolyEdit(source=polys)
    bounds = (0, 0, w, h)
    plot.object = (hv.RGB(cyt, bounds=bounds) * hv.RGB(he, bounds=bounds).apply.opts(alpha=opacity) * polys).opts(
        data_aspect=1, frame_width=900, xaxis=None, yaxis=None, active_tools=["poly_draw", "wheel_zoom"])
    status.object = (f"**Section {k + 1} of {len(order)}**: click to add corners, double-click to close a polygon; "
                     "drag to move it; select and press Backspace to delete")


def next_section(saved: bool) -> None:
    if saved:
        data, h = state["draw"].data, state["h"]
        polygons = []
        for xs, ys in zip(data["xs"], data["ys"]):
            points = [[round(float(x), 1), round(h - float(y), 1)] for x, y in zip(xs, ys)]
            if len(points) > 1 and points[0] == points[-1]:  # drop the closing repeat of the first corner
                points.pop()
            if len(points) >= 3:
                polygons.append(points)
        path = annotation_path(state["sample_id"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"polygons": polygons}))
    state["k"] += 1
    show()


save.on_click(lambda _: next_section(True))
skip.on_click(lambda _: next_section(False))
show()
pn.Column(status, pn.Row(opacity, save, skip), plot).servable(title="Detachment annotation")
