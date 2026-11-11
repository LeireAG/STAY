#!/usr/bin/env python3
"""Blinded manual annotation of detached tissue on the CytAssist images, in the browser (Panel/HoloViews).

Shows the STAY and standard sections one at a time in a new random order per browser session, labelled only
"Section k of n" (never the sample ID or cohort): the CytAssist image with the H&E warped onto it and a polygon layer.
Drawing (Bokeh's polygon tool): press and hold to place the first corner, click to add corners, press and hold (or
Esc) to place the last one; drag a polygon to move it; click it and press Backspace to delete it.
The H&E opacity slider fades between the two images; "Show H&E" shows the H&E fully. The next section is loaded in
the background. Keys: Enter or S save and next, hold H to show the H&E, hold Space and drag to pan, Z undo the last
polygon, B back, K skip.
Polygons are saved in CytAssist image pixels (x right, y down) to
data/analysis/{cohort}/detachment_annotation/<sample_id>.json, which run_detachment_pipeline.py reads; a section saved
with no polygons counts as annotated without detachment.

  panel serve code/figures/figure_2/annotate_detachment.py --port 5006   # on the server
  ssh -N -L 5006:localhost:5006 <user>@<server>                           # on the laptop, then open
  http://localhost:5006/annotate_detachment                               # (add ?review=1 to revisit all)
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import cv2
import holoviews as hv
import numpy as np
import panel as pn
import param
from holoviews import streams

if __name__ == "__main__":  # plain `python` only builds the page; `panel serve` serves it
    sys.exit("Start it with: panel serve code/figures/figure_2/annotate_detachment.py --port 5006")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_detachment_pipeline import (SAMPLES, annotation_path, cohort_dir, read_cytassist_coords,
                                     read_cytassist_image, read_image, read_spots, warp_image)

hv.extension("bokeh")
DISPLAY_PX = 1200  # images are sent to the browser at this size; polygons stay in full-size pixels
KEYS = ("Keys: **Enter**/**S** save and next · hold **H** show H&E · hold **Space** + drag to pan · "
        "**Z** undo last polygon · **B** back · **K** skip")


class Keys(pn.reactive.ReactiveHTML):
    """Shortcut keys anywhere on the page (not while typing in a text field). Space switches Bokeh's pan tool on while
    held, in the browser. H is held, so it gets its own on/off field: a release can never be lost behind the next key."""

    h = param.Boolean()
    pressed = param.String()  # "<key>:<time>" for the other keys, so repeating a key still counts as a change
    _template = "<div></div>"
    _scripts = {"render": """
const report = (e, down) => {
  const key = e.key === " " ? "space" : e.key.toLowerCase()
  if (!["enter", "s", "h", "space", "z", "b", "k"].includes(key) || e.ctrlKey || e.metaKey
      || ["INPUT", "TEXTAREA"].includes(e.target.tagName)) return
  if (key === "space") e.preventDefault()  // no page scrolling
  if (e.repeat) return
  if (key === "space") {
    const tool = down ? "PanTool" : "PolyDrawTool"
    for (const doc of Bokeh.documents) for (const m of doc.all_models) if (m.constructor.__name__ === tool) m.active = true
  } else if (key === "h") data.h = down
  else if (down) data.pressed = key + ":" + Date.now()
}
document.addEventListener("keydown", (e) => report(e, true))
document.addEventListener("keyup", (e) => report(e, false))"""}


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
state = {"k": 0, "saved": 0}
pool, loaded = ThreadPoolExecutor(max_workers=1), {}


def images(k: int):
    """Section k's images, loading them in the background if not requested yet."""
    if k not in loaded:
        loaded[k] = pool.submit(load, order[k])
    return loaded[k]


opacity = pn.widgets.FloatSlider(name="H&E opacity", start=0, end=1, step=0.05, value=0)
show_he = pn.widgets.Toggle(name="Show H&E (H)")
he_alpha = pn.bind(lambda show, value: 1.0 if show else value, show_he, opacity)
save, skip, back = (pn.widgets.Button(name=n, button_type=t) for n, t in
                    [("Save and next (Enter)", "primary"), ("Skip (K)", "default"), ("Back (B)", "default")])
status, plot, keys = pn.pane.Markdown(), pn.pane.HoloViews(), Keys(width=0, height=0)


def polygons() -> list[list[list[float]]]:
    """The polygons drawn so far, in image pixels (y down), without the closing repeat of the first corner."""
    data, h = state["draw"].data or {"xs": [], "ys": []}, state["h"]
    out = []
    for xs, ys in zip(data["xs"], data["ys"]):
        points = [[round(float(x), 1), round(h - float(y), 1)] for x, y in zip(xs, ys)]
        if len(points) > 1 and points[0] == points[-1]:
            points.pop()
        if len(points) >= 3:
            out.append(points)
    return out


def show(drawn: list | None = None) -> None:
    """Display section k; with drawn, redisplay it with those polygons instead of the saved ones."""
    k = state["k"]
    if k >= len(order):
        status.object, plot.object = f"**All {len(order)} sections done.**", None
        for button in (save, skip):
            button.disabled = True
        return
    h, w, cyt, he = images(k).result()
    if k + 1 < len(order):
        images(k + 1)  # preload the next section
    path = annotation_path(order[k])
    if drawn is None:
        drawn = json.loads(path.read_text())["polygons"] if path.exists() else []
    # the plot's y axis points up, image pixels point down: y_plot = h - y_pixel
    polys = hv.Polygons([{"x": [x for x, _ in p], "y": [h - y for _, y in p]} for p in drawn]).opts(
        fill_alpha=0.2, fill_color="yellow", line_color="yellow", line_width=2)
    state.update(h=h, draw=streams.PolyDraw(source=polys, drag=True, show_vertices=True))
    streams.PolyEdit(source=polys)
    bounds = (0, 0, w, h)
    plot.object = (hv.RGB(cyt, bounds=bounds) * hv.RGB(he, bounds=bounds).apply.opts(alpha=he_alpha) * polys).opts(
        data_aspect=1, frame_width=1000, xaxis=None, yaxis=None, active_tools=["poly_draw", "wheel_zoom"])
    done = "already saved, " if path.exists() else ""
    status.object = (f"**Section {k + 1} of {len(order)}** ({done}{state['saved']} saved this session). "
                     f"**Hold** to start a polygon, click corners, **hold** (or Esc) to finish; drag to move, click + Backspace "
                     f"to delete.  \n{KEYS}")
    for button in (save, skip):
        button.disabled = False
    back.disabled = k == 0


def go(step: int, save_current: bool = False) -> None:
    if save_current:
        path = annotation_path(order[state["k"]])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"polygons": polygons()}))
        state["saved"] += 1
    state["k"] = max(0, state["k"] + step)
    show()


def on_key(event) -> None:
    key = event.new.split(":")[0]
    if key in ("enter", "s") and not save.disabled:
        go(+1, save_current=True)
    elif key == "k" and not skip.disabled:
        go(+1)
    elif key == "b" and not back.disabled:
        go(-1)
    elif key == "z" and state["k"] < len(order):
        show(drawn=polygons()[:-1])


save.on_click(lambda _: go(+1, save_current=True))
skip.on_click(lambda _: go(+1))
back.on_click(lambda _: go(-1))
keys.param.watch(on_key, "pressed")
keys.param.watch(lambda event: setattr(show_he, "value", event.new), "h")
show()
pn.Column(status, pn.Row(opacity, show_he, save, skip, back, keys), plot).servable(title="Detachment annotation")
