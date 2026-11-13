#!/usr/bin/env python3
"""Blinded manual annotation of detached tissue on the CytAssist images, in the browser (Panel/HoloViews).

Shows the STAY and standard sections one at a time in a new random order per browser session, labelled only
"Section k of n" (never the sample ID or cohort): the CytAssist image with the H&E warped onto it and a polygon layer.
Drawing (Bokeh's freehand tool): click and drag around a detached area; releasing closes the outline. Click an outline
and press Backspace to delete it; the vertex-edit tool in the toolbar adjusts its corners.
The H&E opacity slider fades between the two images; "Show H&E" shows the H&E fully. The next section is loaded in
the background. Keys: Enter or S save and next, hold H to show the H&E, Shift + scroll to change the H&E opacity, hold
Space and drag to pan, Z undo the last polygon, B back, K skip.
Polygons are saved in CytAssist image pixels (x right, y down) to
data/analysis/{cohort}/detachment_annotation/<sample_id>.json, which run_detachment_pipeline.py reads; a section saved
with no polygons counts as annotated without detachment.

  panel serve code/figures/figure_2/annotate_detachment.py --port 5006   # on the server
  ssh -N -L 5006:localhost:5006 <user>@<server>                           # on the laptop, then open
  http://localhost:5006/annotate_detachment                               # (add ?review=1 to revisit all)
"""

from __future__ import annotations

import base64
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
JPEG_QUALITY = 95  # images go to the browser at native resolution as JPEG, ~1-3 MB instead of ~36 MB raw
KEYS = ("Keys: **Enter**/**S** save and next · hold **H** show H&E · **Shift + scroll** H&E opacity · "
        "hold **Space** + drag to pan · **Z** undo last polygon · **B** back · **K** skip")


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
    const tool = down ? "PanTool" : "FreehandDrawTool"
    for (const doc of Bokeh.documents) for (const m of doc.all_models) if (m.constructor.__name__ === tool) m.active = true
  } else if (key === "h") data.h = down
  else if (down) data.pressed = key + ":" + Date.now()
}
document.addEventListener("keydown", (e) => report(e, true))
document.addEventListener("keyup", (e) => report(e, false))
document.addEventListener("wheel", (e) => {  // Shift + scroll: H&E opacity slider, instead of zooming
  if (!e.shiftKey) return
  e.preventDefault(); e.stopPropagation()
  for (const doc of Bokeh.documents) for (const m of doc.all_models)
    if (m.constructor.__name__ === "Slider" && m.title === "H&E opacity") {
      const step = (e.deltaY || e.deltaX) < 0 ? 0.05 : -0.05
      m.value = Math.round(Math.min(1, Math.max(0, m.value + step)) * 20) / 20
    }
}, {capture: true, passive: false})"""}


def jpeg_url(rgb: np.ndarray) -> str:
    _, jpeg = cv2.imencode(".jpg", cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, JPEG_QUALITY])
    return "data:image/jpeg;base64," + base64.b64encode(jpeg).decode()


def load(sample_id: str) -> tuple[int, int, str, str]:
    """Height, width and JPEG data URLs of the native CytAssist image and the hires H&E warped onto it."""
    sample_dir = cohort_dir("spaceranger", sample_id) / sample_id
    cyt = read_cytassist_image(sample_id)
    spots = read_cytassist_coords(sample_id, read_spots(sample_dir))
    he = warp_image(read_image(sample_dir / "spatial" / "tissue_hires_image.png"),
                    spots[["hires_x", "hires_y"]].to_numpy(float), spots[["cyt_x", "cyt_y"]].to_numpy(float), cyt.shape)
    return *cyt.shape[:2], jpeg_url(cyt), jpeg_url(he)


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
save, skip, back = (pn.widgets.Button(name=n, button_type=t) for n, t in
                    [("Save and next (Enter)", "primary"), ("Skip (K)", "default"), ("Back (B)", "default")])
status, count, plot, keys = pn.pane.Markdown(), pn.pane.Markdown(), pn.pane.HoloViews(), Keys(width=0, height=0)


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
        status.object, count.object, plot.object = f"**All {len(order)} sections done.**", "", None
        for button in (save, skip):
            button.disabled = True
        return
    h, w, cyt_url, he_url = images(k).result()
    if k + 1 < len(order):
        images(k + 1)  # preload the next section
    path = annotation_path(order[k])
    if drawn is None:
        drawn = json.loads(path.read_text())["polygons"] if path.exists() else []
    # the plot's y axis points up, image pixels point down: y_plot = h - y_pixel
    polys = hv.Polygons([{"x": [x for x, _ in p], "y": [h - y for _, y in p]} for p in drawn]).opts(
        fill_alpha=0.35, fill_color="yellow", line_color="yellow", line_width=3)
    corners = {"size": 7, "fill_color": "red", "line_color": "white"}  # shown by the vertex-edit tool
    state.update(h=h, draw=streams.FreehandDraw(source=polys))
    streams.PolyEdit(source=polys, vertex_style=corners)
    state["draw"].add_subscriber(lambda **_: count.param.update(object=f"**Polygons on this section: {len(polygons())}**"))
    count.object = f"**Polygons on this section: {len(drawn)}**"
    def add_images(plot, _) -> None:  # the two images, under the polygons, as JPEGs the browser decodes itself
        fig = plot.state
        layers = [fig.image_url(url=[url], x=0, y=h, w=w, h=h, anchor="top_left") for url in (cyt_url, he_url)]
        fig.renderers = layers + [r for r in fig.renderers if r not in layers]
        state["he_glyph"] = layers[1].glyph
        set_he_alpha()

    plot.object = polys.opts(xlim=(0, w), ylim=(0, h), data_aspect=1, frame_width=1000, xaxis=None, yaxis=None,
                             active_tools=["freehand_draw", "wheel_zoom"], hooks=[add_images])
    done = "already saved, " if path.exists() else ""
    status.object = (f"**Section {k + 1} of {len(order)}** ({done}{state['saved']} saved this session). "
                     f"**Click and drag** around a detached area to outline it; click an outline + Backspace to delete.  "
                     f"\n{KEYS}")
    for button in (save, skip):
        button.disabled = False
    back.disabled = k == 0


def set_he_alpha(*_) -> None:
    if "he_glyph" in state:
        state["he_glyph"].global_alpha = 1.0 if show_he.value else opacity.value


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
opacity.param.watch(set_he_alpha, "value")
show_he.param.watch(set_he_alpha, "value")
keys.param.watch(lambda event: setattr(show_he, "value", event.new), "h")
show()
pn.Column(status, pn.Row(opacity, show_he, save, skip, back, count, keys), plot).servable(title="Detachment annotation")
