# Intermediate analysis results: STAY

STAY cohort of the detachment analysis (Figure 2).

## Layout

- `detachment_annotation_input/<sample_id>_cytassist.jpg`, `<sample_id>_he.jpg`: CytAssist image and H&E warped onto it, downscaled 2×, from `prepare_annotation.py` (tracked by git).
- `detachment_annotation/<sample_id>.json`: manually annotated detached areas, as polygons in CytAssist image pixels, from `annotate_detachment.py`.
- `detachment/<sample_id>.csv`: per-spot calls (`barcode`, `in_tissue`, `detached`, `tissue`), and `detachment/detached_area.csv`: spot counts and detached area per section, from `run_detachment_pipeline.py`.

Sample IDs: `P01` ... `P32` (`P29_1`, `P29_2`: two sections of P29). See `resources/sample_list.csv`.

## Source

- https://doi.org/10.5281/zenodo.23558413
