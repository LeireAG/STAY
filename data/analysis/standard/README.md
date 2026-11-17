# Intermediate analysis results: standard

Standard-protocol cohort of the detachment analysis (Figure 2).

## Layout

- `detachment_annotation/<sample_id>.json`: manually annotated detached areas, as polygons in CytAssist image pixels, from `annotate_detachment.py`.
- `detachment/<sample_id>.csv`: per-spot calls (`barcode`, `in_tissue`, `detached`, `tissue`), and `detachment/detached_area.csv`: spot counts and detached area per section, from `run_detachment_pipeline.py`.

Sample IDs: `S01` ... `S20`. See `resources/sample_list.csv`.

## Source

- https://doi.org/10.5281/zenodo.23558413
