# Intermediate analysis results: paired

Paired STAY vs standard sections of four patients (preprocessing, Figures 3, 6 and 7).

## Layout

`<sample_id>/`: H&E image cropped to the capture area and shifted spot coordinates, from `01_calculate_crop_coords.py`.

`tls/<sample_id>.csv`: TLS score per spot (barcode, TLS_score), from `code/figures/tls_paired/tls_paired.py`.

Sample IDs: `Patient_<n>_STAY` / `Patient_<n>_standard` (Patient 1-4). See `resources/sample_list.csv`.

## Source

- https://doi.org/10.5281/zenodo.23558413
