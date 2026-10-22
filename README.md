# STAY

Code for the study **"Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings"**. STAY stands for Spatial Transcriptomics with Advanced tissue preparation for preserved data Yield: Visium CytAssist FFPE on in-house breast cancer sections.

## Repository structure

```
code/
  preprocessing/     # paired cohort: cropping, object creation, annotation, QC
  figures/           # one script or folder per figure
data/                # not tracked by git: download the data into the matching folders
  raw_data/                    # FASTQs, H&E/CytAssist images, Loupe alignments
  spaceranger/                 # Space Ranger outs/ per sample
  histopathology_annotation/   # QuPath GeoJSON + spot-level CSVs
  analysis/                    # intermediate results (detachment scores, cropped images)
  r_objects/                   # Seurat/semla .rds objects
environment/         # conda env (Python) and renv (R)
resources/           # sample list and Figure 2 config
results/             # figure outputs, created by the scripts
```

Each data folder is split into the four cohorts:

| Cohort | Description | Sample IDs | Figures |
|---|---|---|---|
| `paired` | Paired STAY vs standard sections of 4 patients (4 in 6.5 x 6.5 mm, 4 in 11 x 11 mm capture areas) | `Patient_<n>_STAY`, `Patient_<n>_standard` | 3, 6, 7 |
| `lymph_node` | Lymph node sections of 4 patients (Patient_1' ... Patient_4') | Visium slide ID, e.g. `V53M13-062_A1` | 8 |
| `STAY` | STAY sections of the detachment analysis | `P01` ... `P32` (`P29_1`, `P29_2`) | 2 |
| `standard` | Standard-protocol sections of the detachment analysis | `S01` ... `S20` | 2 |

Data goes in `data/<type>/<cohort>/<sample_id>/`. Each cohort folder has a `README.md` with its layout and download source. `resources/sample_list.csv` lists every sample with its cohort, patient ID and (STAY) subtype. The scripts read their samples from it.

## Data availability

- **Processed data** (Space Ranger output, annotations, analysis results, Seurat objects): https://doi.org/10.5281/zenodo.23558413
- **Raw data** (`data/raw_data/`: FASTQs, images, Loupe alignments): <RAW_DATA_ACCESSION>
- **STAY cohort and lymph node:** https://doi.org/10.6019/S-BIAD3080 and https://doi.org/10.5281/zenodo.18917100

The final R objects are provided, so the figure scripts can run without redoing the preprocessing:

- `data/r_objects/paired/paired_object.rds`: pathology annotations, normalized data, NMF and deconvolution.
- `data/r_objects/lymph_node/lymph_node_object.rds`: QC-filtered lymph node object with Harmony and NMF.

## Setup

```bash
conda env create -f environment/detachment-score.yaml    # Python (Figure 2)
Rscript -e 'renv::restore(project = "environment")'      # R packages
```

In R, load the environment with `renv::load("environment")` before running a notebook. All paths are relative to the repository root. Run the Python scripts and `02_make_semla_object.R` from the root, and knit the notebooks from their own folder (the RStudio default).

## Running the analysis

**Preprocessing (paired)**, in order, in `code/preprocessing/`:

1. `01_calculate_crop_coords.py`: crops H&E images to the capture area and shifts the spot coordinates. Its output is provided in `data/analysis/paired/`.
2. `02_make_semla_object.R`: builds the semla object from the cropped images, with sections in `sample_list.csv` order.
3. `qupath_to_loupe.py`: converts QuPath annotations to spot-level CSVs (`--cohort`, default `paired`). Its output is provided in `data/histopathology_annotation/`.
4. `03_loading_reviewed_path_annotation copy.Rmd`: adds the pathology annotations (Figure 3 C/F, Figure 6 J/K, Supplementary Data 1).
5. `04_QC_normalization_scaling.Rmd`: QC, normalization and scaling (Figure 6F). Run it on `paired_object.rds` to create `vis_mod_norm.rds`, which the Figure 6 notebook needs.

**Figure 2 (STAY vs standard)**: settings are in `resources/detachment.yaml`.

```bash
conda activate detachment-score
python code/figures/figure_2/run_detachment_pipeline.py   # per-spot detachment scores
python code/figures/figure_2/plot_detachment_cohort.py    # cohort plots and statistics
python code/figures/figure_2/plot_detachment_spatial.py   # per-sample spatial panels
```

**Other figures**, in `code/figures/`:

| Script | Figure |
|---|---|
| `figure_6_on_reviewed_annotations_copy.Rmd` | 6 A-H, 6C, 3 G |
| `scatter_plots_figure_6E.Rmd` | 6E |
| `figure_7/figure_7_copy.Rmd` | 7 A, B, D-I; Supplementary Figures 4-11 (NMF) |
| `figure_7/semla_deconvolution.rmd` | deconvolution used in Figure 7 |
| `figure_7/cell_proportion_barplot_figure_7K.Rmd` | 7K |
| `lymph_node_analysis.qmd` | 8 |

`semla_deconvolution.rmd` also needs an external single-cell reference (`data/r_objects/paired/miniatlas.rds`), which is not included. Its results are already in `paired_object.rds`. `lymph_node_analysis.qmd` rebuilds `lymph_node_object.rds` from the Space Ranger output.

Outputs go to `results/<figure>/`, and to `results/preprocessing/` for the preprocessing plots.

## Confidentiality

This repository is shared for peer review. Please do not use or redistribute the STAY protocol or code before publication.

## Contact

First author: Leire Alonso Galicia [leire.alonso@scilifelab.se](mailto:leire.alonso@scilifelab.se)
