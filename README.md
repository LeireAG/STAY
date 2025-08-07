# STAY
# breast_cancer_visium_analysis
This repository contains code for Visium data processing, downstream analysis, and figure generation for the study titled "Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings". We refer to this protocol as Spatial Transcriptomics with Advanced tissue preparation for preserved data Yield (STAY). 
# Data analysed
In-house FFPE breast cancer tissue sections. 
4 tissue sections analysed in 6.5x6.5 mm capture areas with Visium CytAssist FFPE chemistry
4 tissue sections analysed in 11x11 mm capture areas with Visium CytAssist FFPE chemistry
# Overview of directories 
# code/
figs/ 
- renv
- renv.lock
- figure_4/
- figure_1.Rmd
- figure_3.Rmd
- figure_3_on_reviewed_annotations_copy.Rmd
- supplmentary_figs_4_to_12.Rmd
preprocessing/
- 01_calculate_crop_coords.py
- 02_make_semla_object.R
- 03_loading_reviewed_path_annotation copy.Rmd
- 04_QC_normalization_scaling.Rmd
- qupath_to_loupe.py

