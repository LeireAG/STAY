# STAY
# breast_cancer_visium_analysis
This repository contains code for Visium data preprocessing, downstream analysis, and figure generation for the study titled "Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings". We refer to this protocol as Spatial Transcriptomics with Advanced tissue preparation for preserved data Yield (STAY). 
# Data analysed
In-house FFPE breast cancer tissue sections. 
4 tissue sections analysed in 6.5x6.5 mm capture areas with Visium CytAssist FFPE chemistry
4 tissue sections analysed in 11x11 mm capture areas with Visium CytAssist FFPE chemistry
# Overview of directories 
# code/
Contains the scrips to generate the analysis and figures in this study. Data paths in the scripts may need to be adjusted. By default, they assume a local directory structure matching this repository. 
 # figs/ 
- renv - requirements of R packages 
- renv.lock
- figure_4/
- figure_1.Rmd
- figure_3.Rmd
- figure_3_on_reviewed_annotations_copy.Rmd
- supplmentary_figs_4_to_12.Rmd

 # preprocessing/
- 01_calculate_crop_coords.py
- 02_make_semla_object.R
- 03_loading_reviewed_path_annotation copy.Rmd
- 04_QC_normalization_scaling.Rmd
- qupath_to_loupe.py

# Data access 
Full resolution hematoxylin-eosin images and Seurat objects will be uploaded to the KTH Data Repository. 

# Confidentiality 
This repository is shared in the context of peer review. Please do not use or redistribute the STAY protocol or code prior to publication.

# Contact 
First author: Leire Alonso Galicia leire.alonso@scilifelab.se 
