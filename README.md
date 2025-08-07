# STAY
# breast_cancer_visium_analysis
This repository contains code for Visium data preprocessing, downstream analysis, and figure generation for the study titled "Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings". We refer to this protocol as Spatial Transcriptomics with Advanced tissue preparation for preserved data Yield (STAY). 
# Data analysed
In-house FFPE breast cancer tissue sections. 
4 tissue sections analysed in 6.5x6.5 mm capture areas with Visium CytAssist FFPE chemistry
4 tissue sections analysed in 11x11 mm capture areas with Visium CytAssist FFPE chemistry
# Overview of directories 
# code/
Contains the scrips to generate the analysis and figures in this study. Data paths in the scripts may need to be adjusted. By default, they assume a local directory structure matching this repository. Data needs to be downloaded from KTH Data Repository at insert_DOI. 

 # preprocessing/
- 01_calculate_crop_coords.py - Crops H&E images to the 6.5x6.5 mm or 11x11 mm area and adjusts spot coordinates.  
- 02_make_semla_object.R - The cropped images and adjusted spot coordinates are used to create a Semla object all_samples_cropped_img.rds.
- 03_loading_reviewed_path_annotation copy.Rmd - Loads pathology annotations, generates spot annotation plots in Figure 1 C1-4, F1-4, violin plots in Figure 3 J, K and stats in Supplementary Data 1. 
- 04_QC_normalization_scaling.Rmd - Steps of quality control, scaling and normalization of the data. 
- qupath_to_loupe.py - Produces spot based pathology annotation in .csv format to be used in the 03_loading_reviewed_path_annotation copy.Rmd script

 # figs/ 
- renv - Copy of Rstudio local environment where we performed the analysis. Contains the versions of R packages. 
- renv.lock - Lists exact package versions in the project-specific environment 
- figure_4/
- figure_1.Rmd
- figure_3.Rmd
- figure_3_on_reviewed_annotations_copy.Rmd
- supplmentary_figs_4_to_12.Rmd

# Data access 
Full resolution hematoxylin-eosin images and Seurat objects will be uploaded to the KTH Data Repository. 

# Confidentiality 
This repository is shared in the context of peer review. Please do not use or redistribute the STAY protocol or code prior to publication.

# Contact 
First author: Leire Alonso Galicia leire.alonso@scilifelab.se 
