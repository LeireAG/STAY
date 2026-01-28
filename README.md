# STAY
# breast_cancer_visium_analysis
This repository contains code for Visium data preprocessing, downstream analysis, and figure generation for the study titled "Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings". We refer to this protocol as Spatial Transcriptomics with Advanced tissue preparation for preserved data Yield (STAY). 
# Data analysed
In-house FFPE breast cancer tissue sections. 
4 tissue sections analysed in 6.5x6.5 mm capture areas with Visium CytAssist FFPE chemistry
4 tissue sections analysed in 11x11 mm capture areas with Visium CytAssist FFPE chemistry
# Overview of directories 
# code/
Contains the scripts to generate the analysis and figures in this study. Data paths in the scripts may need to be adjusted. By default, they assume a local directory structure matching this repository. Data needs to be downloaded from KTH Data Repository at insert_DOI. Please check Data availability section. 

 # preprocessing/
- 01_calculate_crop_coords.py - Crops H&E images to the 6.5x6.5 mm or 11x11 mm area and adjusts spot coordinates.  
- 02_make_semla_object.R - The cropped images in crop_images directory and adjusted spot coordinates are used to create a seurat object all_samples_cropped_img.rds.
- 03_loading_reviewed_path_annotation copy.Rmd - Loads pathology annotations in the previous seurat object vis_mod_anno.rds, generates spot annotation plots in Figure 1 C1-4, F1-4, violin plots in Figure 3 J, K and stats in Supplementary Data 1. Contains the code to generate spot path annotation plots in Figure 1 and Suppl Figs. 
- 04_QC_normalization_scaling.Rmd - Steps of quality control, scaling and normalization of the data. Generates a normalized Seurat object vis_mod_norm used in figure_3_on_reviewed_annotations_copy.Rmd. 
Contains the code to generate the UMAP in Figure 3F.
- qupath_to_loupe.py - Produces spot based pathology annotation in .csv format in reworked_annotations directory to be used in the 03_loading_reviewed_path_annotation copy.Rmd script

 # figs/ 
- renv - Copy of Rstudio local environment where we performed the analysis. Contains the versions of R packages. 
- renv.lock - Lists exact package versions in the project-specific environment 
- figure_4/ - generates a normalized with corrected number of spots Seurat object including nmf and deconvolution vis_mod_nmf_decon.rds
* figure_4_copy.Rmd contains the code to run Non-Negative_matrix Factorization and saving that Seurat object. Generates fig. 4 A, B, D1-I1 and D2-I2. Used for ploting factors spatially for supplementary figs 3-12. 
* semla_deconvolution.rmd contains the code to run the NNLS deconvolution method included in the Semla package on the previous object containing NMF factors
* cell_proportion_barplot_figure_4K.Rmd contains the code to generate stacked barplots of cell-type proportions for Patients 3 and 4.  
- figure_3_on_reviewed_annotations_copy.Rmd generates spatial plots of normalized unique genes per spot (figure 3 A-H), spot number barplots in Figure 1 G1-4 and violin plots in Figure 3 I. 
- scatter_plots_figure_3E.Rmd contains the code to generate correlation plots for all patients

# Data availability
Full resolution hematoxylin-eosin images and Seurat objects will be uploaded to the KTH Data Repository insert DOI. 

Although intermediate Seurat objects were generated during analysis (one per major processing step), only the final object is required to run the scripts here. The final seurat object contains all relevant information including pathology annotations, normalized data, NMF results, and deconvolution. 
The object can be downloaded from the KTH Data Repository at: insert_DOI_here
Note: In the scripts, variable names for the Seurat object may differ between steps, reflecting the sequential workflow used during development.

# Confidentiality 
This repository is shared in the context of peer review. Please do not use or redistribute the STAY protocol or code prior to publication.

# Contact 
First author: Leire Alonso Galicia leire.alonso@scilifelab.se 
