# Aim: Replace default images in semla object with cropped images
# Images cropped to Visium capture area

# --- Activate the renv
# To run interactively - activate r-conda environment first, then radian in terminal
project_dir = getwd()  # run from the repo root
renv_dir = file.path(project_dir, "environment")
renv::load(renv_dir)

# --- Load packages

library(semla)


# --- Make a new object with the cropped images

# Have cropped the post spaceranger images to be just the visium capture area + a boundary
# Also translated the spot coordinates to fit these images
# Make a new object from scratch containing the cropped images

# Paired sections, in the order of resources/sample_list.csv (= section numbers 1-8 used in the later scripts)
sample_list <- read.csv(file.path(project_dir, "resources", "sample_list.csv"))
paired <- sample_list[sample_list$cohort == "paired", ]
sample_id_list <- paired$sample_id

spaceranger_dirs <- file.path(project_dir, "data/spaceranger/paired", sample_id_list)
cropped_dirs <- file.path(project_dir, "data/analysis/paired", sample_id_list)

# Get the required info
samples <- file.path(spaceranger_dirs, "filtered_feature_bc_matrix.h5")
imgs <- file.path(cropped_dirs, "tissue_hires_image.png")
spotfiles <- file.path(cropped_dirs, "tissue_positions_shifted.csv")
json <- file.path(spaceranger_dirs, "spatial", "scalefactors_json.json")

# Put everything into infotable; extra columns become metadata (sample_id = Patient_<n>_<method>)
infoTable <- tibble(samples, imgs, spotfiles, json, # Add required columns
                    sample_id = sample_id_list,
                    patient = paired$patient_id,
                    method = sub(".*_", "", sample_id_list))
infoTable$patient_method <- paste(infoTable$patient, infoTable$method, sep = "_")

# Make seurat object
se <- ReadVisiumData(infoTable)
se

# Get the spatial data so can plot some stuff
spatial_data <- GetStaffli(se)

# Load images
se <- LoadImages(se)


# --- Make some plots to check

out_dir_plots <- file.path(project_dir, "results/preprocessing/crop_images_check")
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(out_dir_plots, "images_in_object.pdf"), width = , height = 5)
ImagePlot(se)
dev.off()

pdf(file.path(out_dir_plots, "spatial_plot_with_HE.pdf"), width = 20, height = 10)
MapFeatures(se, features = "nFeature_Spatial", slot="counts", override_plot_dims = FALSE, image_use = "raw", pt_stroke = 0, ncol=4)
dev.off()

pdf(file.path(out_dir_plots, "spatial_plot_without_HE.pdf"), width = 20, height = 10)
MapFeatures(se, features = "nFeature_Spatial", slot="counts", override_plot_dims = FALSE, ncol=4)
dev.off()

# --- Save the object

out_dir <- file.path(project_dir, "data/r_objects/paired")
filename <- "all_samples_cropped_img.rds"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(se, file.path(out_dir, filename))