# Aim: Replace default images in semla object with cropped images
# Images cropped to Visium capture area

# --- Activate the renv
# To run interactively - activate r-conda environment first, then radian in terminal
project_dir = getwd()
renv_dir = file.path(project_dir, "environment/semla")
renv::load(renv_dir)

# --- Load packages

library(semla)


# --- Make a new object with the cropped images

# Have cropped the post spaceranger images to be just the visium capture area + a boundary
# Also translated the spot coordinates to fit these images
# Make a new object from scratch containing the cropped images

data_root_directory <- file.path(project_dir, "data/processed/spaceranger_outs", "*")
cropped_dir <- file.path(project_dir, "results/misc/crop_images/per_sample", "*")

# List directories at the top level only
all_dirs <- list.dirs(path = file.path(project_dir, "results/misc/crop_images/per_sample"), full.names = TRUE, recursive = FALSE)

# Remove the parent directory (if you want to exclude it)
directories <- all_dirs[all_dirs != file.path(project_dir, "results/misc/crop_images/per_sample")]

# Extract the names of the directories (correspond to the sample ids)
sample_id_list <- basename(directories)

# Get the required info
samples <- Sys.glob(paths = file.path(data_root_directory, "outs",
                                      "filtered_feature_bc_matrix.h5"))

imgs <- Sys.glob(paths = file.path(cropped_dir, 
                                   "tissue_hires_image.png"))

spotfiles <- Sys.glob(paths = file.path(cropped_dir, 
                                        "tissue_positions_shifted.csv"))

json <- Sys.glob(paths = file.path(data_root_directory, "outs",
                                   "spatial", "scalefactors_json.json"))

# Put everything into infotable
infoTable <- tibble(samples, imgs, spotfiles, json, # Add required columns
                    sample_id = sample_id_list) # Add additional column

# Make seurat object
se <- ReadVisiumData(infoTable)
se

# Get the spatial data so can plot some stuff
spatial_data <- GetStaffli(se)

# Load images
se <- LoadImages(se)


# --- Make some plots to check

out_dir_plots <- file.path(project_dir, "results/misc/crop_images/plots_to_check")
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

out_dir <- file.path(project_dir, "data/processed/semla_objects")
filename <- "all_samples_cropped_img.rds"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(se, file.path(out_dir, filename))