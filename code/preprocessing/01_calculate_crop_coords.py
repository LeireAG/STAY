#!/usr/bin/env python3

# --- Aim
# Calculate crop coords to crop the H&E image to be just the area covered by Visium spots

# For some sections the H&E image is the image of the entire tissue
# For others the H&E image was cropped to fit the capture area

# This is not an issue with scanpy as it automatically crops the H&E image to the capture area when plotting, but appears to be an issue when working in R
# As apparently it's harder to crop the area plotted just to the Visium spots rather than plotting on the entire tissue

# --- Define environment
# env: hard-vs-soft

# --- Load packages
import scanpy as sc
import anndata as ad
import pandas as pd

import numpy as np
import os

import matplotlib.pyplot as plt

import argparse

# sc.logging.print_versions()
# sc.set_figure_params(facecolor="white", figsize=(8, 8))
# sc.settings.verbosity = 3

# Custom functions
from sophiesfunctions.auto_crop import auto_crop_scanpy  # import my custom functions that are packaged locally as sophiesfunctions
from sophiesfunctions.save_multi_image import save_multi_image


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & do basic QC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--res_dir", help="Path to results directory")
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--min_spots", help="Minimal spots per gene", type=int)

args = vars(parser.parse_args())

# --- Arguments
res_dir = args["res_dir"]
h5ad_dir = args["h5ad_dir"]
sample_id = args["sample_id"]
project_dir = args["project_dir"]

# --- Get some info from the sample sheet

# Extract the external_id
sample_dir = os.path.join(project_dir, 'config')
sample_sheet_file = 'sample_sheet.csv'

sample_sheet = pd.read_csv(os.path.join(sample_dir, sample_sheet_file))
external_id = sample_sheet[sample_sheet['section_id'] == sample_id]['external_id'].values[0]
tissue_id = sample_sheet[sample_sheet['section_id'] == sample_id]['tissue_id'].values[0]
slide_size = sample_sheet[sample_sheet['section_id'] == sample_id]['slide_size'].values[0]

joined_id = f'{sample_id}-{external_id}'

# --- Also get some tissue metadata

resource_dir = os.path.join(project_dir, 'resources')
metadata_file = 'tissue_metadata.csv'

tissue_metadata = pd.read_csv(os.path.join(resource_dir, metadata_file))

# Extract subtype
subtype = tissue_metadata[tissue_metadata['tissue_id'] == tissue_id]['clinical_subtype'].values[0]


# --- Load the adata object

use_data = 'log_norm'
adata = ad.read_h5ad(os.path.join(h5ad_dir, use_data, f'{sample_id}_{use_data}.h5ad'))

# Specify which layer to use
use_layer = 'log_norm'

# --- Calculate crop coords for optimal presentation
coords = auto_crop_scanpy(adata, border_size=0.1)


# --- Extract the hires image

img = adata.uns['spatial'][sample_id]['images']['hires']

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.imshow(img)
plt.suptitle(f'{sample_id} ({subtype}) - Hires image in object', y=1, fontweight='bold')


# --- Shift the spot coordinates

# Get the scalefactor for the hires image
scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

# Adjust coordinates for the hires image
coords_hires = tuple(int(x * scalef) for x in coords)

# Function to crop image, introducing white padding if crop is beyond borders

def crop_image_with_padding(img, coords):
    # Calculate the crop coordinates
    x_min = coords[0]
    x_max = coords[1]
    y_min = coords[2]
    y_max = coords[3]
    # Determine the number of channels (3 for RGB, 4 for RGBA)
    num_channels = img.shape[2] if len(img.shape) > 2 else 1
    # Determine the valid cropping area within the image
    x_min_valid = max(x_min, 0)
    x_max_valid = min(x_max, img.shape[1])
    y_min_valid = max(y_min, 0)
    y_max_valid = min(y_max, img.shape[0])
    # Determine the corresponding area in the crop image
    crop_x_min = x_min_valid - x_min
    crop_y_min = y_min_valid - y_min
    crop_x_max = crop_x_min + (x_max_valid - x_min_valid)
    crop_y_max = crop_y_min + (y_max_valid - y_min_valid)
    # Initialize the output image with white background based on data type
    if np.issubdtype(img.dtype, np.floating):
        white_value = 1.0  # For float images (0..1)
    else:
        white_value = 255  # For uint8 images (0..255)
    crop_height = y_max - y_min
    crop_width = x_max - x_min
    crop_img = np.ones((crop_height, crop_width, num_channels), dtype=img.dtype) * white_value
    # Copy the valid image area to the output image
    crop_img[crop_y_min:crop_y_max, crop_x_min:crop_x_max] = img[y_min_valid:y_max_valid, x_min_valid:x_max_valid]
    return crop_img, (x_min_valid, x_max_valid, y_min_valid, y_max_valid)

# Crop the H&E image
crop_img, adjusted_coords = crop_image_with_padding(img, coords_hires)

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.imshow(crop_img)
plt.suptitle(f'{sample_id} ({subtype}) - Hires image after cropping', y=1, fontweight='bold')


# --- Shift the spot coordinates after cropping the image

# Make a copy of anndata & add the crop in the hires slot
adata_crop = adata.copy()
adata_crop.uns['spatial'][sample_id]['images']['hires'] = crop_img

# Shift the spot coordinates so they're still in the same location even after cropping the image
spot_coords = adata.obsm['spatial']
shifted_coords = spot_coords - np.array([coords[0], coords[2]])

adata_crop.obsm['spatial'] = shifted_coords

# Also shift the original coords
# Identify the minimum x and y
x_min_shift = coords[0]
y_min_shift = coords[2]

# Shift the coordinates
transformed_coords = (
    coords[0] - x_min_shift,  # new x_min = 0
    coords[1] - x_min_shift,  # new x_max
    coords[2] - y_min_shift,  # new y_min = 0
    coords[3] - y_min_shift   # new y_max
)


# --- Plot the result after cropping
# Check spots are in the same location before after

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8*2, 6),)

sc.pl.spatial(adata, color="in_tissue", crop_coord=coords, ax=axs[0], title='Spots pre cropping (cropping with scanpy)', img_key='hires')

sc.pl.spatial(adata_crop, color="in_tissue", crop_coord=transformed_coords, ax=axs[1], title='Spots post cropping (no cropping with scanpy)', img_key='hires')

plt.suptitle(f'{sample_id} ({subtype}) - Spots on hires image before + after cropping', y=1, fontweight='bold')

# --- Also shift the original tissue_positions coordinates
# In case we want to make an object from scratch
# Load original tissue_positions file
spatial_dir = os.path.join(project_dir, 'data', 'processed', 'spaceranger_outs', joined_id, 'outs', 'spatial')
tissue_positions_file = 'tissue_positions.csv'

tissue_positions = pd.read_csv(os.path.join(spatial_dir, tissue_positions_file))

# Shift the coordinates
tissue_positions_shifted = tissue_positions.copy()

tissue_positions_shifted['pxl_row_in_fullres'] = tissue_positions['pxl_row_in_fullres'] - coords[2]  # Shift rows by coords[2]
tissue_positions_shifted['pxl_col_in_fullres'] = tissue_positions['pxl_col_in_fullres'] - coords[0]  # Shift columns by coords[0]


# --- Save the output

out_dir = os.path.join(project_dir, 'results/misc/crop_images/per_sample', sample_id)
os.makedirs(out_dir, exist_ok=True)

# -- Save the cropped image
filename = 'tissue_hires_image.png'
plt.imsave(os.path.join(out_dir, filename), crop_img)

# -- Save the shifted spot coordinates
# Convert the array to a DataFrame
shifted_coords_df = pd.DataFrame(shifted_coords, columns=['x', 'y'])

# Add obs_names as index
shifted_coords_df.index = adata.obs_names

filename = 'shifted_spot_coordinates.csv'
shifted_coords_df.to_csv(os.path.join(out_dir, filename))

# -- Save the shifted tissue_positions (as if it came from spaceranger)

filename = 'tissue_positions_shifted.csv'
tissue_positions_shifted.to_csv(os.path.join(out_dir, filename), index=False)

# -- Also save the original coordinates
original_coords_df = pd.DataFrame(spot_coords, columns=['x', 'y'])

# Add obs_names as index
original_coords_df.index = adata.obs_names

filename = 'original_spot_coordinates.csv'
original_coords_df.to_csv(os.path.join(out_dir, filename))

# -- Save the plots
filename = f'{sample_id}_HE_cropping_hires_plots.pdf'
save_multi_image(os.path.join(out_dir, filename))

# Close images
plt.close('all')

# --- Also crop the lowres image
# Originally I thought semla would probably use the hires image (default for scanpy)
# But from the tutorial it seems like it might be using the lowres only

# --- Extract the lowres image

img = adata.uns['spatial'][sample_id]['images']['lowres']

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.imshow(img)
plt.suptitle(f'{sample_id} ({subtype}) - Lowres image in object', y=1, fontweight='bold')

plt.savefig('test.pdf', bbox_inches='tight')


# --- Shift the spot coordinates

# Get the scalefactor for the hires image
scalef = adata.uns['spatial'][sample_id]['scalefactors']['tissue_lowres_scalef']

# Adjust coordinates for the lowres image
coords_lowres = tuple(int(x * scalef) for x in coords)

# Crop the H&E image
crop_img, adjusted_coords = crop_image_with_padding(img, coords_lowres)

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.imshow(crop_img)
plt.suptitle(f'{sample_id} ({subtype}) - Lowres image after cropping', y=1, fontweight='bold')


# --- Shift the spot coordinates after cropping the image

# Make a copy of anndata & add the crop in the lowres slot
adata_crop = adata.copy()
adata_crop.uns['spatial'][sample_id]['images']['lowres'] = crop_img

# Shift the spot coordinates so they're still in the same location even after cropping the image
spot_coords = adata.obsm['spatial']
shifted_coords = spot_coords - np.array([coords[0], coords[2]])

adata_crop.obsm['spatial'] = shifted_coords

# Also shift the original coords
# Identify the minimum x and y
x_min_shift = coords[0]
y_min_shift = coords[2]

# Shift the coordinates
transformed_coords = (
    coords[0] - x_min_shift,  # new x_min = 0
    coords[1] - x_min_shift,  # new x_max
    coords[2] - y_min_shift,  # new y_min = 0
    coords[3] - y_min_shift   # new y_max
)


# --- Plot the result after cropping
# Check spots are in the same location before after

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8*2, 6),)

sc.pl.spatial(adata, color="in_tissue", crop_coord=coords, ax=axs[0], title='Spots pre cropping (cropping with scanpy)', img_key='lowres')

sc.pl.spatial(adata_crop, color="in_tissue", crop_coord=transformed_coords, ax=axs[1], title='Spots post cropping (no cropping with scanpy)', img_key='lowres')

plt.suptitle(f'{sample_id} ({subtype}) - Spots on lowres image before + after cropping', y=1, fontweight='bold')


# -- Save the cropped image
filename = 'tissue_lowres_image.png'
plt.imsave(os.path.join(out_dir, filename), crop_img)

# -- Save the plots
filename = f'{sample_id}_HE_cropping_lowres_plots.pdf'
save_multi_image(os.path.join(out_dir, filename))


# --- Finish
print('Script successfully completed')