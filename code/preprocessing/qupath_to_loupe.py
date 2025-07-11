#!/usr/bin/env python3
#
#
#
# The MIT License (MIT)

# Copyright (c) 2018 10x Genomics, Inc.

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
# THE SOFTWARE.

__version_info__ = ('2024','02','12')
__version__ = '-'.join(__version_info__)

import argparse
import numpy as np
import polars as pl
import json
import os
from shapely.geometry import shape, GeometryCollection, Point, box
import pandas as pd
from collections import defaultdict

def _load_geojson(qupath_geojson_file):
    """Loads the geometries (should be Polygon and MultiPolygon) and their associated names from QuPath GeoJSON
    This could also be replaced by using GeoPandas to do the reading.
    """
    print(f"loading QuPath GeoJSON file from {qupath_geojson_file}")
    with open(qupath_geojson_file) as f:
        features = json.load(f)["features"]
    # Deal with annotations that do not have a label (and therefore no classification)
    # It's quite common to get a file where at least one annotation doesn't have a label due to oversight
    # Previously this would break the script
    names = [
        feature["properties"]["classification"]["name"]
        if "classification" in feature["properties"] and "name" in feature["properties"]["classification"]
        else 'no_label'
        for feature in features
    ]
    geoms = [shape(feature["geometry"]).buffer(0) for feature in features]
    colors = [
        feature["properties"]["classification"]["color"]
        if "classification" in feature["properties"] and "color" in feature["properties"]["classification"]
        else [0, 0, 0]
        for feature in features
    ]
    return names, geoms, colors


def _load_barcode_positions(space_ranger_outs,):
    """Takes the location of a Space Ranger outs folder and loads tissue positions.  Original Visium only
    for this version, but will be updated for Visium HD shortly.
    """
    pathname=os.path.join(space_ranger_outs, "spatial", "tissue_positions.csv")
    coords = pl.read_csv(pathname)
    print(f"loading Space Ranger tissue positions from {pathname}")
    return coords


def _load_spot_diameter(space_ranger_outs,):
    """Takes the location of a Space Ranger outs folder and loads the spot diameter in the fullres image.
    """
    pathname=os.path.join(space_ranger_outs, "spatial", "scalefactors_json.json")
    print(f"loading Space Ranger spot size from {pathname}")
    with open(pathname, 'r') as json_file:
                data = json.load(json_file)
    spot_diameter = data['spot_diameter_fullres']
    return spot_diameter


def _optional_prefix(prefix, name):
    if prefix is not None:
        return os.path.join(prefix, name)
    else:
        return name

def convert_qupath_to_loupe(qupath_geojson_file, space_ranger_outs, output_prefix, consider_spot_diameter,):
    names, geoms, colors = _load_geojson(qupath_geojson_file)
    coords = _load_barcode_positions(space_ranger_outs,)
    spot_diameter = _load_spot_diameter(space_ranger_outs,)
    radius = spot_diameter / 2

    output = coords[("barcode",)]
    for geom, name in zip(geoms, names):
        print(f"computing on {name}")
        if consider_spot_diameter:
            print('Considering the whole spot size (shape: circle)')
        s = []
        if consider_spot_diameter:
            for _, _, _, _, r, c in coords.rows():  # same order as "coords" and "output"
                circle = Point(c, r).buffer(radius)
                s.append(geom.intersects(circle))  # Check if geom intersects with the circle
        else:
            for _, _, _, _, r, c in coords.rows():  # same order as "coords" and "output"
                s.append(geom.contains(Point(c, r)))
        # Deal with multiple separate annotations of the same category (same name)    
        if name in output.columns:
            existing_series = output[name].to_list()
            combined_series = [existing or new for existing, new in zip(existing_series, s)]
            output = output.with_columns(pl.Series(name, combined_series))
        else:
            output = output.with_columns(pl.Series(name, s))

    print("writing output files...")
    # write output with one column per annotation and true/false for membership
    output.write_csv(_optional_prefix(output_prefix, "all_annotations.csv"))

    # write Loupe(tm) compatible output which CHOOSES ONE annotation per barcode
    output.melt(id_vars="barcode", variable_name="annotation").filter(
        pl.col("value") == True
    )[["barcode", "annotation"]].unique(subset="barcode").write_csv(
        _optional_prefix(output_prefix, "loupe-one-per-barcode.csv")
    )

    # write Loupe(tm) compatible output with all combinations of annotations
    # per barcode
    new_output = output[["barcode"]]
    for name in names:
        new_output = new_output.with_columns(
            output[name].replace({True: name, False: pl.Null}).alias(name)
        )
    new_output = new_output.with_columns(
    pl.concat_str(names, separator='+', ignore_nulls=True).alias("annotation")
    ).filter(pl.col("annotation") != "")
    # Function to remove duplicates in concatenated string
    def remove_duplicates(annotation):
        return '+'.join(sorted(set(annotation.split('+'))))
    # Remove duplicates from annotation column using map_elements
    new_output = new_output.with_columns(
        pl.col("annotation").map_elements(remove_duplicates, return_dtype=pl.Utf8).alias("annotation")
    )
    new_output[["barcode", "annotation"]].write_csv(
        _optional_prefix(output_prefix, "loupe-all-per-barcode.csv")
    )

    # -- Identify spots with overlapping labels and assign dominant label only (or if multiple dominant labels: mixed)
    # Select only the boolean columns (excluding 'barcode')
    bool_cols = output.columns[1:]
    # Filter rows where at least 2 columns are True
    filtered_output = output.filter(output.select(bool_cols).sum(axis=1) >= 2)
    results = []

    # Loop through selected rows
    for row in filtered_output.iter_rows(named=True):
        barcode = row["barcode"]
        true_cols = [col for col in bool_cols if row[col]]  # Extract columns with True
        # Get coordinates
        coord_row = coords.filter(coords["barcode"] == barcode).row(0, named=True)
        pxl_col, pxl_row = coord_row["pxl_col_in_fullres"], coord_row["pxl_row_in_fullres"]
        # Create a circle
        circle = Point(pxl_col, pxl_row).buffer(radius)  # Creates a circle with the specified radius
        spot = circle
        spot_area = spot.area
        # Step 1: Quick filtering using bounding box intersection
        filtered_geoms = [
            (geom, name) for geom, name in zip(geoms, names)
            if geom.bounds[2] >= spot.bounds[0] and  # geom max_x >= square min_x
            geom.bounds[0] <= spot.bounds[2] and  # geom min_x <= square max_x
            geom.bounds[3] >= spot.bounds[1] and  # geom max_y >= square min_y
            geom.bounds[1] <= spot.bounds[3]      # geom min_y <= square max_y
        ]
        # Step 2: Exact intersection check
        overlapping_geoms = [(geom, name) for geom, name in filtered_geoms if geom.intersects(spot)]
        # Compute overlap with each True category
        for category in true_cols:
            for geom, name in overlapping_geoms:
                if name == category:  # Match annotation name with category
                    overlap_area = spot.intersection(geom).area
                    percent_overlap = (overlap_area / spot_area) * 100
                    results.append({
                        "barcode": barcode,
                        "category": category,
                        "percent_overlap": percent_overlap
                    })
    # Group by barcode and store categories with their overlap percentages
    barcode_to_overlaps = defaultdict(list)

    for entry in results:
        barcode = entry['barcode']
        barcode_to_overlaps[barcode].append((entry['category'], entry['percent_overlap']))

    # Determine the highest overlap label for each barcode
    final_labels = {}

    for barcode, overlaps in barcode_to_overlaps.items():
        # Sort by overlap percentage in descending order
        overlaps.sort(key=lambda x: x[1], reverse=True)
        # Find the maximum overlap percentage
        max_overlap = overlaps[0][1]
        # Check for ties: find all categories that have the max overlap
        top_categories = [cat for cat, perc in overlaps if perc == max_overlap]
        # Make sure to remove any duplicate labels
        top_categories = list(set(top_categories))
        if len(top_categories) == 1:
            final_labels[barcode] = top_categories[0]  # Assign the single highest category
        else:
            # Sort to ensure uniform order
            top_categories = sorted(top_categories)
            final_labels[barcode] = f"mixed_{'_'.join(top_categories)}"  # Assign mixed category

    # Convert the result into a DataFrame
    final_labels_df = pd.DataFrame(list(final_labels.items()), columns=["barcode", "label"])
    # Add to original df
    df = pd.DataFrame(new_output[["barcode", "annotation"]], columns=['barcode', 'annotation'])
    # Merge the final labels DataFrame with your original DataFrame on 'barcode'
    df_updated = df.merge(final_labels_df, on='barcode', how='left')
    # Overwrite the 'annotation' column with the 'label' from final_labels_df
    df_updated['annotation'] = np.where(df_updated['label'].notna(),
                                        df_updated['label'],
                                        df_updated['annotation'])
    # Drop the 'label' column as it's no longer needed
    df_updated = df_updated.drop(columns=['label'])
    # Save to csv
    df_updated.to_csv(
        _optional_prefix(output_prefix, "loupe-dominant-per-barcode.csv"),
        index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="qupath_to_loupe",
        description="Convert QuPath annotations and Space Ranger output to Loupe CSVs for input",
    )
    parser.add_argument('-v', '--version', action='version', version="%(prog)s ("+__version__+")")
    # parser.add_argument("qupath_geojson_file", help="GeoJSON-format annotations on image used as input to Space Ranger")
    # parser.add_argument("space_ranger_outs", help="pathname of the spaceranger outs directory")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        type=str,
        help="Prefix pathname to prepend to output filenames",
    )

    parser.add_argument("--sample_id", help="Sample ID")
    parser.add_argument("--res_dir", help="Results directory")
    parser.add_argument("--input_dir", help="General path to input")
    parser.add_argument("--spaceranger_dir", help="Path to spaceranger outs (if not provided, will be determined from input_dir)")
    parser.add_argument("--qupath_dir", help="Path to qupath annotations (if not provided, will be determined from input_dir)")
    parser.add_argument("--consider_spot_diameter", help="Whether to consider the full size of the spot for label transfer, otherwise considers only the center coordinate (default: True)", 
                        action='store_false', default=True)

    args = parser.parse_args()

    sample_id = args.sample_id
    input_dir = args.input_dir
    res_dir = args.res_dir
    spaceranger_dir = args.spaceranger_dir
    qupath_dir = args.qupath_dir
    consider_spot_diameter = args.consider_spot_diameter

    # Path to qupath file
    if qupath_dir is None:
        qupath_dir = os.path.join(input_dir, sample_id, 'qupath_annotations')
    else:
        # If not absolute, join with input_dir
        if not os.path.isabs(qupath_dir):
            qupath_dir = os.path.join(input_dir, qupath_dir)
        # Only append sample_id if not already present at the end
        if not qupath_dir.rstrip(os.sep).endswith(sample_id):
            qupath_dir = os.path.join(qupath_dir, sample_id)
    
    # List all files in the directory
    files = os.listdir(qupath_dir)
    qupath_files = [file for file in files if file.endswith('.geojson')]

    print(f'Qupath annotations: {qupath_files}')

    if spaceranger_dir is None:
        spaceranger_dir = os.path.join(input_dir, sample_id, 'spaceranger_outs', 'outs')
    else:
        # If not absolute, join with input_dir
        if not os.path.isabs(spaceranger_dir):
            spaceranger_dir = os.path.join(input_dir, spaceranger_dir)
        # Only append sample_id/outs if not already present at the end
        expected_suffix = os.path.join(sample_id, 'outs')
        if not spaceranger_dir.rstrip(os.sep).endswith(expected_suffix):
            spaceranger_dir = os.path.join(spaceranger_dir, sample_id, 'outs')

    print(f'Using spaceranger outs: {spaceranger_dir}')

    # Results dir
    res_dir = os.path.join(res_dir, sample_id)
    # Make directory
    os.makedirs(res_dir, exist_ok=True)

    for qupath_file in qupath_files:
        # Make a results directory per file
        # Split filename on last occurrence of dot
        filename = qupath_file.rsplit('.', 1)[0]
        res_dir_file = os.path.join(res_dir, filename)
        os.makedirs(res_dir_file, exist_ok=True)

        # Convert annotations
        convert_qupath_to_loupe(
            os.path.join(qupath_dir, qupath_file), spaceranger_dir, res_dir_file, consider_spot_diameter,
        )