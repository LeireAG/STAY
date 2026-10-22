Processed data (https://doi.org/10.5281/zenodo.23558413) for the STAY study: "Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings".

Code: <GITHUB_URL>

## Files

Each zip holds one data type for one cohort. Inside, paths start at the repository root (`data/<type>/<cohort>/...`).

| Zip | Contents |
|---|---|
| `spaceranger_paired.zip` | Space Ranger output, paired cohort |
| `spaceranger_standard.zip` | Space Ranger output, standard cohort |
| `histopathology_annotation_paired.zip` | QuPath annotations and spot-level annotations, paired cohort |
| `histopathology_annotation_standard.zip` | QuPath annotations and spot-level epithelial annotations, standard cohort |
| `analysis_paired.zip` | Cropped H&E images and shifted spot coordinates, paired cohort |
| `analysis_STAY.zip` | Per-spot detachment scores, STAY cohort |
| `analysis_standard.zip` | Per-spot detachment scores, standard cohort |
| `r_objects_paired.zip` | `paired_object.rds`: final Seurat/semla object, paired cohort |
| `r_objects_lymph_node.zip` | `lymph_node_object.rds`: final Seurat/semla object, lymph node cohort |

Raw data (`data/raw_data/`) is deposited separately: <RAW_DATA_ACCESSION>

STAY cohort and lymph node data (raw data, Space Ranger output, annotations): https://doi.org/10.6019/S-BIAD3080 and https://doi.org/10.5281/zenodo.18917100

## How to use

```bash
git clone <GITHUB_URL> STAY
cd STAY
# download the zips into this folder, then:
for f in *.zip; do unzip -n "$f"; done
```

Every file lands in its correct `data/` folder. Then follow the repository README to set up the environments and run the analysis.
