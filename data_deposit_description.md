Processed data for the STAY study: "Adapting spatial transcriptomics to FFPE tissue sections prepared with routine histology workflows in clinical settings".

Code: <GITHUB_URL>

## Files

Each zip holds one data type for one cohort (`<type>_<cohort>.zip`). Inside, paths start at the repository root (`data/<type>/<cohort>/<sample_id>/...`).

| Zip | Contents |
|---|---|
| `spaceranger_<cohort>.zip` | Space Ranger output per sample |
| `histopathology_annotation_<cohort>.zip` | QuPath annotations (GeoJSON) and spot-level annotations (CSV) |
| `analysis_<cohort>.zip` | Intermediate results (detachment scores, cropped images) |
| `r_objects_<cohort>.zip` | Seurat/semla objects (.rds) |

Cohorts: `paired`, `lymph_node`, `STAY`, `standard`.

Raw data (`data/raw_data/`) is deposited separately: <RAW_DATA_ACCESSION>

## How to use

```bash
git clone <GITHUB_URL> STAY
cd STAY
# download the zips into this folder, then:
for f in *.zip; do unzip -n "$f"; done
```

Every file lands in its correct `data/` folder. Then follow the repository README to set up the environments and run the analysis.
