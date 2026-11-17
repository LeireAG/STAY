# Methods: Figure 2 (STAY vs standard detachment)

## Detachment annotation

Tissue detachment was annotated manually on the CytAssist images of all STAY (n = 33) and standard (n = 20) sections. For each section, the CytAssist image from the Space Ranger output (**[Space Ranger version]**) was displayed together with the high-resolution H&E image. The H&E was registered onto the CytAssist image by an affine transformation, estimated with RANSAC from the positions of the Visium spots in both images: the H&E spot coordinates from Space Ranger, and the CytAssist spot coordinates from the manual Loupe Browser alignment. This registration was used only to display the two images together. The H&E could be overlaid on the CytAssist image with adjustable opacity, so tissue present before processing (H&E) could be compared directly with the tissue remaining on the slide (CytAssist).

Annotation was done in a browser-based viewer built with Panel and HoloViews (Bokeh backend), served from the analysis server (`annotate_detachment.py`). Areas were classified as detached when **[criteria, e.g. tissue visible in the H&E was absent from, or displaced in, the CytAssist image]** and were outlined as polygons. To avoid bias, sections were presented in random order and identified only by their position in the sequence; sample identifiers and cohort were hidden throughout. All sections were annotated by **[one annotator / initials]**.

## Detachment quantification

Each polygon was rasterised onto the CytAssist image (`run_detachment_pipeline.py`). A spot was classified as detached if its centre, in CytAssist coordinates from the Loupe alignment, fell inside an annotated polygon. Because Space Ranger's tissue detection can exclude detached tissue, tissue spots were defined as spots called as tissue by Space Ranger, plus all detached spots. All spots were taken from the unfiltered feature-barcode matrix.

The detachment rate of a section was the percentage of its tissue spots classified as detached. The detached area was the number of annotated pixels multiplied by the pixel area. The CytAssist pixel size was derived from the 100 µm centre-to-centre distance between neighbouring Visium spots, using the median distance between neighbouring spots of the same row.

Detachment rates were also computed separately for epithelial and non-epithelial spots. Epithelial spots were taken from the pathologist's annotations **[annotated in QuPath, converted to spot level]**.

## RNA quality metrics

For each section, the median number of detected genes, the median UMI count, and the median percentage of reads from mitochondrial genes (genes starting with "MT-") were computed across its tissue spots (`plot_detachment_cohort.py`).

## Statistics

Sections were the unit of analysis. The two sections from one patient (P29) were stained independently and were treated as separate sections. STAY and standard sections were compared with two-sided Mann–Whitney U tests.

As the effect size, the Hodges–Lehmann estimate of the median difference (STAY − standard) is reported: the median of all pairwise differences between STAY and standard sections. Its 95% confidence interval was obtained by percentile bootstrap (10,000 resamples of sections within each cohort, fixed random seed).

The primary outcome was the detachment rate. Detached area, compartment-specific detachment rates and RNA quality metrics were secondary outcomes; their P values are reported without adjustment for multiple comparisons.

## Software

Analyses were done in Python **[3.11]** with NumPy, pandas, SciPy (**[version]**), OpenCV, h5py, Panel and HoloViews; the environment is `environment/stay-python.yaml`. The code is available at **[repository URL / DOI]**.
