# Methods: Figure 2 (STAY vs standard detachment)

## Detachment annotation

Detachment was annotated manually in all STAY (n = 33) and standard (n = 20) sections, using the CytAssist and high-resolution H&E images from the Space Ranger output (v2.0.1; 10x Genomics, 2022). Processing used Python 3.11.18 with NumPy 2.4.6 (Harris et al., 2020) and pandas 2.3.3 (McKinney, 2010).

The H&E image was registered to the CytAssist image using the Visium spots as landmarks: CytAssist spot positions came from the manual Loupe Browser alignment (v8.0.0; 10x Genomics, 2024), which records each spot's array row, column and pixel coordinates, and H&E spot positions from the Space Ranger spot coordinates scaled by the high-resolution scale factor. Spots were paired by array row and column, and an affine H&E-to-CytAssist transformation was estimated by RANSAC (Fischler & Bolles, 1981; reprojection threshold 5 pixels) in OpenCV 4.12.0 (Bradski, 2000) and applied to the H&E image. The registration served only for display; all measurements were made in CytAssist coordinates.

Annotation was performed in a browser-based viewer built with Panel 1.9.4 and HoloViews 1.23.2 (Stevens et al., 2015), with the registered H&E image overlaid on the CytAssist image at adjustable opacity. Detachment regions on the CytAssist image were outlined as detached polygons while cross-checking with the H&E image. Sections were presented in random order with sample identity and cohort hidden. 

## Detachment quantification

Polygons were rasterised onto the CytAssist grid, and a spot was classified as detached if its centre lay within a polygon. As Space Ranger tissue detection can exclude detached tissue, tissue spots were defined as spots classified as tissue by Space Ranger plus detached spots, taken from the unfiltered feature-barcode matrix (read with h5py 3.16.0).

Per section, the detachment rate was the percentage of tissue spots classified as detached, and the detached area the area enclosed by the polygons, calculated as the number of CytAssist pixels inside them (overlaps counted once) multiplied by the area of one pixel. The pixel size (µm per pixel) was the Visium centre-to-centre spot distance (100 µm) divided by the median distance in pixels between adjacent spots in an array row. Detachment rates were also calculated separately for epithelial and non-epithelial spots, with epithelial spots defined by pathologist annotation.

## RNA quality metrics

The number of detected genes, the number of UMIs and the percentage of UMIs from mitochondrial genes (names beginning with "MT-") were calculated per tissue spot from the unfiltered matrix using SciPy 1.17.1 (Virtanen et al., 2020) and summarised as the median per section.

## Statistics

Cohorts were compared by two-sided Mann–Whitney U tests (Mann & Whitney, 1947) in SciPy. Effect sizes are Hodges–Lehmann estimates of the difference in location (STAY − standard), the median of all pairwise differences between STAY and standard sections (Hodges & Lehmann, 1963), with 95% confidence intervals from a percentile bootstrap (Efron, 1979; 10,000 resamples of sections within each cohort). 

## References

- 10x Genomics (2022). Space Ranger v2.0.1. https://www.10xgenomics.com/support/software/space-ranger
- 10x Genomics (2024). Loupe Browser v8.0.0. https://www.10xgenomics.com/support/software/loupe-browser
- Bradski, G. (2000). The OpenCV Library. *Dr. Dobb's Journal of Software Tools*.
- Efron, B. (1979). Bootstrap methods: another look at the jackknife. *The Annals of Statistics*, 7(1), 1–26.
- Fischler, M. A. & Bolles, R. C. (1981). Random sample consensus: a paradigm for model fitting with applications to image analysis and automated cartography. *Communications of the ACM*, 24(6), 381–395.
- Harris, C. R. et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.
- Hodges, J. L. & Lehmann, E. L. (1963). Estimates of location based on rank tests. *The Annals of Mathematical Statistics*, 34(2), 598–611.
- Mann, H. B. & Whitney, D. R. (1947). On a test of whether one of two random variables is stochastically larger than the other. *The Annals of Mathematical Statistics*, 18(1), 50–60.
- McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61.
- Stevens, J.-L. R., Rudiger, P. & Bednar, J. A. (2015). HoloViews: building complex visualizations easily for reproducible science. *Proceedings of the 14th Python in Science Conference*, 59–66.
- Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
