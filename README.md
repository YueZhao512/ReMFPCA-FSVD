# ReMFPCA via functional SVD: MotionSense example

R code for the MotionSense application in **A Functional SVD Framework for Regularized Multivariate Functional PCA with Dual Penalization**.

`code.R` compares non-regularized, smooth, sparse, and smooth-and-sparse multivariate functional principal component analysis. It also applies k-means clustering to the estimated principal component scores. This repository contains the real-data example; the simulation studies are not included.

## Requirements

Install the packages used by the script from CRAN:

```r
install.packages(c(
  "ReMFPCA", "expm", "fda", "ggplot2", "scales",
  "factoextra", "cluster", "aricode", "mclust", "gridExtra"
))
```

The script uses the smoothing and sparsity options in `ReMFPCA` 2.0.0. Check the installed version with:

```r
packageVersion("ReMFPCA")
```

## Data

The cleaned data are included in the [ReMFPCA package](https://CRAN.R-project.org/package=ReMFPCA) as `motion_sense_data`. The script uses two variables, `user_acceleration` and `pitch_attitude`, for 96 observations: 24 each for jogging, standing, walking, and sitting. Each variable is stored as a 200-by-96 matrix, with observations in columns.

No separate data download or local data path is required. The original MotionSense dataset is available from [mmalekzadeh/motion-sense](https://github.com/mmalekzadeh/motion-sense).

## Running the analysis

Download the repository and open `code.R` in R or RStudio. Run the script from top to bottom: later plotting and clustering sections use objects created in the earlier sections. Parameter selection evaluates multiple candidate values and can take some time.

To run the script from the R console, set the working directory to the downloaded repository and use:

```r
source("code.R", print.eval = TRUE)
```

`print.eval = TRUE` also displays the standalone ggplot objects near the end of the clustering section. In RStudio, plots appear in the Plots pane; use its history arrows to view earlier plots.

The first line of `code.R` installs `ReMFPCA`. If the required version is already installed, that line can be skipped when running the script interactively.

## Analysis and outputs

The script:

1. Represents the two functional variables using 120 B-spline basis functions, evaluates them on a grid of 1,000 points, and weights each variable by its total variation.
2. Estimates three principal components under five settings: no regularization, sparsity only, sequential smoothing with sparsity, sequential smoothing, and joint smoothing.
3. Plots the estimated component functions and principal component scores.
4. Examines silhouette widths and runs k-means clustering using the first score, the second score, and both scores together.
5. Plots component perturbations around the mean and compares trajectories grouped by estimated clusters and recorded activities.

The section comments in `code.R` retain the figure labels from the analysis script: Figures 7.1, 7.2, 7.3, C.4, and C.6. Figure numbering may differ in the manuscript version being read.

Fitted models, scores, and cluster assignments remain in the R session. The script draws plots on the active graphics device and does not explicitly export figure files or result tables. To retain a plot, export it from RStudio or run its plotting section inside a graphics device such as `pdf()` and close the device with `dev.off()`.

## Reproducibility notes

The script sets `set.seed(120)` before each of the three reported k-means fits. It does not set a seed before model fitting or silhouette calculations. Cross-validation uses random partitions, so selected tuning parameters and subsequent results may vary between runs. Package versions can also affect numerical results or plot appearance. Record `sessionInfo()` when saving results.

The analysis code is provided as supplied, without changes to its calculations or plotting settings.

## Checked environment

All analysis and plotting sections completed in R 4.3.1 with the package versions below. Packages were already installed, so the installation command on the first line was skipped for this check.

| Package | Version |
| --- | --- |
| ReMFPCA | 2.0.0 |
| expm | 1.0.0 |
| fda | 6.3.0 |
| ggplot2 | 3.5.2 |
| scales | 1.4.0 |
| factoextra | 1.0.7 |
| cluster | 2.1.8.2 |
| aricode | 1.0.3 |
| mclust | 6.1.1 |
| gridExtra | 2.3 |
