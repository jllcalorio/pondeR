# Perform Principal Component Analysis

Performs Principal Component Analysis (PCA) on preprocessed metabolomics
data. Optionally applies scaling and transformation before PCA if the
input data has not been pre-scaled. Supports sample exclusion by group.

## Usage

``` r
run_pca(
  x,
  metadata,
  scale_method = NULL,
  transform_method = NULL,
  group = "Group",
  exclude = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns. Can be output from
  [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md)
  or raw data (which will be scaled if `scale_method` is specified).

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x). Must
  contain columns for grouping, batch information, and sample
  identifiers.

- scale_method:

  Character or NULL. Scaling method to apply before PCA if `x` is not
  already scaled. Options: "mean", "auto", "pareto". If NULL (default)
  and `x` is not from
  [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md),
  a warning is issued and PCA proceeds without additional scaling.

- transform_method:

  Character or NULL. Transformation method to apply before scaling and
  PCA. Options: "log2", "log10", "sqrt", "cbrt", "vsn", "glog". Default:
  NULL (no transformation).

- group:

  Character. Name of column in `metadata` containing group labels for
  sample exclusion. Required if `exclude` is specified. Default:
  "Group".

- exclude:

  Character vector or NULL. Group labels (from the `group` column) to
  exclude from PCA analysis. For example, c("QC", "Blank") excludes all
  samples with these group labels. Default: NULL (include all samples).

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_pca" containing:

- scores:

  Matrix of PC scores (samples × PCs)

- loadings:

  Matrix of PC loadings (features × PCs)

- variance_explained:

  Numeric vector of variance explained (%) per PC

- eigenvalues:

  Numeric vector of eigenvalues per PC

- cumulative_variance:

  Numeric vector of cumulative variance explained (%) per PC

- center_values:

  Numeric vector of centering values used (if any)

- scale_values:

  Numeric vector of scaling values used (if any)

- pca_object:

  Original prcomp object from stats::prcomp

- data_used:

  Matrix of data used for PCA (after any transformation/scaling)

- metadata:

  Data frame of metadata aligned with PCA scores

- n_samples:

  Integer. Number of samples used in PCA

- n_samples_excluded:

  Integer. Number of samples excluded

- n_features:

  Integer. Number of features

- n_pcs:

  Integer. Number of PCs computed

- excluded_groups:

  Character vector of excluded group labels

- preprocessing:

  List documenting preprocessing steps applied

- parameters:

  List of parameters used

## Details

**Sample Exclusion:**

The `exclude` parameter allows removal of specific sample groups before
PCA:

- Common use: Exclude QC samples (`exclude = c("QC", "SQC", "EQC")`)

- Multiple groups can be excluded:
  `exclude = c("QC", "Blank", "Failed")`

- Exclusion is performed before any preprocessing steps

- Excluded samples are not included in transformation/scaling
  calculations

**Input Data Requirements:**

The function accepts data in two forms:

1.  **Pre-scaled data** from
    [`run_scale()`](https://jllcalorio.github.io/pondeR/reference/run_scale.md):
    PCA is performed directly without additional preprocessing. This is
    the recommended workflow.

2.  **Raw data**: If `scale_method` is specified, the function will
    first apply optional transformation (via `transform_method`), then
    scaling (via `scale_method`), before performing PCA.

**PCA Methodology:**

PCA is performed using singular value decomposition (SVD) via
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html). The function
assumes data is appropriately preprocessed (normalized, transformed, and
scaled) before analysis. PCA is performed without additional centering
or scaling (`center = FALSE, scale. = FALSE`) as these should have been
applied during the scaling step.

**Interpretation:**

- **PC Scores**: Coordinates of samples in the principal component
  space. Used for visualization and identifying patterns/clusters.

- **PC Loadings**: Contribution of each feature to each principal
  component. Features with high absolute loadings drive separation along
  that PC.

- **Variance Explained**: Percentage of total variance captured by each
  PC. First few PCs typically capture most variance in metabolomics
  data.

- **Eigenvalues**: Variance captured by each PC in original units.
  Eigenvalue \> 1 (Kaiser criterion) often used to determine significant
  PCs.

**Scaling Recommendations:**

- **Auto-scaling** ("auto"): Recommended for PCA. Gives equal weight to
  all features regardless of magnitude.

- **Pareto-scaling** ("pareto"): Alternative that partially preserves
  variance structure. Can be useful when feature magnitudes carry
  biological meaning.

- **Mean-centering** ("mean"): Only centers data. Use when features are
  already on comparable scales.

## References

Becker, R.A., Chambers, J.M., & Wilks, A.R. (1988). The New S Language.
Wadsworth & Brooks/Cole.
[doi:10.1201/9781351074988](https://doi.org/10.1201/9781351074988)

Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). Multivariate Analysis.
London: Academic Press.

Venables, W.N., & Ripley, B.D. (2002). Modern Applied Statistics with S.
Springer-Verlag.

van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A., Smilde, A.K., &
van der Werf, M.J. (2006). Centering, scaling, and transformations:
improving the biological information content of metabolomics data. BMC
Genomics, 7, 142.
[doi:10.1186/1471-2164-7-142](https://doi.org/10.1186/1471-2164-7-142)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(rnorm(100 * 50, mean = 100, sd = 20),
            nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)
rownames(x) <- paste0("Sample", 1:100)

metadata <- data.frame(
  Sample = paste0("Sample", 1:100),
  Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
  Batch = rep(1:4, each = 25),
  Time = rep(c(0, 6, 12, 24), 25)
)

# PCA on pre-scaled data (recommended)
scaled_result <- run_scale(x, method = "auto")
pca_result <- run_pca(scaled_result$data, metadata)

# PCA excluding QC samples
pca_result2 <- run_pca(scaled_result$data, metadata, 
                       group = "Group",
                       exclude = "QC")

# PCA with automatic scaling, excluding multiple groups
pca_result3 <- run_pca(x, metadata, 
                       scale_method = "auto",
                       group = "Group",
                       exclude = c("QC", "Blank"))

# PCA with transformation and scaling
pca_result4 <- run_pca(x, metadata, 
                       transform_method = "log2",
                       scale_method = "pareto",
                       group = "Group",
                       exclude = "QC")

# View results
print(pca_result)
summary(pca_result)
} # }
```
