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
  verbose = FALSE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns. Can be output from
  [`run_scale()`](https://jllcalorio.github.io/dimsprepr/reference/run_scale.html)
  or raw data (which will be scaled if `scale_method` is specified).

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x). Must
  contain columns for grouping, batch information, and sample
  identifiers.

- scale_method:

  Character or NULL. Scaling method to apply before PCA if `x` is not
  already scaled. Options: "mean", "auto", "pareto". If NULL (default)
  and `x` is not from
  [`run_scale()`](https://jllcalorio.github.io/dimsprepr/reference/run_scale.html),
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

  Logical. Print progress messages. Default: FALSE.

## Value

A list containing:

- scores:

  Matrix of PC scores (samples x PCs)

- loadings:

  Matrix of PC loadings (features x PCs)

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
    [`run_scale()`](https://jllcalorio.github.io/dimsprepr/reference/run_scale.html):
    PCA is performed directly without additional preprocessing. This is
    the recommended workflow.

2.  **Raw data**: If `scale_method` is specified, the function will
    first apply optional transformation (via `transform_method`), then
    scaling (via `scale_method`), before performing PCA.

**PCA Methodology:**

PCA is performed using singular value decomposition (SVD) via
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html). The function
assumes data is appropriately preprocessed (normalized, transformed, and
scaled) before analysis. Columns are always mean-centered (a no-op for
data that is already centered) and never rescaled here
(`center = TRUE, scale. = FALSE`): scaling belongs to the preprocessing
step.

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

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html)

## Author

John Lennon L. Calorio

## Examples

``` r
meta <- data.frame(Sample = rownames(iris), Species = iris$Species)

# PCA on pre-scaled data (base scale() output is detected as pre-scaled)
pca <- run_pca(scale(iris[, 1:4]), meta, group = "Species", verbose = FALSE)
print(pca)
#> === Principal Component Analysis ===
#> Samples:      150 
#> Features:     4 
#> PCs computed: 4 
#> 
#> Variance explained by first 5 PCs:
#>   PC1: 73.0% (cumulative: 73.0%)
#>   PC2: 22.9% (cumulative: 95.8%)
#>   PC3: 3.7% (cumulative: 99.5%)
#>   PC4: 0.5% (cumulative: 100.0%)
#> 
#> Use plot_scree() to visualize variance explained
#> Use plot_score() to visualize sample positions
summary(pca)
#> $n_samples
#> [1] 150
#> 
#> $n_samples_excluded
#> [1] 0
#> 
#> $excluded_groups
#> character(0)
#> 
#> $n_features
#> [1] 4
#> 
#> $n_pcs
#> [1] 4
#> 
#> $variance_table
#>    PC Variance explained (%) Cumulative (%) Eigenvalue
#> 1 PC1             72.9624454       72.96245 2.91849782
#> 2 PC2             22.8507618       95.81321 0.91403047
#> 3 PC3              3.6689219       99.48213 0.14675688
#> 4 PC4              0.5178709      100.00000 0.02071484
#> 
#> $preprocessing
#> $preprocessing$came_from_run_scale
#> [1] TRUE
#> 
#> $preprocessing$transformation_applied
#> [1] FALSE
#> 
#> $preprocessing$scaling_applied
#> [1] FALSE
#> 
#> $preprocessing$transform_method
#> [1] NA
#> 
#> $preprocessing$scale_method_used
#> [1] NA
#> 
#> 

# Exclude a group before fitting
pca2 <- run_pca(scale(iris[, 1:4]), meta, group = "Species",
                exclude = "setosa", verbose = FALSE)

# In-function scaling/transformation needs the dimsprepr package
if (requireNamespace("dimsprepr", quietly = TRUE)) {
  pca3 <- run_pca(iris[, 1:4], meta, group = "Species",
                  scale_method = "pareto", verbose = FALSE)
}
```
