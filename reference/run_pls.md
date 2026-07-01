# Perform Partial Least Squares (PLS) Analysis

Performs Partial Least Squares Discriminant Analysis (PLS-DA),
Orthogonal PLS-DA (OPLS-DA), or sparse PLS-DA (sPLS-DA) on preprocessed
metabolomics data. Optionally applies scaling and transformation before
analysis if the input data has not been pre-scaled. Designed to match
the structure of
[`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md).

## Usage

``` r
run_pls(
  x,
  metadata,
  method = "oplsda",
  group = "Group",
  exclude = NULL,
  scale_method = NULL,
  transform_method = NULL,
  ncomp = 2,
  orthoI = NA,
  crossvalI = 10,
  permI = 20,
  keepX = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns. Can be output from `run_scale()`, raw data (which will be
  scaled if `scale_method`, or from `run_DIpreprocess`. is specified).

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x).
  Ignored, even if supplied, when x is `run_DIpreprocess` class. Must
  contain columns for grouping and sample identifiers.

- method:

  Character. The PLS method to use. Options are `"oplsda"` (default),
  `"plsda"`, or `"splsda"`.

- group:

  Character. Name of column in `metadata` containing group labels for
  classification and sample exclusion. Required parameter. Default:
  "Group".

- exclude:

  Character vector or NULL. Group labels (from the `group` column) to
  exclude from the analysis. For example, c("QC", "Blank") excludes all
  samples with these group labels. Default: NULL (include all samples).

- scale_method:

  Character or NULL. Scaling method to apply before PLS if `x` is not
  already scaled. Options: "mean", "auto", "pareto". If NULL (default),
  the function assumes the data is appropriately scaled. Pareto scaling
  is generally recommended for PLS methods.

- transform_method:

  Character or NULL. Transformation method to apply before scaling and
  PLS. Options: "log2", "log10", "sqrt", "cbrt", "vsn", "glog". Default:
  NULL (no transformation).

- ncomp:

  Integer. Number of predictive components to compute for PLS-DA and
  sPLS-DA. Default: 2.

- orthoI:

  Integer or NA. Number of orthogonal components for OPLS-DA. If NA
  (default), the optimal number is determined automatically by
  cross-validation.

- crossvalI:

  Integer. Number of cross-validation folds. Default: 10.

- permI:

  Integer. Number of permutations for model validation. Default: 20.

- keepX:

  Integer vector. For sPLS-DA only: specifies the number of variables to
  keep on each component. Default: NULL (automatic selection or keeps
  all).

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list containing:

- scores:

  Matrix of scores (samples × components), aliased as PC1, PC2... for
  `plot_score` compatibility

- loadings:

  Matrix of feature loadings (features × components)

- variance_explained:

  Numeric vector of variance explained (%) per component

- vip_scores:

  Numeric vector of Variable Importance in Projection (VIP) scores
  (PLS-DA/OPLS-DA only)

- pls_object:

  Original model object from `ropls` or `mixOmics`

- data_used:

  Matrix of data used for PLS (after any transformation/scaling)

- metadata:

  Data frame of metadata aligned with scores

- n_samples:

  Integer. Number of samples used

- n_samples_excluded:

  Integer. Number of samples excluded

- n_features:

  Integer. Number of features

- n_pcs:

  Integer. Number of components computed

- excluded_groups:

  Character vector of excluded group labels

- preprocessing:

  List documenting preprocessing steps applied

- parameters:

  List of parameters used

## Details

**Compatibility:**

The output of this function is specifically designed to be compatible
with
[`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md).
Component names in the `scores` matrix are aliased as `PC1`, `PC2`,
etc., to allow seamless plotting. For OPLS-DA, "PC1" represents the
Predictive Component, and "PC2" represents the first Orthogonal
Component.

**Methods:**

- **oplsda**: Orthogonal PLS-DA separates predictive variation
  (correlated to class) from orthogonal variation (uncorrelated).
  Excellent for biomarker discovery. Utilizes
  [`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html).

- **plsda**: Standard PLS-DA. Maximizes covariance between features and
  class labels. Utilizes
  [`ropls::opls()`](https://rdrr.io/pkg/ropls/man/opls.html).

- **splsda**: Sparse PLS-DA incorporates L1 penalization to perform
  feature selection during model building. Utilizes
  [`mixOmics::splsda()`](https://rdrr.io/pkg/mixOmics/man/splsda.html).

**Sample Exclusion:**

QC samples should typically be excluded from discriminant analysis model
building. Use the `exclude` parameter (e.g., `exclude = "QC"`) to drop
them prior to scaling and fitting.

## References

Thevenot, E. A., Roux, A., Xu, Y., Ezan, E., & Junot, C. (2015).
Analysis of the Human Adult Urinary Metabolome Variations with Age, Body
Mass Index, and Gender by Implementing a Comprehensive Workflow for
Univariate and OPLS Statistical Analyses. Journal of Proteome Research,
14(8), 3322-3335.

Rohart F, Gautier B, Singh A, Lê Cao K-A. (2017) mixOmics: An R package
for 'omics feature selection and multiple data integration. PLoS Comput
Biol 13(11): e1005752.

## See also

`run_DIpreprocess`

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(rnorm(100 * 50, mean = 100, sd = 20), nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)
rownames(x) <- paste0("Sample", 1:100)

metadata <- data.frame(
  Sample = paste0("Sample", 1:100),
  Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20))
)

# Run PLS-DA excluding QC samples, with automatic pareto scaling
pls_result <- run_pls(x, metadata, 
                      method = "plsda",
                      scale_method = "pareto",
                      group = "Group",
                      exclude = "QC")

# Visualize with plot_score
plot_score(pls_result, color_by = "Group", points_from = "Sample")

# Run OPLS-DA
opls_result <- run_pls(x, metadata, 
                       method = "oplsda",
                       scale_method = "pareto",
                       group = "Group",
                       exclude = "QC")
} # }
```
