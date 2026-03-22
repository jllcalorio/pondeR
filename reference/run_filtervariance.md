# Filter Features by Low Variance

Removes features with low variance across samples. Features with minimal
variation provide little information for downstream analysis and can be
safely removed to reduce dimensionality.

## Usage

``` r
run_filtervariance(x, percentile = 10, verbose = TRUE)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- percentile:

  Numeric. Percentile threshold (0-100) for variance filtering. Features
  in the bottom X percentile of variance are removed. For example, 10
  removes the 10% of features with lowest variance. Default: 10.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_filtervariance" containing:

- data:

  Matrix or data frame of filtered data (same class as input `x`)

- features_removed:

  Character vector of removed feature names

- features_kept:

  Character vector of retained feature names

- variance_values:

  Numeric vector of variance values for all original features

- variance_threshold:

  Numeric. The variance value corresponding to the percentile cutoff

- n_features_before:

  Integer. Number of features before filtering

- n_features_after:

  Integer. Number of features after filtering

- n_features_removed:

  Integer. Number of features removed

- parameters:

  List of parameters used

## Details

**Why Filter by Variance:**

Features with very low variance across samples:

- Provide minimal discriminatory power for classification

- Contribute little to multivariate models (PCA, PLS-DA)

- May represent technical noise or detection artifacts

- Increase computational burden without adding information

**Choosing the Percentile:**

- **Conservative** (5-10%): Remove only the least variable features

- **Moderate** (10-20%): Standard approach for most studies

- **Aggressive** (20-30%): When high dimensionality is a concern

The default (10%) removes the bottom 10% of features by variance, which
is a commonly used threshold in metabolomics.

**Important Notes:**

1.  Variance filtering should be applied **after** normalization and
    transformation, as these steps affect feature variance.

2.  The variance calculation uses all samples (including QCs if
    present). Consider removing QC samples before variance filtering if
    they have artificially low variance.

3.  For scaled data (auto or Pareto), this filter identifies features
    with inherently low biological variation across the study groups.

4.  Very high percentiles (\>30%) risk removing biologically relevant
    features with naturally stable expression levels.

## References

Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control,
Batch Correction, Exploration & Data Cleaning (Version 2.1) Zendono.
[doi:10.5281/zenodo.16824822](https://doi.org/10.5281/zenodo.16824822) .
Retrieved from <https://github.com/broadhurstdavid/QC-MXP>.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with varying levels of variance
set.seed(123)
n_samples <- 100
n_features <- 50

# Create features with different variance levels
x <- matrix(nrow = n_samples, ncol = n_features)
for (i in 1:n_features) {
  # Variance increases with feature index
  sd_val <- i / 2
  x[, i] <- rnorm(n_samples, mean = 100, sd = sd_val)
}
colnames(x) <- paste0("Feature", 1:n_features)

# Remove bottom 10% (5 features with lowest variance)
result1 <- run_filtervariance(x, percentile = 10)

# More aggressive filtering (bottom 20%)
result2 <- run_filtervariance(x, percentile = 20)
} # }
```
