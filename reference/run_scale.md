# Scale Metabolomics Data

Applies various scaling methods to standardize features and prepare data
for multivariate analysis.

## Usage

``` r
run_scale(x, method = "auto", verbose = TRUE)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- method:

  Character. Scaling method to apply. Options:

  - `"mean"`: Mean-centering only (centers distribution at zero)

  - `"auto"`: Auto-scaling (mean-centering + unit variance scaling)

  - `"pareto"`: Pareto-scaling (mean-centering + sqrt(SD) scaling)

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_scale" containing:

- data:

  Matrix or data frame of scaled data (same class as input `x`)

- scaling_factors:

  Numeric vector of scaling factors applied to each feature

- center_values:

  Numeric vector of centering values (means) for each feature

- method_used:

  Character describing scaling method applied

- parameters:

  List of parameters used

## Details

**Scaling Methods:**

- **mean**: Mean-centering only. Each feature is centered to have mean =
  0 but retains original variance. Useful when features are already on
  similar scales and you want to preserve relative variance information.

- **auto**: Auto-scaling (also called standardization or unit variance
  scaling). Each feature is mean-centered and divided by its standard
  deviation, resulting in mean = 0 and SD = 1 for all features. This
  gives equal weight to all features regardless of their original
  variance.

  **Recommended for:**

  - Principal Component Analysis (PCA)

  - Hierarchical clustering

  - Correlation analysis

  - t-tests and ANOVA

  - Most univariate and multivariate methods

  Auto-scaling can amplify noise in low-variance features.

- **pareto**: Pareto-scaling. Each feature is mean-centered and divided
  by the square root of its standard deviation. This provides a
  compromise between no scaling and auto-scaling, reducing the influence
  of large values while preserving some of the original variance
  structure.

  **Recommended for:**

  - PLS-DA (Partial Least Squares Discriminant Analysis)

  - OPLS-DA (Orthogonal PLS-DA)

  - sPLS-DA (Sparse PLS-DA)

  - Other PLS-type methods

  Pareto-scaling is generally preferred for discriminant analysis as it
  balances the importance of large and small features better than
  auto-scaling.

**Choosing a Scaling Method:**

The choice depends on your downstream analysis:

- For exploratory data analysis (PCA, clustering): **auto-scaling**

- For classification/discrimination (PLS-DA): **Pareto-scaling**

- When feature magnitudes carry meaning: **mean-centering**

## References

van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A., Smilde, A.K., &
van der Werf, M.J. (2006). Centering, scaling, and transformations:
improving the biological information content of metabolomics data. BMC
Genomics, 7, 142.
[doi:10.1186/1471-2164-7-142](https://doi.org/10.1186/1471-2164-7-142)

Becker, R.A., Chambers, J.M., & Wilks, A.R. (1988). The New S Language.
Wadsworth & Brooks/Cole.
[doi:10.1201/9781351074988](https://doi.org/10.1201/9781351074988)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(rnorm(100 * 50, mean = 100, sd = 50),
            nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)

# Auto-scaling for PCA
result1 <- run_scale(x, method = "auto")

# Pareto-scaling for PLS-DA
result2 <- run_scale(x, method = "pareto")

# Mean-centering only
result3 <- run_scale(x, method = "mean")
} # }
```
