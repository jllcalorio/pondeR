# Impute Missing Values in Metabolomics Data

Replaces missing values (NAs) using either a simple deterministic method
(fraction of minimum) or a statistical imputation method (quantile
regression).

## Usage

``` r
run_mvimpute(x, method = 0.2, tune_sigma = 1, verbose = TRUE)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns. Missing values should be represented as NA.

- method:

  Numeric or character. If numeric, represents the fraction of the
  minimum positive value per feature to use for imputation (e.g., 0.2 =
  1/5th of minimum). If character, must be one of the available
  imputation methods. Currently supported: "quantileregression" (uses
  QRILC algorithm). Default: 0.2 (1/5th of minimum).

- tune_sigma:

  Numeric. Only used when `method = "quantileregression"`. Controls the
  standard deviation of the left-censored missing value distribution.
  Smaller values (e.g., 0.5) assume tighter distribution; larger values
  (e.g., 2) allow more spread. Default: 1.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_mvimpute" containing:

- data:

  Matrix or data frame with imputed values (same class as input `x`)

- n_missing_before:

  Integer. Number of missing values before imputation

- n_missing_after:

  Integer. Number of missing values after imputation (should be 0)

- imputed_summary:

  Data frame with per-feature imputation statistics

- parameters:

  List of parameters used

- method_used:

  Character describing the imputation method applied

## Details

**Missing Value Imputation Strategies:**

1.  **Deterministic (fraction of minimum)**: When `method` is numeric,
    each missing value in a feature is replaced by the minimum positive
    value in that feature multiplied by the specified fraction. This
    assumes missing values represent low abundance (Missing Not At
    Random - MNAR) and is common in metabolomics for values below
    detection limit.

    - **Recommended range**: 0.01 to 0.5 (1% to 50% of minimum)

    - Values \> 1 will trigger a warning as they exceed the observed
      minimum

    - Default (0.2) replaces NAs with 1/5th of the feature's minimum

2.  **Quantile Regression (QRILC)**: When
    `method = "quantileregression"`, uses the Quantile Regression
    Imputation of Left-Censored data algorithm. This method assumes
    missing values are Missing Not At Random (MNAR) due to low abundance
    and models them from a left-censored distribution.

    - The `tune_sigma` parameter controls the spread of imputed values

    - Requires the `imputeLCMD` package

    - More sophisticated than simple fraction-of-minimum

## References

Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
Missing Value Imputation Approach for Mass Spectrometry-based
Metabolomics Data. Scientific Reports, 8(1), 663.
[doi:10.1038/s41598-017-19120-0](https://doi.org/10.1038/s41598-017-19120-0)

Lazar, C., Gatto, L., Ferro, M., Bruley, C., & Burger, T. (2016).
Accounting for the Multiple Natures of Missing Values in Label-Free
Quantitative Proteomics Data Sets to Compare Imputation Strategies.
Journal of Proteome Research, 15(4), 1116-1125.
[doi:10.1021/acs.jproteome.5b00981](https://doi.org/10.1021/acs.jproteome.5b00981)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with missing values
set.seed(123)
x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
x[sample(length(x), 200)] <- NA
colnames(x) <- paste0("Feature", 1:50)

# Simple imputation (1/5th of minimum)
result1 <- run_mvimpute(x, method = 0.2)

# Quantile regression imputation
result2 <- run_mvimpute(x, method = "quantileregression", tune_sigma = 1)
} # }
```
