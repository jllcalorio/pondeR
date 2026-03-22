# Filter Features by Missing Value Threshold

Removes features (columns) that exceed a specified missing value
threshold. Can assess missingness globally or by group, with optional
inclusion of QC samples. Supports treating zero values as missing, which
is common in metabolomics data.

## Usage

``` r
run_filtermissing(
  x,
  metadata,
  threshold = 0.2,
  filter_by_group = TRUE,
  include_QC = FALSE,
  group_col = "Group",
  qc_types = c("QC", "SQC", "EQC"),
  zero_as_missing = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns. Missing values should be represented as NA. Zero values
  can optionally be treated as missing via the `zero_as_missing`
  parameter.

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x). Must
  contain a column specified by `group_col` when
  `filter_by_group = TRUE`.

- threshold:

  Numeric. Maximum proportion of missing values allowed (0-1). Features
  with missing proportion \>= threshold are removed. Default: 0.2 (20%).

- filter_by_group:

  Logical. If TRUE, assess missingness within each group separately. A
  feature is kept if it passes the threshold in at least one group.
  Default: TRUE.

- include_QC:

  Logical. If TRUE, include QC samples in missingness calculation.
  Default: FALSE (QC samples are excluded from assessment).

- group_col:

  Character. Name of the column in `metadata` containing group
  information. Required when `filter_by_group = TRUE`. Default: "Group".

- qc_types:

  Character vector. Group labels that identify QC samples. These samples
  are excluded from missingness calculation when `include_QC = FALSE`.
  Default: c("QC", "SQC", "EQC").

- zero_as_missing:

  Logical. If TRUE, treat zero values as missing (NA) before applying
  the filter. This is common in metabolomics where zeros often represent
  values below detection limit rather than true zeros. Default: TRUE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_filtermissing" containing:

- data:

  Matrix or data frame of filtered data (same class as input `x`)

- features_removed:

  Character vector of removed feature names

- features_kept:

  Character vector of retained feature names

- n_features_before:

  Integer. Number of features before filtering

- n_features_after:

  Integer. Number of features after filtering

- n_features_removed:

  Integer. Number of features removed

- n_zeros_converted:

  Integer. Number of zero values converted to NA (if applicable)

- parameters:

  List of parameters used

- missingness_summary:

  Data frame with per-feature missingness statistics

## Details

This function filters features based on the proportion of missing values
(NAs).

**Zero values treated as missing**

In metabolomics data, zero values often represent one of the following:

- Peak intensities below the detection limit

- Failed peak integration

- Metabolites not detected in the sample

When `zero_as_missing = TRUE` (default), zero values are converted to
missing values before filtering. This is recommended for most untargeted
metabolomics datasets.

Set `zero_as_missing = FALSE` only if one of the following is true:

- Your preprocessing pipeline already converted zeros to NA

- Zero represents a true biological absence with quantitative meaning

- Your dataset contains negative values (for example after baseline
  correction)

**Filtering strategies**

Two filtering modes are available:

- **Global filtering** (`filter_by_group = FALSE`): Missingness is
  calculated across all eligible samples. Features with missing
  proportion greater than or equal to `threshold` are removed.

- **Group-wise filtering** (`filter_by_group = TRUE`): Missingness is
  calculated separately within each group. A feature is retained if it
  satisfies the threshold in at least one group. This approach is less
  stringent and helps preserve group-specific features.

**QC sample handling**

When `include_QC = FALSE`, samples whose group labels match `qc_types`
are excluded from missingness calculations. This is recommended because
QC samples often exhibit different missingness patterns compared with
biological samples.

## References

Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control,
Batch Correction, Exploration & Data Cleaning (Version 2.1) Zendono.
[doi:10.5281/zenodo.16824822](https://doi.org/10.5281/zenodo.16824822) .
Retrieved from <https://github.com/broadhurstdavid/QC-MXP>.

Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
Missing Value Imputation Approach for Mass Spectrometry-based
Metabolomics Data. Scientific Reports, 8(1), 663.
[doi:10.1038/s41598-017-19120-0](https://doi.org/10.1038/s41598-017-19120-0)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with zeros and NAs
set.seed(123)
x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
x[sample(length(x), 300)] <- 0    # Add zeros
x[sample(length(x), 200)] <- NA   # Add NAs
colnames(x) <- paste0("Feature", 1:50)

metadata <- data.frame(
  Sample = paste0("S", 1:100),
  Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20))
)

# Treat zeros as missing (default)
result1 <- run_filtermissing(x, metadata, threshold = 0.3)

# Keep zeros as valid values
result2 <- run_filtermissing(x, metadata, threshold = 0.3, 
                              zero_as_missing = FALSE)

# Group-wise filtering with zeros as missing
result3 <- run_filtermissing(x, metadata, threshold = 0.3, 
                              filter_by_group = TRUE,
                              zero_as_missing = TRUE)
} # }
```
