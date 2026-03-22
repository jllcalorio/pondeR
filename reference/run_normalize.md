# Normalize Metabolomics Data

Applies various normalization methods to account for dilution effects
and sample-to-sample variation in metabolomics data.

## Usage

``` r
run_normalize(
  x,
  metadata,
  method = "sum",
  factor_col = "Normalization",
  ref_sample = NULL,
  group_sample = "QC",
  qc_normalize = "median",
  groups = "Group",
  qc_types = c("QC", "SQC", "EQC"),
  reference_method = "mean",
  sample_id_col = "Sample",
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x). May
  contain normalization factors, group labels, and QC sample
  identifiers.

- method:

  Character or numeric vector. Normalization method to apply. Options:

  - `"sum"`: Total sum normalization

  - `"median"`: Median normalization

  - `"specific_factor"`: Use values from a metadata column (specify via
    `factor_col`)

  - `"pqn_global"`: Probabilistic Quotient Normalization using global
    median

  - `"pqn_reference"`: PQN using a specific reference sample (specify
    via `ref_sample`)

  - `"pqn_group"`: PQN using pooled QC samples as reference

  - `"quantile"`: Quantile normalization

  - `"col_rel_abundance"`: Relative abundance per column (feature); each
    value divided by its column sum so each feature sums to 1 across all
    samples

  - `"row_rel_abundance"`: Relative abundance per row (sample); each
    value divided by its row sum so each sample sums to 1 across all
    features

  - Numeric vector: Custom normalization factors (length must equal
    nrow(x) excluding QCs)

- factor_col:

  Character. Name of column in `metadata` containing normalization
  factors when `method = "specific_factor"`. Default: "Normalization".

- ref_sample:

  Character. Sample name to use as reference for
  `method = "pqn_reference"`. Must match a value in the sample
  identifier column of metadata. Default: NULL.

- group_sample:

  Character. Group label for `method = "pqn_group"`. Uses samples with
  this group label as reference. Default: "QC".

- qc_normalize:

  Character. How to normalize QC samples when using specific factors.
  Options: "mean", "median", "none". Default: "median".

- groups:

  Character. Name of column in `metadata` containing group labels.
  Required for PQN and QC identification. Default: "Group".

- qc_types:

  Character vector. Group labels identifying QC samples. Default:
  c("QC", "SQC", "EQC").

- reference_method:

  Character. Method for computing reference in quantile normalization.
  Options: "mean", "median". Default: "mean".

- sample_id_col:

  Character. Name of column in `metadata` containing sample identifiers
  for matching `ref_sample`. Default: "Sample".

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_normalize" containing:

- data:

  Matrix or data frame of normalized data (same class as input `x`)

- normalization_factors:

  Numeric vector of factors used for normalization

- method_used:

  Character describing normalization method applied

- parameters:

  List of parameters used

## Details

**Normalization Methods:**

- **sum**: Divides each sample by its total sum. Assumes similar total
  metabolite abundance across samples.

- **median**: Divides each sample by its median value. More robust to
  outliers than sum.

- **specific_factor**: Uses externally measured factors (e.g.,
  osmolality for urine, protein concentration for plasma) from a
  metadata column. QC samples are normalized separately using
  mean/median of biological sample factors.

- **pqn_global**: Probabilistic Quotient Normalization using the global
  median spectrum as reference. Robust method that accounts for dilution
  effects.

- **pqn_reference**: PQN using a specific reference sample. Useful when
  one sample represents the "ideal" composition.

- **pqn_group**: PQN using pooled QC samples as reference. Recommended
  for metabolomics data with regular QC injections.

- **quantile**: Forces all samples to have identical distributions. Very
  strong normalization that may remove biological variation - use with
  caution.

- **col_rel_abundance**: Computes relative abundance per feature
  (column). Each value is divided by the sum of its column, so each
  feature's values sum to 1 across all samples. Useful for comparing the
  contribution of a feature across samples. Note that the normalization
  factors returned are the column sums.

- **row_rel_abundance**: Computes relative abundance per sample (row).
  Each value is divided by the sum of its row, so each sample's values
  sum to 1 across all features. Equivalent to `"sum"` normalization but
  expressed as proportions. Useful for compositional data analysis. Note
  that the normalization factors returned are the row sums.

**Custom Normalization Factors:**

You can provide a numeric vector as `method` with custom normalization
factors. The vector length must equal the number of biological samples
(excluding QCs). QC samples will be normalized using the mean or median
of these factors.

## References

Dieterle, F., Ross, A., Schlotterbeck, G., & Senn, H. (2006).
Probabilistic Quotient Normalization as Robust Method to Account for
Dilution of Complex Biological Mixtures. Analytical Chemistry, 78(13),
4281-4290. [doi:10.1021/ac051632c](https://doi.org/10.1021/ac051632c)

Jankevics A, Lloyd GR, Weber RJM (2025). pmp: Peak Matrix Processing and
signal batch correction for metabolomics datasets.
[doi:10.18129/B9.bioc.pmp](https://doi.org/10.18129/B9.bioc.pmp) , R
package version 1.20.0, <https://bioconductor.org/packages/pmp>.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(abs(rnorm(100 * 50, mean = 100, sd = 20)),
            nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)

metadata <- data.frame(
  Sample = paste0("S", 1:100),
  Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
  Normalization = runif(100, 0.8, 1.2)
)

# Sum normalization
result1 <- run_normalize(x, metadata, method = "sum")

# PQN with QC reference
result2 <- run_normalize(x, metadata, method = "pqn_group")

# Specific factors from metadata
result3 <- run_normalize(x, metadata, method = "specific_factor",
                         factor_col = "Normalization")

# Relative abundance per feature (column)
result4 <- run_normalize(x, metadata, method = "col_rel_abundance")

# Relative abundance per sample (row)
result5 <- run_normalize(x, metadata, method = "row_rel_abundance")
} # }
```
