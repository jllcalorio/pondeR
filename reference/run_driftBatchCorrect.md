# Correct Signal Drift and Batch Effects Using QC Samples

Applies Quality Control-based Robust Spline Correction (QCRSC) or batch
correction using ComBat to remove systematic signal drift and batch
effects in metabolomics data.

## Usage

``` r
run_driftBatchCorrect(
  x,
  metadata,
  perform_correction = TRUE,
  batch_corr_only = FALSE,
  injection_sequence = "InjectionSequence",
  batch_numbers = "Batches",
  groups = "Group",
  qc_label = "QC",
  qc_types = c("QC", "SQC", "EQC"),
  spline_smooth_param = 0,
  min_QC = 5,
  spar_limit = c(-1.5, 1.5),
  log_scale = TRUE,
  use_parametric = TRUE,
  display_plots = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x). Must
  contain columns specified by `injection_sequence`, `batch_numbers`,
  and `groups`.

- perform_correction:

  Logical. If TRUE, perform correction; if FALSE, return original data
  unchanged. Default: TRUE.

- batch_corr_only:

  Logical. If TRUE, perform only batch correction using ComBat (no drift
  correction). If FALSE, perform QCRSC (drift + batch correction).
  Default: FALSE.

- injection_sequence:

  Character. Name of column in `metadata` containing injection order.
  Required for QCRSC. Default: "InjectionSequence".

- batch_numbers:

  Character. Name of column in `metadata` containing batch numbers.
  Required for both QCRSC and ComBat. Default: "Batches".

- groups:

  Character. Name of column in `metadata` containing sample group labels
  (e.g., "Sample", "QC", "SQC", "EQC"). Required for QCRSC. Default:
  "Group".

- qc_label:

  Character. Group label identifying QC samples for QCRSC. Default:
  "QC".

- qc_types:

  Character vector. Group labels that should be converted to `qc_label`
  internally. Default: c("QC", "SQC", "EQC").

- spline_smooth_param:

  Numeric. Smoothing parameter for spline fitting in QCRSC (0 = more
  flexible, 1 = more rigid). Default: 0.

- min_QC:

  Integer. Minimum number of QC samples required per batch for QCRSC.
  Batches with fewer QCs will not be corrected. Default: 5.

- spar_limit:

  Numeric vector of length 2. Limits for spline fitting in QCRSC as
  c(min, max). Default: c(-1.5, 1.5).

- log_scale:

  Logical. If TRUE, fit QCRSC spline on log-transformed data. Default:
  TRUE.

- use_parametric:

  Logical. For ComBat only. Use parametric adjustment for batch effects.
  Default: TRUE.

- display_plots:

  Logical. For ComBat only. Display diagnostic plots during correction.
  Default: FALSE.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_driftBatchCorrect" containing:

- data:

  Matrix or data frame of corrected data (same class as input `x`)

- data_before_correction:

  Original data before correction

- correction_applied:

  Logical indicating if correction was performed

- method_used:

  Character describing correction method ("QCRSC" or "ComBat")

- uncorrected_features:

  Character vector of features that could not be corrected

- n_uncorrected:

  Integer. Number of uncorrected features

- parameters:

  List of parameters used

## Details

**QCRSC (Quality Control-based Robust Spline Correction):**

QCRSC is the default and recommended method for metabolomics data with
regular QC injections. It works by:

1.  Fitting a smoothing spline through QC sample intensities across
    injection order

2.  Using this spline to model systematic drift and batch effects

3.  Normalizing all samples (biological + QC) to remove the fitted drift

**When to use QCRSC:**

- Long analytical runs with visible signal drift

- Multiple batches with sufficient QC samples (≥ `min_QC` per batch)

- QC samples show systematic trends across injection order

**ComBat Batch Correction:**

When `batch_corr_only = TRUE`, uses the ComBat algorithm for batch
effect removal only (no drift correction). This is useful when:

- QC samples are insufficient or absent

- Primary concern is batch differences, not within-batch drift

- Data comes from multiple platforms or sites

ComBat uses empirical Bayes methods to adjust for batch effects while
preserving biological variation.

**Uncorrected Features:**

Some features may remain uncorrected due to:

- Insufficient QC samples (\< `min_QC`) in one or more batches

- Poor spline fit or numerical issues

- All QC values identical (no variance to model)

Uncorrected features are returned unchanged with their names listed in
the output.

## References

Kirwan, J.A., Broadhurst, D.I., Davidson, R.L. et al. (2013).
Characterising and correcting batch variation in an automated direct
infusion mass spectrometry (DIMS) metabolomics workflow. Analytical and
Bioanalytical Chemistry, 405, 5147-5157.
[doi:10.1007/s00216-013-6856-7](https://doi.org/10.1007/s00216-013-6856-7)

Johnson, W.E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects
in microarray expression data using empirical Bayes methods.
Biostatistics, 8(1), 118-127.
[doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

Jankevics A, Lloyd GR, Weber RJM (2025). pmp: Peak Matrix Processing and
signal batch correction for metabolomics datasets.
[doi:10.18129/B9.bioc.pmp](https://doi.org/10.18129/B9.bioc.pmp) , R
package version 1.20.0, <https://bioconductor.org/packages/pmp>.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with drift and batch effects
set.seed(123)
n_samples <- 120
n_features <- 50
x <- matrix(rnorm(n_samples * n_features, mean = 100, sd = 20),
            nrow = n_samples, ncol = n_features)
colnames(x) <- paste0("Feature", 1:n_features)

# Add systematic drift
drift <- seq(0, 30, length.out = n_samples)
x <- x + drift

metadata <- data.frame(
  Sample = paste0("S", 1:n_samples),
  InjectionSequence = 1:n_samples,
  Batches = rep(1:3, each = 40),
  Group = rep(c(rep("Sample", 9), "QC"), n_samples / 10)
)

# Apply QCRSC correction
result1 <- run_driftBatchCorrect(x, metadata)

# Apply ComBat correction only
result2 <- run_driftBatchCorrect(x, metadata, batch_corr_only = TRUE)
} # }
```
