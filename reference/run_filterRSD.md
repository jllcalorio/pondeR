# Filter Features by Relative Standard Deviation in QC Samples

Removes features with high relative standard deviation (RSD, also known
as coefficient of variation) in quality control samples. High RSD
indicates poor analytical reproducibility.

## Usage

``` r
run_filterRSD(
  x,
  metadata,
  max_rsd = 0.3,
  groups = "Group",
  qc_type = "EQC",
  qc_types = c("QC", "SQC", "EQC"),
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x). Must
  contain column specified by `groups` to identify QC samples.

- max_rsd:

  Numeric. Maximum allowed RSD (0-1). Features with RSD \>= this
  threshold in QC samples are removed. Default: 0.3 (30%).

- groups:

  Character. Name of column in `metadata` containing sample group
  labels. Default: "Group".

- qc_type:

  Character. QC sample type to use for RSD calculation. Options:

  - `"EQC"`: Extract QC samples (diluted sample pool)

  - `"SQC"`: Sample QC samples (undiluted sample pool)

  - `"QC"`: Any QC sample (use when both EQC and SQC present, or no
    distinction)

  Default: "EQC".

- qc_types:

  Character vector. Group labels that should be treated as QC samples.
  Default: c("QC", "SQC", "EQC").

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_filterRSD" containing:

- data:

  Matrix or data frame of filtered data (same class as input `x`)

- features_removed:

  Character vector of removed feature names

- features_kept:

  Character vector of retained feature names

- rsd_values:

  Numeric vector of RSD values for all original features

- n_features_before:

  Integer. Number of features before filtering

- n_features_after:

  Integer. Number of features after filtering

- n_features_removed:

  Integer. Number of features removed

- parameters:

  List of parameters used

## Details

**Relative Standard Deviation (RSD):**

RSD is calculated as: RSD = (SD / Mean) × 100%

Also known as coefficient of variation (CV), RSD expresses variability
as a percentage of the mean, making it scale-independent and suitable
for comparing features with different magnitudes.

**Why Filter by RSD:**

Quality control samples should have low variability since they represent
the same biological matrix. High RSD in QC samples indicates:

- Poor analytical reproducibility

- Instrument instability

- Ion suppression/enhancement effects

- Features near detection limit

Features with high QC RSD are unreliable and should be removed before
biological interpretation.

**Choosing QC Type:**

- **EQC** (default): Most commonly used. Represents diluted pooled
  samples analyzed regularly throughout the run.

- **SQC**: Undiluted pooled samples. Use when you have a specific SQC
  type distinct from EQC.

- **QC**: Use when you don't distinguish between QC types, or want to
  assess RSD across all QC samples combined.

**Typical RSD Thresholds:**

- Stringent: 20% (0.2) - high-quality data only

- Standard: 30% (0.3) - recommended for most metabolomics studies

- Lenient: 40% (0.4) - when sample size is limited

## References

Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control,
Batch Correction, Exploration & Data Cleaning (Version 2.1) Zendono.
[doi:10.5281/zenodo.16824822](https://doi.org/10.5281/zenodo.16824822) .
Retrieved from <https://github.com/broadhurstdavid/QC-MXP>.

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

# Add high variability to some features in QC samples
qc_rows <- 91:100
x[qc_rows, 1:10] <- x[qc_rows, 1:10] * runif(10 * 10, 0.5, 1.5)

metadata <- data.frame(
  Sample = paste0("S", 1:100),
  Group = c(rep("Control", 40), rep("Treatment", 40), rep("QC", 10), rep("EQC", 10))
)

# Filter by RSD in EQC samples
result <- run_filterRSD(x, metadata, max_rsd = 0.3, qc_type = "EQC")
} # }
```
