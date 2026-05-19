# Direct-Injection Metabolomics Preprocessing Pipeline

Runs a complete, sequential preprocessing pipeline for metabolomics data
by orchestrating the individual `run_*` functions of pondeR. Steps
execute in a fixed, reproducible order regardless of argument
arrangement.

## Usage

``` r
run_DIpreprocess(
  x,
  metadata,
  sample_id_col = "Sample",
  group_col = "Group",
  qc_types = c("QC", "SQC", "EQC"),
  batch_col = "Batch",
  injection_col = "InjectionSequence",
  norm_factor_col = "Normalization",
  subject_id_col = "SubjectID",
  outliers = NULL,
  missing_threshold = 0.2,
  missing_by_group = TRUE,
  missing_include_qc = FALSE,
  impute_fraction = 0.2,
  positive_only = TRUE,
  correct_drift = TRUE,
  remove_uncorrected = FALSE,
  spline_smooth = 0,
  spline_spar_limit = c(-1.5, 1.5),
  correct_on_log = TRUE,
  min_qc_per_batch = 5L,
  normalize_method = "sum",
  normalize_ref_sample = NULL,
  normalize_qc_method = "median",
  transform_method = "log2",
  vsn_cores = 1L,
  scale_nonpls = "auto",
  scale_pls = "pareto",
  rsd_threshold = 0.3,
  rsd_qc_type = "EQC",
  variance_percentile = 10,
  scale_filter_ref = "auto",
  merge_replicates = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric feature data with **samples in rows**
  and **features in columns**. All values should be non-negative raw
  intensities; zeros are treated as missing internally.

- metadata:

  Data frame. Sample metadata with `nrow(metadata) == nrow(x)`. At
  minimum, the group column (default `"Group"`) must be present.

- sample_id_col:

  Character. Column in `metadata` containing unique sample identifiers
  that match `rownames(x)`. Default: `"Sample"`.

- group_col:

  Character. Column in `metadata` containing group labels. Default:
  `"Group"`.

- qc_types:

  Character vector. Group labels that identify QC samples. Default:
  `c("QC", "SQC", "EQC")`.

- batch_col:

  Character. Column in `metadata` containing batch numbers. Required for
  drift correction. Default: `"Batch"`.

- injection_col:

  Character. Column in `metadata` containing injection order (integer).
  Required for drift correction. Default: `"InjectionSequence"`.

- norm_factor_col:

  Character. Column in `metadata` containing external normalization
  factors. Default: `"Normalization"`.

- subject_id_col:

  Character. Column in `metadata` containing subject identifiers for
  technical-replicate merging. Default: `"SubjectID"`.

- outliers:

  Character vector or `NULL`. Row names of `x` to remove before
  processing. Default: `NULL`.

- missing_threshold:

  Numeric \[0, 1\]. Missing value threshold. Default: `0.2`.

- missing_by_group:

  Logical. If `TRUE`, assess per group. Default: `TRUE`.

- missing_include_qc:

  Logical. If `TRUE`, include QC samples. Default: `FALSE`.

- impute_fraction:

  Numeric \> 0. Fraction of the smallest positive value per feature used
  for imputation. Default: `0.2`.

- positive_only:

  Logical. If `TRUE`, imputation only considers positive values.
  Default: `TRUE`.

- correct_drift:

  Logical. If `TRUE`, apply QCRSC correction. Default: `TRUE`.

- remove_uncorrected:

  Logical. If `TRUE`, remove features QCRSC could not correct. Default:
  `FALSE`.

- spline_smooth:

  Numeric \[0, 1\]. Smoothing parameter. Default: `0`.

- spline_spar_limit:

  Numeric vector. Lower/upper bounds for spline. Default:
  `c(-1.5, 1.5)`.

- correct_on_log:

  Logical. If `TRUE`, fit on log-transformed data. Default: `TRUE`.

- min_qc_per_batch:

  Integer. Minimum QC samples per batch. Default: `5`.

- normalize_method:

  Character. Normalization method.

  - `"sum"`: Normalizes by sum.

  - `"median"`: Normalizes by median.

  - `"specific_factor"`: Uses external factors in `norm_factor_col`.

  - `"pqn_global"`: PQN using a global reference.

  - `"pqn_reference"`: PQN using `normalize_ref_sample`.

  - `"pqn_group"`: PQN grouping.

  - `"quantile"`: Quantile normalization.

  - `"none"`: No normalization.

  Default: `"sum"`.

- normalize_ref_sample:

  Character. Reference sample name for PQN. Default: `NULL`.

- normalize_qc_method:

  Character. QC normalization when using factors. One of `"mean"`,
  `"median"`, `"none"`. Default: `"median"`.

- transform_method:

  Character. Transformation method.

  - `"log2"`, `"log10"`, `"sqrt"`, `"cbrt"`, `"clr"`, `"vsn"`, `"glog"`,
    or `"none"`.

  Default: `"log2"`.

- vsn_cores:

  Integer or `"max"`. Cores for VSN. Default: `1`.

- scale_nonpls:

  Character. Scaling for NONPLS branch (`"auto"`, `"pareto"`, `"mean"`).

- scale_pls:

  Character. Scaling for PLS branch (`"pareto"`, `"auto"`, `"mean"`).

- rsd_threshold:

  Numeric \[0, 1\]. QC-RSD threshold. Default: `0.3`.

- rsd_qc_type:

  Character. QC type for RSD (`"EQC"`, `"SQC"`, `"QC"`).

- variance_percentile:

  Numeric \[0, 100\]. Percentile to filter. Default: `10`.

- scale_filter_ref:

  Character. Harmonisation strategy.

  - `"auto"`: Keep intersection of both branches.

  - `"NONPLS"`: Apply NONPLS filters to both.

  - `"PLS"`: Apply PLS filters to both.

- merge_replicates:

  Logical. Average technical replicates. Default: `FALSE`.

- verbose:

  Logical. Print progress. Default: `TRUE`.

## Value

A named list of class `run_DIpreprocess` containing:

- `metadata`: Processed sample metadata.

- `data_raw`: Original feature matrix.

- `data_nonpls`: Final NONPLS-scaled, filtered data.

- `data_pls`: Final PLS-scaled, filtered data.

- `features_final`: Character vector of retained features.

- `dimensions`: Sample and feature counts at each step.

- `elapsed_seconds`: Execution time.

## Details

**Preprocessing Workflow (Fixed Order):**

1.  Input validation

2.  Outlier removal

3.  Missing-value filtering
    ([`run_filtermissing`](https://jllcalorio.github.io/pondeR/reference/run_filtermissing.md))

4.  Missing-value imputation
    ([`run_mvimpute`](https://jllcalorio.github.io/pondeR/reference/run_mvimpute.md))

5.  Signal drift correction
    ([`run_driftBatchCorrect`](https://jllcalorio.github.io/pondeR/reference/run_driftBatchCorrect.md))

6.  Normalization
    ([`run_normalize`](https://jllcalorio.github.io/pondeR/reference/run_normalize.md))

7.  Transformation
    ([`run_transform`](https://jllcalorio.github.io/pondeR/reference/run_transform.md))

8.  Scaling
    ([`run_scale`](https://jllcalorio.github.io/pondeR/reference/run_scale.md))

9.  Quality filtering
    ([`run_filterRSD`](https://jllcalorio.github.io/pondeR/reference/run_filterRSD.md),
    [`run_filtervariance`](https://jllcalorio.github.io/pondeR/reference/run_filtervariance.md))

10. Common-feature harmonisation

## Note

**Zero Handling:** The pipeline automatically converts all `0` values in
`x` to `NA` prior to analysis. This represents a strong, domain-specific
assumption (common in mass spectrometry metabolomics) that zeros signify
structural missingness (e.g., values below the limit of detection)
rather than a true biological absence.

## References

Kirwan, J.A., et al. (2013). *Analytical and Bioanalytical Chemistry*,
405, 5147–5157.

## See also

[`run_filtermissing`](https://jllcalorio.github.io/pondeR/reference/run_filtermissing.md),
[`run_mvimpute`](https://jllcalorio.github.io/pondeR/reference/run_mvimpute.md),
[`run_driftBatchCorrect`](https://jllcalorio.github.io/pondeR/reference/run_driftBatchCorrect.md),
[`run_normalize`](https://jllcalorio.github.io/pondeR/reference/run_normalize.md)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
n_s <- 80; n_f <- 100
x <- matrix(abs(rnorm(n_s * n_f, 500, 150)), nrow = n_s, ncol = n_f)
colnames(x) <- paste0("Feature", seq_len(n_f))
x[sample(length(x), 400)] <- 0

meta <- data.frame(
  Sample          = paste0("S", seq_len(n_s)),
  Group           = rep(c("Control", "Treatment", "QC"), c(30, 30, 20)),
  Batch           = rep(1:2, each = 40),
  InjectionSequence = seq_len(n_s),
  SubjectID       = c(paste0("BIO", seq_len(60)), rep(NA, 20)),
  stringsAsFactors = FALSE
)
rownames(x) <- meta$Sample

result <- run_DIpreprocess(
  x              = x,
  metadata       = meta,
  normalize_method = "pqn_group",
  transform_method = "log2",
  scale_nonpls   = "auto",
  scale_pls      = "pareto",
  correct_drift  = FALSE
)

dim(result$data_nonpls)
result$dimensions
} # }
```
