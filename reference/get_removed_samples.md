# Get Samples Removed by `run_regreg()` Due to Missing Covariate Data

Replicates the NA-detection logic of
[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
to identify which samples would be (or were) removed because of missing
values in one or more columns supplied to `not_penalized`. Returns a
data frame of those samples together with the columns responsible, so
the user can inspect, impute, or exclude them deliberately before
re-running the model.

The function operates entirely on the original inputs — no fitted model
object is needed. This means it can be called before
[`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
to pre-screen data, or after to audit a run.

## Usage

``` r
get_removed_samples(
  x,
  metadata,
  not_penalized,
  pred,
  include_metadata_cols = TRUE,
  row_numbers = TRUE
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, `matrix`, **or a named list** of such
  objects — matching what was passed (or would be passed) to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
  When a list is supplied, each element is processed independently.

- metadata:

  A `data.frame` or `tibble` of sample metadata, identical to the one
  passed to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).

- not_penalized:

  A character vector of column names (from `metadata` and/or `x`) that
  were passed to the `not_penalized` argument of
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).

- pred:

  A single character string naming the dependent variable column in
  `metadata` (same as the `pred` argument of
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)).
  Used only to exclude the response column from the covariate NA check.

- include_metadata_cols:

  Logical. If `TRUE` (default), all columns of `metadata` are appended
  to the returned data frame for context. If `FALSE`, only the columns
  specified in `not_penalized` (plus `taxon` when `x` is a list) are
  returned.

- row_numbers:

  Logical. If `TRUE`, a column `.row` containing the original row index
  of each flagged sample is prepended to the result. Default: `TRUE`.

## Value

A `data.frame` of flagged samples. Columns include:

- `taxon`:

  Only present when `x` is a named list. The name of the list element
  from which the sample originated.

- `.row`:

  Original row index in `metadata` (when `row_numbers = TRUE`).

- `.missing_in`:

  A semicolon-separated string naming the `not_penalized` columns that
  are `NA` for that sample.

- Metadata columns:

  All columns of `metadata` (when `include_metadata_cols = TRUE`).

Returns an empty `data.frame` (zero rows) with a message if no samples
are flagged.

## Details

Identify Samples Removed Due to Missing Values in Unpenalized Covariates

The function applies the same filtering sequence as
[`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md):

1.  The dependent variable column (`pred`) is excluded from
    `not_penalized` (it is never a covariate).

2.  Only columns in `not_penalized` that come from `metadata` (i.e., not
    already in `x`) are checked, because columns already in `x` are
    filtered earlier by the feature-NA check.

3.  Rows with at least one `NA` in those columns are flagged.

When `x` is a named list (e.g., from a `lapply` call over taxonomic
levels), the function iterates over each list element and returns a
combined data frame with a leading `taxon` column.

## See also

[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md),
[`get_groupsizes`](https://jllcalorio.github.io/pondeR/reference/get_groupsizes.md)

## Author

John Lennon L. Calorio

## Examples

``` r
set.seed(1)
feat <- data.frame(f1 = rnorm(20), f2 = rnorm(20))
meta <- data.frame(Group = rep(c("A", "B"), 10), Age = c(NA, 30:47, NA))

# Samples 1 and 20 are dropped because Age is missing
get_removed_samples(feat, meta, not_penalized = "Age", pred = "Group")
#>   .row .missing_in Group Age
#> 1    1         Age     A  NA
#> 2   20         Age     B  NA
```
