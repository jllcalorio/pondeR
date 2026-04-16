# Inspect Group Sample Sizes That May Trigger CV Fold Reduction in `run_regreg()`

Tabulates the number of samples per level of the dependent variable
(`pred`) after applying the same pre-processing steps that
[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
performs: removing rows with `NA` in `pred`, and optionally subsetting
rows to match those that would survive the `not_penalized` NA removal.
This lets the user see exactly why `run_regreg` may warn that `cv_folds`
was reduced — because the smallest group has fewer samples than the
requested number of folds — and which group is the bottleneck.

## Usage

``` r
get_groupsizes(
  metadata,
  pred,
  not_penalized = NULL,
  is_numeric = NULL,
  train_percent = 0.8,
  cv_folds = 10L,
  show_summary_table = FALSE,
  summarytable_args = list()
)
```

## Arguments

- metadata:

  A `data.frame`, `tibble`, **or a named list** of such objects —
  identical to what was (or will be) passed to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
  When a list, names must be set.

- pred:

  A single character string naming the dependent variable column in
  `metadata`.

- not_penalized:

  Optional character vector. When supplied, rows with `NA` in any
  metadata-sourced `not_penalized` column are removed before counting,
  mirroring the exact sample set `run_regreg` would use. Default: `NULL`
  (no additional row removal).

- is_numeric:

  Optional character vector of column names listed in `not_penalized`
  that should be treated as numeric, bypassing the default
  factor/character detection. This mirrors the `is_numeric` argument of
  [`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  and ensures that the NA-filtering step uses the same coercion logic,
  so group sizes are counted on exactly the same sample set. For
  example, if `"Serum albumin (g/L)"` is stored as character but is
  truly numeric, listing it here prevents it from being coerced to
  integer codes — which could introduce spurious NAs or hide real ones.
  Default: `NULL`.

- train_percent:

  A numeric scalar in `(0, 1)` matching the `train_percent` argument
  passed to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
  Used to compute the projected minimum training-group size. Default:
  `0.8`.

- cv_folds:

  Integer. The `cv_folds` value passed to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
  Used to flag groups whose projected training-set size is smaller than
  the requested folds. Default: `10L`.

- show_summary_table:

  Logical. If `TRUE` and `pred` is categorical, a
  [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  summary of `metadata` stratified by `pred` is printed as a
  side-effect. The numeric count data frame is still returned invisibly.
  Default: `FALSE`.

- summarytable_args:

  A named list of additional arguments forwarded to
  [`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  when `show_summary_table = TRUE`. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `data.frame` with the following columns:

- `taxon`:

  Only when `metadata` is a named list.

- `group`:

  Level of the `pred` variable.

- `n_total`:

  Number of samples in that group (after NA removal).

- `n_train_projected`:

  Projected training-set size: `floor(n_total * train_percent)`.

- `folds_requested`:

  The value of `cv_folds`.

- `fold_reduction_expected`:

  Logical; `TRUE` when `n_train_projected < cv_folds`.

- `suggested_cv_folds`:

  The fold count `run_regreg` would actually use:
  `max(3, n_train_projected)`.

## Details

Summarise Group Sizes Relevant to Cross-Validation Fold Reduction

The cross-validation fold reduction warning fires when
`min(table(y_train)) < cv_folds`. Because the training set is a
`train_percent` fraction of the data, the actual minimum group size used
in fold assignment is smaller than what appears in the full dataset.
This function reports both the full-data group sizes and (optionally) a
projected training-set minimum, so the user can anticipate the reduction
before running the model.

An optional publication-ready summary table can be produced via
[`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
when `show_summary_table = TRUE`.

When `metadata` is a named list (parallel to a `lapply` call), each
element is processed and results are stacked with a leading `taxon`
column.

## See also

[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md),
[`get_removed_samples`](https://jllcalorio.github.io/pondeR/reference/get_removed_samples.md),
[`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
## ---- Single metadata data frame ------------------------------------------
gs <- get_groupsizes(
  metadata      = df_metadata_laboratory |>
                    dplyr::filter(`Disease Severity` != "Healthy"),
  pred          = "Disease Severity",
  not_penalized = vars_include_laboratory,
  is_numeric    = "Serum albumin (g/L)",
  train_percent = 0.6,
  cv_folds      = 10L
)
print(gs)

## ---- With a publication-ready summary table ------------------------------
get_groupsizes(
  metadata           = df_metadata_laboratory,
  pred               = "Disease Severity",
  not_penalized      = vars_include_laboratory,
  is_numeric         = "Serum albumin (g/L)",
  train_percent      = 0.6,
  cv_folds           = 10L,
  show_summary_table = TRUE
)

## ---- List of metadata (one per taxon) ------------------------------------
# If each taxon uses the same metadata, just pass it once (not a list).
# If each taxon has its own row-subsetted metadata, wrap them in a list:
gs_list <- get_groupsizes(
  metadata      = setNames(
    lapply(names(df_list_clr), function(nm) df_metadata_laboratory),
    names(df_list_clr)
  ),
  pred          = "Disease Severity",
  not_penalized = vars_include_laboratory,
  train_percent = 0.6,
  cv_folds      = 10L
)
# Taxa where fold reduction is expected
gs_list[gs_list$fold_reduction_expected, ]
} # }
```
