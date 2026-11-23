# Inspect Group Sample Sizes That May Trigger CV Fold Reduction in [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)

Tabulates the number of samples per level of the dependent variable
(`pred`) after applying the same pre-processing steps that
[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
performs: removing rows with `NA` in `pred`, and optionally subsetting
rows to match those that would survive the `not_penalized` NA removal.
This lets the user see exactly why
[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
may warn that `cv_folds` was reduced - because the smallest group has
fewer samples than the requested number of folds - and which group is
the bottleneck.

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

  A `data.frame`, `tibble`, **or a named list** of such objects -
  identical to what was (or will be) passed to
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md).
  When a list, names must be set.

- pred:

  A single character string naming the dependent variable column in
  `metadata`.

- not_penalized:

  Optional character vector. When supplied, rows with `NA` in any
  metadata-sourced `not_penalized` column are removed before counting,
  mirroring the exact sample set
  [`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  would use. Default: `NULL` (no additional row removal).

- is_numeric:

  Optional character vector of column names listed in `not_penalized`
  that should be treated as numeric, bypassing the default
  factor/character detection. This mirrors the `is_numeric` argument of
  [`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  and ensures that the NA-filtering step uses the same coercion logic,
  so group sizes are counted on exactly the same sample set. For
  example, if `"Serum albumin (g/L)"` is stored as character but is
  truly numeric, listing it here prevents it from being coerced to
  integer codes - which could introduce spurious NAs or hide real ones.
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

  Logical. If `TRUE` and `pred` is categorical, the result becomes a
  list holding the `group_sizes` data frame and a
  [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  `summary_table` of `metadata` stratified by `pred`. Default: `FALSE`.

- summarytable_args:

  A named list of additional arguments forwarded to
  [`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  when `show_summary_table = TRUE`. Default:
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A report-ready `data.frame` with the following columns (or a list
containing it; see `show_summary_table`):

- `Taxon`:

  Only when `metadata` is a named list.

- `Group`:

  Level of the `pred` variable.

- `Total n`:

  Number of samples in that group (after NA removal).

- `Projected training n`:

  Projected training-set size: `floor(Total n * train_percent)`.

- `Folds requested`:

  The value of `cv_folds`.

- `Fold reduction expected`:

  Logical; `TRUE` when the projected training size is below `cv_folds`.

- `Suggested CV folds`:

  The fold count
  [`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  would actually use: the projected training size, but at least 3.

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
element is processed and results are stacked with a leading `Taxon`
column.

## References

Kuhn M (2008). Building predictive models in R using the caret package.
*Journal of Statistical Software*, 28(5), 1-26.
[doi:10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05)

## See also

[`run_regreg`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md),
[`get_removed_samples`](https://jllcalorio.github.io/pondeR/reference/get_removed_samples.md),
[`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)

## Author

John Lennon L. Calorio

## Examples

``` r
set.seed(1)
meta <- data.frame(
  Group = sample(c("Case", "Control"), 60, TRUE),
  Age   = round(runif(60, 20, 60)),
  Sex   = sample(c("F", "M"), 60, TRUE)
)
meta$Age[c(3, 9)] <- NA   # rows with missing covariates are dropped

get_groupsizes(meta, pred = "Group", not_penalized = c("Age", "Sex"),
               train_percent = 0.7, cv_folds = 5L)
#>     Group Total n Projected training n Folds requested Fold reduction expected
#> 1    Case      30                   21               5                   FALSE
#> 2 Control      28                   19               5                   FALSE
#>   Suggested CV folds
#> 1                 21
#> 2                 19
```
