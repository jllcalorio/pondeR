# Convert a simplified run_diff summary_table to a gt table

Renders the `df_footnote` and `pvalue_legend` attributes of the
simplified table as source notes, mirroring
[`as_gt.run_correl`](https://jllcalorio.github.io/pondeR/reference/as_gt.md)'s
footnote pattern. Row-subsetting via `` `[` `` (see
[`[.run_diff_summary`](https://jllcalorio.github.io/pondeR/reference/sub-.run_diff_summary.md))
preserves these attributes; operations that bypass it
([`rbind()`](https://rdrr.io/r/base/cbind.html), `dplyr` verbs,
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)) do not,
even though `class(x)` still reads `"run_diff_summary"` - in that case a
generic fallback note is shown instead of nothing.

## Usage

``` r
# S3 method for class 'run_diff_summary'
as_gt(x, ...)
```

## Arguments

- x:

  A `run_diff_summary` data frame (i.e. `$summary_table` from
  `run_diff(..., summary_table = TRUE)` with `simplify = TRUE`).

- ...:

  Additional arguments forwarded to
  [`gt::gt()`](https://gt.rstudio.com/reference/gt.html).

## Value

A `gt_tbl` object.
