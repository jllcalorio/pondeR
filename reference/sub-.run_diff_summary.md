# Subset a run_diff_summary table while keeping its footnote attributes

Base `` `[.data.frame` `` keeps `class` on subsetting but drops other
attributes, so a plain `tbl[rows, ]` (e.g. filtering to significant
rows) silently loses the `df_footnote` / `pvalue_legend` attributes that
[`as_gt.run_diff_summary`](https://jllcalorio.github.io/pondeR/reference/as_gt.run_diff_summary.md)
relies on even though `class(tbl)` still reads `"run_diff_summary"`.
This method re-attaches them after subsetting. It does not help across
[`rbind()`](https://rdrr.io/r/base/cbind.html) / `dplyr` verbs /
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), which
bypass this method entirely; recombine `run_diff_summary` tables *after*
calling
[`as_gt`](https://jllcalorio.github.io/pondeR/reference/as_gt.md) on
each (or re-attach the attributes manually) if you need the footnotes to
survive a merge.

## Usage

``` r
# S3 method for class 'run_diff_summary'
x[...]
```

## Arguments

- x:

  A `run_diff_summary` data frame.

- ...:

  Forwarded to `` `[.data.frame` ``.

## Value

The subset, re-classed as `run_diff_summary` when the result is still a
data frame.
