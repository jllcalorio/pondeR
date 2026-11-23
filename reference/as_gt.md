# Convert a pondeR Result to a gt Table

Renders pondeR tables as
[`gt`](https://gt.rstudio.com/reference/gt.html) tables. Any footnote or
p-value legend stored as an attribute (`"footnote"`, `"df_footnote"`,
`"pvalue_legend"`) is appended as a source note via
[`tab_source_note`](https://gt.rstudio.com/reference/tab_source_note.html).

## Usage

``` r
as_gt(x, ...)

# S3 method for class 'gtsummary'
as_gt(x, ...)

# S3 method for class 'data.frame'
as_gt(x, ...)

# Default S3 method
as_gt(x, ...)

# S3 method for class 'run_correl'
as_gt(x, format_p = TRUE, ...)
```

## Arguments

- x:

  An object to convert: a
  [`run_correl`](https://jllcalorio.github.io/pondeR/reference/run_correl.md)
  result, a `run_diff_summary` table from
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md),
  a `gtsummary` table, a data frame/tibble (or anything coercible with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)), or a
  `gt_tbl` (returned unchanged).

- ...:

  Additional arguments forwarded to
  [`gt`](https://gt.rstudio.com/reference/gt.html) (or to
  [`as_gt`](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html)
  for gtsummary objects).

- format_p:

  Logical. If `TRUE` (default), the p-value columns (when present) are
  displayed via [`format.pval`](https://rdrr.io/r/base/format.pval.html)
  so that values indistinguishable from zero at double precision render
  as `"< 2.2e-16"`, matching
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html)'s convention.
  Cosmetic only; does not affect the source `x`.

## Value

A `gt_tbl` object.

## Details

pondeR's generic masks
[`as_gt`](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html)
when both packages are attached, so `gtsummary` objects are forwarded
unchanged to
[`gtsummary::as_gt()`](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html).
A `gt_tbl` is returned as is; any other object is coerced with
[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html).

## References

Iannone R, Cheng J, Schloerke B, et al. (2024). *gt: Easily Create
Presentation-Ready Display Tables*. R package.
<https://CRAN.R-project.org/package=gt>

## See also

[`run_correl`](https://jllcalorio.github.io/pondeR/reference/run_correl.md),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md),
[`run_summarytable`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)

## Examples

``` r
if (requireNamespace("gt", quietly = TRUE)) {
  as_gt(head(mtcars))
  as_gt(run_correl(mtcars[, 1:4], y = "mpg"))
}


  

Variable 1
```
