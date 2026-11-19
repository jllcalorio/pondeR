# Print a run_correl object

Print a run_correl object

## Usage

``` r
# S3 method for class 'run_correl'
print(x, format_p = TRUE, ...)
```

## Arguments

- x:

  A
  [`run_correl`](https://jllcalorio.github.io/pondeR/reference/run_correl.md)
  object.

- format_p:

  Logical. If `TRUE` (default), the `p_value` and `p_value_adj` columns
  (when present) are displayed via
  [`format.pval`](https://rdrr.io/r/base/format.pval.html) so that
  values indistinguishable from zero at double precision print as
  `"< 2.2e-16"`, matching
  [`stats::cor.test`](https://rdrr.io/r/stats/cor.test.html)'s
  convention. Purely cosmetic: the object returned (invisibly) is
  unchanged, so p-values remain numeric for downstream filtering (e.g.
  `filter(p_value_adj < .05)`).

- ...:

  Additional arguments passed to `print.data.frame`.
