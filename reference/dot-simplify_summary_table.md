# Simplify a run_diff \$summary_table into a compact, presentation-ready table

Called by
[`.build_summary_table`](https://jllcalorio.github.io/pondeR/reference/dot-build_summary_table.md)
when `simplify = TRUE`. Strips verbose `df` labels (meaning moved to a
footnote), formats p-values / numeric columns, compacts effect-size
columns into a single `"ES (95% CI)"` column with symbolic metric
notation, drops `test_used` in favour of superscript markers on
`p_value` / `p_value_adj`, and line-breaks `interpretation` at each
`"; "`. The table's footnotes are attached as the `"df_footnote"` and
`"pvalue_legend"` attributes, rendered by
[`as_gt.run_diff_summary`](https://jllcalorio.github.io/pondeR/reference/as_gt.run_diff_summary.md).

## Usage

``` r
.simplify_summary_table(
  tbl,
  n_digits_pvalues = 3,
  n_digits = 2,
  zero_as_exp = TRUE
)
```

## Arguments

- tbl:

  A `summary_table` data frame as built by
  [`.build_summary_table`](https://jllcalorio.github.io/pondeR/reference/dot-build_summary_table.md)
  (pre-simplification column order).
