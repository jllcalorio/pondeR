# Build a summary_table data frame from a list of run_diff objects

Used by the multi-outcome dispatch in
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
to construct the `$summary_table` element and subgroup-level summary
tables.

## Usage

``` r
.build_summary_table(
  results_list,
  outcome_names,
  test_alpha,
  p_adjust_method = "BH",
  simplify = TRUE,
  n_digits_pvalues = 3,
  n_digits = 2,
  zero_as_exp = TRUE
)
```

## Arguments

- results_list:

  Named list of `"run_diff"` objects (one per outcome).

- outcome_names:

  Character vector of outcome names (defines row order).

- test_alpha:

  Significance threshold for the `significant` column.

- simplify:

  Logical. When `TRUE` (default), the table is post-processed into a
  compact, presentation-ready form via
  [`.simplify_summary_table`](https://jllcalorio.github.io/pondeR/reference/dot-simplify_summary_table.md):
  `df` drops its `name = ` labels (meanings moved to a footnote),
  `p_value` / `p_value_adj` are formatted per `n_digits_pvalues` and
  carry test-identity superscripts (in place of the dropped `test_used`
  column), numeric columns are formatted per `n_digits` / `zero_as_exp`,
  effect-size columns are compacted into a single `"ES (95% CI)"` column
  with symbolic metric notation, and `interpretation` line-breaks at
  each `"; "`.

- n_digits_pvalues:

  Passed to
  [`.simplify_summary_table`](https://jllcalorio.github.io/pondeR/reference/dot-simplify_summary_table.md).

- n_digits:

  Passed to
  [`.simplify_summary_table`](https://jllcalorio.github.io/pondeR/reference/dot-simplify_summary_table.md).

- zero_as_exp:

  Passed to
  [`.simplify_summary_table`](https://jllcalorio.github.io/pondeR/reference/dot-simplify_summary_table.md).

## Value

A data frame with one row per outcome.
