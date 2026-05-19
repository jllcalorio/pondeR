# Build a summary_table data frame from a list of run_diff objects

Used by the multi-outcome dispatch in `run_diff` to construct the
`$summary_table` element and subgroup-level summary tables.

## Usage

``` r
.build_summary_table(results_list, outcome_names, test_alpha)
```

## Arguments

- results_list:

  Named list of `"run_diff"` objects (one per outcome).

- outcome_names:

  Character vector of outcome names (defines row order).

- test_alpha:

  Significance threshold for the `significant` column.

## Value

A data frame with one row per outcome.
