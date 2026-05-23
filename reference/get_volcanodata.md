# Get Volcano Plot Data

Extracts and combines results from
[`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
and
[`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
into a single data frame. This function automatically aligns
features/outcomes, retrieves fold-change and significance thresholds
from the provided objects, and classifies each feature as "Up", "Down",
or "NS" (Non-significant).

## Usage

``` r
get_volcanodata(..., up = NULL, down = NULL, pval = NULL, filter = TRUE)
```

## Arguments

- ...:

  Objects of class `run_foldchange` and/or `run_diff`.

- up, down, pval:

  Optional numeric overrides for thresholds. If `NULL` (default),
  thresholds are retrieved from the objects' parameters.

- filter:

  Logical. If `TRUE` (default), the data frame is filtered to include
  only "Up" or "Down" regulated features.

## Value

A `data.frame` containing merged data from both objects, including
columns for fold change, p-values, and a `.regulation` classification
column.

## Details

The function automatically detects the thresholds:

- `up` and `down` are retrieved from the `run_foldchange` parameters.

- Significance threshold (`pval`) is retrieved from `run_diff`
  parameters (defaulting to 0.05 if not found).

- It automatically prefers `adj_p_value` over `p_value` if available in
  the `run_diff` summary table.

## See also

[`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md),
[`plot_volcano`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md),
[`run_DIpreprocess`](https://jllcalorio.github.io/pondeR/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming fc_res and diff_res are outputs from pondeR functions
volcano_data <- get_volcanodata(fc_res, diff_res, filter = TRUE)
} # }
```
