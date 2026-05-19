# Compute the appropriate effect size for a classical run_diff result

Compute the appropriate effect size for a classical run_diff result

## Usage

``` r
.compute_effect_size(
  y,
  grp,
  groups,
  test_family,
  parametric_used,
  variance_violated,
  fitted_model,
  warnings_list
)
```

## Value

A list with `estimate`, `ci_low`, `ci_high`, `magnitude`, `metric`, and
`interpretation`.
