# Execute RM or mixed ANOVA and return a full run_diff-compatible result list

Called exclusively from the single-outcome path of `run_diff` when
`within` is non-NULL.

## Usage

``` r
.run_rm_anova(
  x,
  outcome,
  group,
  within,
  subject_id,
  test_family,
  normality_method,
  alpha_normality,
  alpha_sphericity,
  type_sosquares,
  rm_effect_size,
  correction,
  test_alpha,
  perform_posthoc,
  p_adjust_method,
  verbose,
  warnings_list,
  parameters
)
```

## Value

A list of class `"run_diff"`.
