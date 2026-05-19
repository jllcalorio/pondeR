# Internal single-outcome run_diff for plot_diff

A lighter version of `run_diff` that returns a rich list suitable for
downstream use by `plot_diff`. Not intended for direct use.

## Usage

``` r
.run_diff_single(
  x,
  outcome,
  group,
  paired = FALSE,
  test_type = c("auto", "parametric", "nonparametric"),
  normality_method = c("auto", "shapiro", "lilliefors"),
  alpha_normality = 0.05,
  alpha_variance = 0.05,
  test_alpha = 0.05,
  min_n_threshold = 3,
  calculate_effect_size = TRUE,
  perform_posthoc = TRUE,
  p_adjust_method = "BH",
  group_order = NULL,
  verbose = FALSE,
  num_cores = 1L
)
```
