# Internal single-outcome run_diff for plot_diff

A lighter version of `run_diff` that returns a rich list suitable for
downstream use by `plot_diff`. Not intended for direct use.

## Usage

``` r
.run_diff_single(
  x,
  metadata = NULL,
  outcome,
  group = NULL,
  filter = NULL,
  within = NULL,
  subject_id = NULL,
  paired = FALSE,
  test_type = c("auto", "parametric", "nonparametric"),
  normality_method = c("auto", "shapiro", "lilliefors"),
  alpha_normality = 0.05,
  alpha_variance = 0.05,
  alpha_sphericity = 0.05,
  type_sosquares = 2,
  rm_effect_size = "ges",
  correction = "auto",
  test_alpha = 0.05,
  min_n_threshold = 3,
  calculate_effect_size = TRUE,
  perform_posthoc = TRUE,
  p_adjust_method = "BH",
  group_order = NULL,
  subgroup = NULL,
  verbose = FALSE,
  num_cores = 1L
)
```
