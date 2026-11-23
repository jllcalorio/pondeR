# pondeR 1.0.0

First CRAN release.

## Breaking changes

* Plot text now has one control. `base_size` defaults to 11, which suits a
  default 7 x 5 inch Quarto figure, and scales all text. `global_font_size` and
  `zoom` are removed. `title_size`, `subtitle_size`, `xlab_size`, `ylab_size`,
  `axis_text_size`, `legend_size`, `label_size` and similar arguments default
  to `NULL`, which follows `base_size`.
* Duplicate or unused arguments are removed: `plot_bar(values_position,
  axis_title_size)`, `plot_diff(plot_vars, posthoc)` (use `outcome`),
  `run_auc(compute_auc)` and `plot_dend(leaf_label_size)`.
* Output tables have report-ready column names, such as `Adjusted p-value`,
  `95% CI lower` and `Fold change`, so code that reads the old names needs
  updating. `run_foldchange()` keeps per-feature values under short names in
  `$feature_stats`.
* `run_randomize()`, `run_subset()` and `run_regreg()` set a random seed only
  when `seed` is supplied.
* `save_excel()` and `save_plots()` no longer write to the working directory by
  default. Supply `filename`, `folder`, or both.

## New features

* Functions that take data accept `dimsprepr::run_DIpreprocess()` results
  directly, including custom group columns and merged replicates. Results from
  several preprocessing methods give a clear error asking you to pick one.
* `plot_score()` works without `color_by` and `points_from`, using the fitted
  groups and row names.
* `run_lefse()` gains a `print()` method.

## Bug fixes

* `plot_meanmap()`: rows no longer appear mirrored against the dendrogram, and
  cluster labels are no longer clipped.
* `plot_diff()` now reads formatted p-values correctly when given a `run_diff()`
  summary table.
* `run_correl()`: automatic selection no longer rejects Pearson correlation
  because of unequal variances, and Spearman and Kendall confidence intervals
  use the Fieller (1957) standard errors.
* `run_pca()` always centers the data, and `run_pca()` and `run_pls()` now
  record which groups were excluded.
* `run_pls()` no longer draws diagnostic plots while fitting. These wrote
  `Rplots.pdf` to the working directory in non-interactive sessions and added
  unwanted figures to reports.
* `run_diff(calculate_effect_size = FALSE)` is now honoured.
* `print()` and `summary()` work for multi-outcome `run_diff()` results.
* `run_relabel()` recodes all labels at once, so chained labels such as
  `A -> B, B -> C` are correct, and several labels can merge into one.
* `run_reduce()` results now have their class, so they print correctly.
* `run_lefse()` handles missing values.
* `get_groupsizes(show_summary_table = TRUE)` no longer fails.
* `plot_relabund()` honours `ignore_normalization`, no longer prints from inside
  the function, and no longer needs a dimsprepr method that was removed.
* `run_auc()`: the legend title follows `ci_level`, and labels stay correct when
  confidence intervals are turned off.
* `plot_multi_y()`: default colours are distinct across all series.
* `plot_dend()`: long leaf labels are no longer cut off.
* `plot_heatmap()`: non-significant crosses no longer hide the coefficients.
* `run_randomize()`: the default `per = "none"` samples intact rows. The previous
  default scrambled every cell across the table.
* `run_lefse()`: LDA effect sizes follow Segata et al. (2011).
* `run_pca()` and `run_pls()` need 'dimsprepr' only when `scale_method` or
  `transform_method` is used.
* `print()` and `as_gt()` for `run_correl()` results show the correlation table.
* `as_gt()` forwards gtsummary tables to `gtsummary::as_gt()`.
* `save_excel(list())` no longer calls a missing internal function.
* `plot_bar()` value labels follow dodged, stacked and filled bars.
* `run_summarytable()` no longer warns about `mark_posthoc_significance` on
  default calls.
* `plot_scree()` no longer warns when `max_pc` exceeds the number of components.

## Package infrastructure

* Examples run without `\dontrun{}`, and slow ones use `\donttest{}`.
* Unit tests cover every exported function.
* References with DOIs are given for the statistical methods.
* Fewer required dependencies: 'magrittr', 'moments', 'purrr', 'reshape2',
  'scales', 'tibble' and 'viridis' are no longer used, and 'lme4' and 'lmerTest'
  are optional.
