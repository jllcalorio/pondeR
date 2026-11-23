# Changelog

## pondeR 1.0.0

First CRAN release.

### Breaking changes

- Plot text now has one control. `base_size` defaults to 11, which suits
  a default 7 x 5 inch Quarto figure, and scales all text.
  `global_font_size` and `zoom` are removed. `title_size`,
  `subtitle_size`, `xlab_size`, `ylab_size`, `axis_text_size`,
  `legend_size`, `label_size` and similar arguments default to `NULL`,
  which follows `base_size`.
- Duplicate or unused arguments are removed:
  `plot_bar(values_position, axis_title_size)`,
  `plot_diff(plot_vars, posthoc)` (use `outcome`),
  `run_auc(compute_auc)` and `plot_dend(leaf_label_size)`.
- Output tables have report-ready column names, such as
  `Adjusted p-value`, `95% CI lower` and `Fold change`, so code that
  reads the old names needs updating.
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  keeps per-feature values under short names in `$feature_stats`.
- [`run_randomize()`](https://jllcalorio.github.io/pondeR/reference/run_randomize.md),
  [`run_subset()`](https://jllcalorio.github.io/pondeR/reference/run_subset.md)
  and
  [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md)
  set a random seed only when `seed` is supplied.
- [`save_excel()`](https://jllcalorio.github.io/pondeR/reference/save_excel.md)
  and
  [`save_plots()`](https://jllcalorio.github.io/pondeR/reference/save_plots.md)
  no longer write to the working directory by default. Supply
  `filename`, `folder`, or both.

### New features

- Functions that take data accept
  [`dimsprepr::run_DIpreprocess()`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html)
  results directly, including custom group columns and merged
  replicates. Results from several preprocessing methods give a clear
  error asking you to pick one.
- [`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
  works without `color_by` and `points_from`, using the fitted groups
  and row names.
- [`run_lefse()`](https://jllcalorio.github.io/pondeR/reference/run_lefse.md)
  gains a [`print()`](https://rdrr.io/r/base/print.html) method.

### Bug fixes

- [`plot_meanmap()`](https://jllcalorio.github.io/pondeR/reference/plot_meanmap.md):
  rows no longer appear mirrored against the dendrogram, and cluster
  labels are no longer clipped.
- [`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md)
  now reads formatted p-values correctly when given a
  [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  summary table.
- [`run_correl()`](https://jllcalorio.github.io/pondeR/reference/run_correl.md):
  automatic selection no longer rejects Pearson correlation because of
  unequal variances, and Spearman and Kendall confidence intervals use
  the Fieller (1957) standard errors.
- [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  always centers the data, and
  [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  and
  [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  now record which groups were excluded.
- [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  no longer draws diagnostic plots while fitting. These wrote
  `Rplots.pdf` to the working directory in non-interactive sessions and
  added unwanted figures to reports.
- `run_diff(calculate_effect_size = FALSE)` is now honoured.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) work for
  multi-outcome
  [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  results.
- [`run_relabel()`](https://jllcalorio.github.io/pondeR/reference/run_relabel.md)
  recodes all labels at once, so chained labels such as `A -> B, B -> C`
  are correct, and several labels can merge into one.
- [`run_reduce()`](https://jllcalorio.github.io/pondeR/reference/run_reduce.md)
  results now have their class, so they print correctly.
- [`run_lefse()`](https://jllcalorio.github.io/pondeR/reference/run_lefse.md)
  handles missing values.
- `get_groupsizes(show_summary_table = TRUE)` no longer fails.
- [`plot_relabund()`](https://jllcalorio.github.io/pondeR/reference/plot_relabund.md)
  honours `ignore_normalization`, no longer prints from inside the
  function, and no longer needs a dimsprepr method that was removed.
- [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md):
  the legend title follows `ci_level`, and labels stay correct when
  confidence intervals are turned off.
- [`plot_multi_y()`](https://jllcalorio.github.io/pondeR/reference/plot_multi_y.md):
  default colours are distinct across all series.
- [`plot_dend()`](https://jllcalorio.github.io/pondeR/reference/plot_dend.md):
  long leaf labels are no longer cut off.
- [`plot_heatmap()`](https://jllcalorio.github.io/pondeR/reference/plot_heatmap.md):
  non-significant crosses no longer hide the coefficients.
- [`run_randomize()`](https://jllcalorio.github.io/pondeR/reference/run_randomize.md):
  the default `per = "none"` samples intact rows. The previous default
  scrambled every cell across the table.
- [`run_lefse()`](https://jllcalorio.github.io/pondeR/reference/run_lefse.md):
  LDA effect sizes follow Segata et al. (2011).
- [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)
  and
  [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  need ‘dimsprepr’ only when `scale_method` or `transform_method` is
  used.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`as_gt()`](https://jllcalorio.github.io/pondeR/reference/as_gt.md)
  for
  [`run_correl()`](https://jllcalorio.github.io/pondeR/reference/run_correl.md)
  results show the correlation table.
- [`as_gt()`](https://jllcalorio.github.io/pondeR/reference/as_gt.md)
  forwards gtsummary tables to
  [`gtsummary::as_gt()`](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html).
- `save_excel(list())` no longer calls a missing internal function.
- [`plot_bar()`](https://jllcalorio.github.io/pondeR/reference/plot_bar.md)
  value labels follow dodged, stacked and filled bars.
- [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md)
  no longer warns about `mark_posthoc_significance` on default calls.
- [`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md)
  no longer warns when `max_pc` exceeds the number of components.

### Package infrastructure

- Examples run without `\dontrun{}`, and slow ones use `\donttest{}`.
- Unit tests cover every exported function.
- References with DOIs are given for the statistical methods.
- Fewer required dependencies: ‘magrittr’, ‘moments’, ‘purrr’,
  ‘reshape2’, ‘scales’, ‘tibble’ and ‘viridis’ are no longer used, and
  ‘lme4’ and ‘lmerTest’ are optional.
