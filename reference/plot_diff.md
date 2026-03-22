# Create Box or Violin Plots with Statistical Comparisons

This function generates box plots or violin plots with statistical tests
and pairwise comparisons. It automatically selects appropriate
statistical tests based on assumptions, or allows manual specification.
Statistical testing is performed using the companion function
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md).

## Usage

``` r
plot_diff(
  x,
  outcome,
  group,
  plot_vars = NULL,
  plot_type = "boxplot",
  test_type = NULL,
  posthoc = NULL,
  p_adjust_method = "BH",
  hide_ns = TRUE,
  show_p_numeric = FALSE,
  paired = FALSE,
  test_alpha = 0.05,
  colors = NULL,
  plot_title = NULL,
  plot_subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  show_global_p = FALSE,
  global_p_position = NULL,
  global_p_x = NULL,
  global_font_size = NULL,
  title_size = NULL,
  subtitle_size = NULL,
  xlab_size = NULL,
  ylab_size = NULL,
  pvalue_size = NULL,
  axis_text_size = NULL,
  format_p_numeric = "auto",
  group_order = NULL,
  label_points = NULL,
  sample_id = NULL,
  show_mean = FALSE,
  show_agg_val = TRUE,
  horizontal = FALSE,
  bracket_spacing = 0.08,
  subgroup = NULL,
  ...
)
```

## Arguments

- x:

  A data frame containing the variables to plot.

- outcome:

  Character string specifying the column name for the numeric outcome
  variable. Used as the default for `plot_vars` if it is not specified.

- group:

  Character string specifying the column name for the grouping variable.

- plot_vars:

  Character vector specifying column name(s) for numeric variables to
  plot. If `NULL` (default), uses the value of `outcome`.

- plot_type:

  Character string specifying plot type: `"boxplot"` or `"violin"`.
  Default is `"boxplot"`.

- test_type:

  Character string or NULL. Statistical test strategy: `"auto"`,
  `"parametric"`, or `"nonparametric"`. See
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  for details. If NULL (default), uses `"auto"`.

- posthoc:

  Character string or NULL. Post-hoc test for multiple groups: `"auto"`,
  `"none"`, `"tukey"`, `"games-howell"`, `"dunn"`, `"t_test"`, or
  `"wilcox_test"`. If NULL (default), uses `"auto"`.

- p_adjust_method:

  Character string. P-value adjustment method for multiple comparisons.
  Options: `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`,
  `"BY"`, `"fdr"`, `"none"`. Default is `"BH"`.

- hide_ns:

  Logical. If TRUE, hide non-significant p-values from plot brackets.
  Default is TRUE.

- show_p_numeric:

  Logical. If TRUE, show exact p-values as numbers instead of
  significance symbols (`*`, `**`, `***`). Default is FALSE.

- paired:

  Logical. If TRUE, use paired tests. Default is FALSE.

- test_alpha:

  Numeric. Significance level for determining significance on the plot.
  Default is 0.05.

- colors:

  Character vector of colors (names or HEX codes) or NULL for automatic
  soft colors. Default is NULL.

- plot_title:

  Character string or NULL. Custom title for the plot. If NULL,
  generates an automatic title. Default is NULL.

- plot_subtitle:

  Character string or NULL. Custom subtitle for the plot. If NULL, the
  name of the statistical test used is shown automatically. Default is
  NULL.

- xlab:

  Character string or NULL. X-axis label. If NULL, uses the group
  variable name. Default is NULL.

- ylab:

  Character string or NULL. Y-axis label. If NULL, uses the plot
  variable name. Default is NULL.

- show_global_p:

  Logical. If TRUE, annotates the plot with the global (omnibus) p-value
  from the main statistical test (e.g., Kruskal-Wallis, ANOVA). Shown
  for all group counts when TRUE. Default is FALSE.

- global_p_position:

  Numeric or NULL. Y-axis position (data coordinates) for the global
  p-value annotation. If NULL (default), positioned automatically near
  the top of the plot area.

- global_p_x:

  Numeric or NULL. X-axis position for the global p-value annotation. If
  NULL (default), positioned at 95% of the x-axis range (right side).

- global_font_size:

  Numeric or NULL. Base font size used to derive all text sizes
  proportionally. When provided, `title_size`, `subtitle_size`,
  `xlab_size`, `ylab_size`, `pvalue_size`, `axis_text_size`, and
  annotation text sizes (mean/median values, global p) are all scaled
  relative to this value. Individual size parameters override
  `global_font_size` only for their specific element when set to a
  non-NULL numeric. Default is NULL (each parameter uses its own
  explicit default).

- title_size:

  Numeric or NULL. Font size for the plot title. If NULL and
  `global_font_size` is set, derived proportionally (1.2x). Default is
  NULL.

- subtitle_size:

  Numeric or NULL. Font size for the plot subtitle. If NULL and
  `global_font_size` is set, derived proportionally (0.9x). Default is
  NULL.

- xlab_size:

  Numeric or NULL. Font size for the x-axis label. If NULL and
  `global_font_size` is set, derived proportionally (1.0x). Default is
  NULL.

- ylab_size:

  Numeric or NULL. Font size for the y-axis label. If NULL and
  `global_font_size` is set, derived proportionally (1.0x). Default is
  NULL.

- pvalue_size:

  Numeric or NULL. Font size (in ggplot2 `size` units, i.e. mm) for
  p-value annotations on brackets and the global p-value. If NULL and
  `global_font_size` is set, derived proportionally (0.35x). Default is
  NULL.

- axis_text_size:

  Numeric or NULL. Font size for axis tick labels (both x and y axes).
  If NULL and `global_font_size` is set, derived proportionally (0.85x).
  Default is NULL.

- format_p_numeric:

  Character string or numeric. Controls how p-values are formatted when
  `show_p_numeric = TRUE`. `"auto"` (default) rounds to 3 decimal
  places, showing `"<.001"` if below that threshold and `">.999"` if
  above. A numeric value specifies the exact number of decimal places to
  display.

- group_order:

  Character vector or NULL. Specifies the order of groups on the x-axis.
  If NULL (default), uses factor levels or alphabetical order.

- label_points:

  Character string or NULL. Point labeling: `"all"` (label all jitter
  points), `"outliers"` (label only outliers), or NULL (no labels).
  Requires `sample_id`. Default is NULL.

- sample_id:

  Character string or NULL. Column name containing sample identifiers
  for point labeling. Required when `label_points` is not NULL.

- show_mean:

  Logical. If TRUE, adds a red dot (mean indicator) to each group.
  Default is FALSE.

- show_agg_val:

  Logical. If TRUE (default), shows the numeric mean (parametric) or
  median (non-parametric) value as text on the plot.

- horizontal:

  Logical. If TRUE, flips the plot to horizontal orientation. Default is
  FALSE.

- bracket_spacing:

  Numeric. Multiplier controlling the vertical spacing between
  comparison brackets. Default is 0.08.

- subgroup:

  Character vector or NULL. Column name(s) in `x` for subgroup analysis.
  Each must be categorical (factor or character). When specified, the
  full analysis is repeated for each level of each subgroup column.
  Default is NULL.

- ...:

  Additional arguments passed to `ggboxplot` or `ggviolin` from ggpubr.

## Value

A list containing:

- plots:

  List of all generated plots.

- significant_plots:

  List of plots with significant results only.

- statistics:

  List of
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  results for each variable.

- subgroup_analysis:

  Nested list of results per subgroup variable and level, or NULL if no
  subgroup was specified.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Show mean dot AND mean value
results_anova <- plot_diff(
  x            = PlantGrowth,
  outcome      = "weight",
  group        = "group",
  plot_type    = "boxplot",
  show_mean    = TRUE,
  show_agg_val = TRUE
)
print(results_anova$plots$weight)

# Example 2: Global font size controls all text proportionally
results_anova2 <- plot_diff(
  x               = PlantGrowth,
  outcome         = "weight",
  group           = "group",
  global_font_size = 16,
  show_global_p   = TRUE  # show the Kruskal-Wallis / ANOVA p-value on the plot
)
print(results_anova2$plots$weight)

# Example 3: Paired t-test with custom subtitle, formatted p-values
sleep_data     <- sleep
sleep_data$ID  <- rep(1:10, 2)
results_paired <- plot_diff(
  x              = sleep_data,
  outcome        = "extra",
  group          = "group",
  paired         = TRUE,
  plot_subtitle  = "Paired comparison: Drug 1 vs Drug 2",
  show_p_numeric = TRUE,
  format_p_numeric = 4  # 4 decimal places
)
print(results_paired$plots$extra)
} # }
```
