# Create Box or Violin Plots with Statistical Comparisons

Generates box plots or violin plots with statistical tests and pairwise
comparisons. Automatically selects appropriate statistical tests based
on assumptions, or allows manual specification. Statistical testing is
performed internally via `.run_diff_single`, the unexported
single-outcome engine that powers
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

  Character string or `NULL`. Statistical test strategy: `"auto"`,
  `"parametric"`, or `"nonparametric"`. See
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  for details. If `NULL` (default), uses `"auto"`.

- posthoc:

  Character string or `NULL`. Retained for API compatibility; post-hoc
  selection is handled automatically by the internal engine. Default is
  `NULL`.

- p_adjust_method:

  Character string. P-value adjustment method for post-hoc comparisons:
  `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`,
  `"fdr"`, `"none"`. Default is `"BH"`.

- hide_ns:

  Logical. If `TRUE`, hide non-significant p-values from plot brackets.
  Default is `TRUE`.

- show_p_numeric:

  Logical. If `TRUE`, show exact p-values instead of significance
  symbols. Default is `FALSE`.

- paired:

  Logical. If `TRUE`, use paired tests. Default is `FALSE`.

- test_alpha:

  Numeric. Significance level for determining significance on the plot.
  Default is 0.05.

- colors:

  Character vector of colors or `NULL` for automatic soft colors.
  Default is `NULL`.

- plot_title:

  Character string or `NULL`. Custom plot title. Default is `NULL`.

- plot_subtitle:

  Character string or `NULL`. Custom subtitle. If `NULL`, the test name
  is shown. Default is `NULL`.

- xlab:

  Character string or `NULL`. X-axis label. Default is `NULL`.

- ylab:

  Character string or `NULL`. Y-axis label. Default is `NULL`.

- show_global_p:

  Logical. If `TRUE`, annotates the plot with the omnibus p-value.
  Default is `FALSE`.

- global_p_position:

  Numeric or `NULL`. Y-axis position for the global p-value annotation.
  Default is `NULL` (automatic).

- global_p_x:

  Numeric or `NULL`. X-axis position for the global p-value annotation.
  Default is `NULL` (automatic).

- global_font_size:

  Numeric or `NULL`. Base font size for proportional scaling of all text
  elements. Default is `NULL`.

- title_size:

  Numeric or `NULL`. Title font size. Default is `NULL`.

- subtitle_size:

  Numeric or `NULL`. Subtitle font size. Default is `NULL`.

- xlab_size:

  Numeric or `NULL`. X-axis label font size. Default is `NULL`.

- ylab_size:

  Numeric or `NULL`. Y-axis label font size. Default is `NULL`.

- pvalue_size:

  Numeric or `NULL`. Font size for p-value annotations (ggplot2 `size`
  units). Default is `NULL`.

- axis_text_size:

  Numeric or `NULL`. Axis tick label font size. Default is `NULL`.

- format_p_numeric:

  Character string or numeric. Controls p-value formatting when
  `show_p_numeric = TRUE`. `"auto"` (default) uses 3 decimal places.

- group_order:

  Character vector or `NULL`. Group order on x-axis. Default is `NULL`.

- label_points:

  Character string or `NULL`. `"all"` labels all jitter points;
  `"outliers"` labels only outliers. Requires `sample_id`. Default is
  `NULL`.

- sample_id:

  Character string or `NULL`. Column name for sample identifiers used in
  point labeling. Default is `NULL`.

- show_mean:

  Logical. If `TRUE`, adds a red dot (mean indicator) to each group.
  Default is `FALSE`.

- show_agg_val:

  Logical. If `TRUE` (default), shows the numeric mean (parametric) or
  median (non-parametric) value as text on the plot.

- horizontal:

  Logical. If `TRUE`, flips to horizontal orientation. Default is
  `FALSE`.

- bracket_spacing:

  Numeric. Multiplier for vertical spacing between comparison brackets.
  Default is 0.08.

- subgroup:

  Character vector or `NULL`. Column name(s) for subgroup analysis.
  Default is `NULL`.

- ...:

  Additional arguments passed to `ggboxplot` or `ggviolin` from ggpubr.

## Value

A list containing:

- plots:

  Named list of all generated `ggplot` objects.

- significant_plots:

  Named list of plots for variables with significant results only.

- statistics:

  Named list of internal `.run_diff_single` result objects for each
  variable, providing access to raw test objects, descriptive
  statistics, and post-hoc tables.

- subgroup_analysis:

  Nested list of results per subgroup variable and level, or `NULL` if
  no subgroup was specified.

## Details

Create Box or Violin Plots with Statistical Comparisons

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Box plot with mean dot and mean value annotation
results <- plot_diff(
  x            = PlantGrowth,
  outcome      = "weight",
  group        = "group",
  show_mean    = TRUE,
  show_agg_val = TRUE
)
print(results$plots$weight)

# Example 2: Global font size scaling and omnibus p-value annotation
results2 <- plot_diff(
  x                = PlantGrowth,
  outcome          = "weight",
  group            = "group",
  global_font_size = 16,
  show_global_p    = TRUE
)
print(results2$plots$weight)

# Example 3: Paired comparison with numeric p-values
sleep_data    <- sleep
sleep_data$ID <- rep(1:10, 2)
results3 <- plot_diff(
  x                = sleep_data,
  outcome          = "extra",
  group            = "group",
  paired           = TRUE,
  show_p_numeric   = TRUE,
  format_p_numeric = 4
)
print(results3$plots$extra)
} # }
```
