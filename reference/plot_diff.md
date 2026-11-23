# Create Box or Violin Plots with Statistical Comparisons

Generates box plots or violin plots with statistical tests and pairwise
comparisons. Automatically selects appropriate statistical tests based
on assumptions, or allows manual specification. Statistical testing uses
the same engine as
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md).

## Usage

``` r
plot_diff(
  x,
  outcome = NULL,
  group = NULL,
  metadata = NULL,
  plot_type = "boxplot",
  test_type = NULL,
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
  linewidth = 0.5,
  show_global_p = FALSE,
  global_p_position = NULL,
  global_p_x = NULL,
  base_size = 11,
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

  A data frame containing the variables to plot, OR a pre-computed
  multi-outcome `"run_diff"` object (i.e. the result of
  `run_diff(..., summary_table = TRUE)` with more than one `outcome`).
  When `x` is a `"run_diff"` object, `plot_diff` reuses its statistics
  instead of recomputing them, and prefers the `Adjusted p-value` column
  of `x$summary_table` (falling back to `p-value` when adjustment was
  not performed, i.e. a single outcome) for determining significance.
  `subgroup` and `sample_id`/`label_points` are not supported in this
  mode; run `run_diff(..., subgroup = ...)` beforehand, or pass the
  underlying data frame instead.

- outcome:

  Character vector of numeric outcome column(s) to plot. When `x` is a
  `"run_diff"` object, `NULL` (default) plots every outcome it contains.

- group:

  Character string specifying the column name for the grouping variable.
  Ignored when `x` is a `"run_diff"` object (the grouping variable used
  by
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  is reused automatically).

- metadata:

  A data frame containing metadata (e.g., grouping variables). If
  `NULL`, `x` is used as metadata. Default is `NULL`. Not used when `x`
  is a `"run_diff"` object.

- plot_type:

  Character string specifying plot type: `"boxplot"` or `"violin"`.
  Default is `"boxplot"`.

- test_type:

  Character string or `NULL`. Statistical test strategy: `"auto"`,
  `"parametric"`, or `"nonparametric"`. See
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  for details. If `NULL` (default), uses `"auto"`.

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

  Character vector or `NULL`. Custom plot title(s). If a vector, it
  should match the length of variables being plotted. Default is `NULL`.

- plot_subtitle:

  Character string or `NULL`. Custom subtitle. If `NULL`, the test name
  is shown. Default is `NULL`.

- xlab:

  Character string or `NULL`. X-axis label. Default is `NULL`.

- ylab:

  Character string or `NULL`. Y-axis label. Default is `NULL`.

- linewidth:

  Numeric. The width of the lines for the box/violin outlines and
  statistical brackets. Default is 0.5.

- show_global_p:

  Logical. If `TRUE`, annotates the plot with the omnibus p-value.
  Default is `FALSE`.

- global_p_position:

  Numeric or `NULL`. Y-axis position for the global p-value annotation.
  Default is `NULL` (automatic).

- global_p_x:

  Numeric or `NULL`. X-axis position for the global p-value annotation.
  Default is `NULL` (automatic).

- base_size:

  Base font size in points; every text element scales from it. Default
  is `11`, which suits a default 7 x 5 inch Quarto figure.

- title_size, subtitle_size, xlab_size, ylab_size, axis_text_size,
  pvalue_size:

  Font sizes in points for the title, subtitle, x-axis label, y-axis
  label, tick labels and p-value annotations. Default `NULL` derives
  each from `base_size`.

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

  Named list of test results for each variable, with raw test objects,
  descriptive statistics, and post-hoc tables.

- subgroup_analysis:

  Nested list of results per subgroup variable and level, or `NULL` if
  no subgroup was specified.

## Details

Create Box or Violin Plots with Statistical Comparisons

## References

Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

Mann HB, Whitney DR (1947). On a test of whether one of two random
variables is stochastically larger than the other. *The Annals of
Mathematical Statistics*, 18(1), 50-60.
[doi:10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)

Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance
analysis. *Journal of the American Statistical Association*, 47(260),
583-621.
[doi:10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)

Welch BL (1947). The generalization of 'Student's' problem when several
different population variances are involved. *Biometrika*, 34(1-2),
28-35.
[doi:10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28)

Dunn OJ (1964). Multiple comparisons using rank sums. *Technometrics*,
6(3), 241-252.
[doi:10.1080/00401706.1964.10490181](https://doi.org/10.1080/00401706.1964.10490181)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# Example 1: Box plot with mean dot and mean value annotation
results <- plot_diff(
  x            = PlantGrowth,
  outcome      = "weight",
  group        = "group",
  show_mean    = TRUE,
  show_agg_val = TRUE
)
print(results$plots$weight)


# Example 2: Larger text and omnibus p-value annotation
results2 <- plot_diff(
  x             = PlantGrowth,
  outcome       = "weight",
  group         = "group",
  base_size     = 14,
  show_global_p = TRUE
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

```
