# Fold Change vs. Feature Metric Plot (MA-style)

Plots the fold change of each feature (y-axis) against a per-feature
metric (x-axis): mean abundance (the classic MA plot), effect size,
variance, standard deviation, or VIP score. Unlike
[`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)
it needs no p-values, so it remains useful when
[`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
found nothing significant. With `metric = "abundance"` it shows whether
large fold changes come from low-intensity "noise" features whose
baseline is near zero. Styling arguments match
[`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md).

## Usage

``` r
plot_fcmetric(
  x,
  y = NULL,
  z = NULL,
  group = NULL,
  features = NULL,
  metric = c("abundance", "effect_size", "variance", "sd", "vip"),
  log2 = TRUE,
  log10_x = NULL,
  x_cut = NULL,
  up = 1.5,
  down = 0.5,
  annotate = NULL,
  use_ggrepel = TRUE,
  up_color = "red",
  down_color = "blue",
  ns_color = "grey60",
  h_color = NULL,
  v_color = NULL,
  size = 2,
  alpha = 0.6,
  theme = "nature",
  plot_title = NULL,
  plot_subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  base_size = 11,
  title_size = NULL,
  subtitle_size = NULL,
  xlab_size = NULL,
  ylab_size = NULL,
  axis_text_size = NULL,
  label_size = NULL,
  legend_size = NULL,
  legend_title = NULL,
  legend_nrow = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  result, or a named list holding it plus the object needed by `metric`,
  e.g. `list(fc = fc, diff = diff, pls = pls)`; **or** a
  `data.frame`/`tibble`/named `matrix` with one row per feature. One
  plot is made per comparison.

  The data stage (normalized, transformed, scaled, ...) is chosen
  upstream, e.g. `run_foldchange(from = "normalized")`. For a
  `data.frame`, supply columns computed from whichever stage you want.

- y:

  A single character string naming the column of **raw** fold changes
  (positive numeric). Default `NULL` (`"fold_change"` for
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  input).

- z:

  A single character string naming the column holding the x-axis metric.
  Default `NULL` (derived from `x` and `metric` for
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  input; required for a `data.frame`, where `metric` then only sets the
  axis label and default scale).

- group:

  An optional column name in `x`; one plot is made per level. Default
  `NULL`.

- features:

  A single character string naming the feature identifier column.
  Required when `annotate` is used. Default `NULL`.

- metric:

  A single string, the x-axis metric. See *Metrics*. Default
  `"abundance"`.

- log2:

  Logical. If `TRUE` (default) the y-axis is \\\log_2\\(fold change); if
  `FALSE` it is the raw fold change.

- log10_x:

  Logical or `NULL`. If `TRUE` the x-axis is \\\log\_{10}\\ of the
  metric. Default `NULL` chooses automatically: `TRUE` for
  `"abundance"`, `"variance"` and `"sd"`; `FALSE` for `"effect_size"`
  and `"vip"`. If the metric has non-positive values it is set to
  `FALSE` with a warning.

- x_cut:

  Optional numeric vector of reference values (on the original metric
  scale, not log10) drawn as vertical dashed lines, e.g. `1` for VIP or
  `c(0.2, 0.5, 0.8)` for Cohen's d. Default `NULL`.

- up, down:

  Raw fold-change thresholds (`up > 1`, `0 < down < 1`). Features beyond
  them are coloured Up/Down. Inherited from a
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  object when missing; otherwise `1.5`/`0.5`.

- annotate:

  Optional (named) character vector of feature identifiers to label; see
  [`plot_volcano`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md).
  Default `NULL`.

- use_ggrepel, up_color, down_color, ns_color, h_color, v_color, size,
  alpha, theme, plot_title, plot_subtitle, xlab, ylab, base_size,
  title_size, subtitle_size, xlab_size, ylab_size, axis_text_size,
  label_size, legend_size, legend_title, legend_nrow, ...:

  As in
  [`plot_volcano`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md).
  Here `h_color` colours the horizontal `up`/`down` fold-change lines
  and `v_color` the vertical `x_cut` lines. `plot_title = NULL` gives a
  title naming the metric (plus the comparison when there are several).

## Value

A list of class `plot_fcmetric` with `plots` (named list of `ggplot`
objects), `classified` (a table with `Comparison`, `Feature`,
`Fold change`, `log2 FC`, the metric and `Regulation`) and `params`.

## Metrics

- `"abundance"`:

  Mean of the two group means of each comparison. Needs only
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md).

- `"variance"`, `"sd"`:

  Pooled within-group variance / SD of the two groups in each
  comparison, \\s_p^2 = ((n_1-1)s_1^2 + (n_2-1)s_2^2)/(n_1+n_2-2)\\, on
  the same data scale as the fold change. Needs only
  [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md).
  Raw-intensity variance typically grows with abundance; consider
  transformed data for a variance-stabilised view.

- `"effect_size"`:

  Per-feature estimate from
  [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  (multi-outcome; e.g. Cohen's d, rank-biserial, eta-squared depending
  on the test chosen). Supply `list(fc = ..., diff = ...)`. The estimate
  is one value per feature, so with more than two groups it is omnibus
  and repeats across comparisons. Estimates from different effect-size
  metrics are not directly comparable; a warning is given when they are
  mixed.

- `"vip"`:

  VIP score from
  [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  (PLS-DA/OPLS-DA; sPLS-DA yields none). Supply
  `list(fc = ..., pls = ...)`. One value per feature, repeated across
  comparisons.

## References

Cohen J (1988). *Statistical Power Analysis for the Behavioral Sciences*
(2nd ed.). Lawrence Erlbaum Associates.

Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York.
[doi:10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)

## See also

[`plot_volcano`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md),
[`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md),
[`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)

## Author

John Lennon L. Calorio

## Examples

``` r
x    <- iris[iris$Species != "setosa", 1:4]
meta <- droplevels(iris[iris$Species != "setosa", "Species", drop = FALSE])
fc   <- run_foldchange(x, metadata = meta, group = "Species")

plot_fcmetric(fc)                                  # MA-style plot
#> -- plot_fcmetric results -----------------------------------
#>   Metric     : abundance
#>   Features   : 4  (Up: 0 | Down: 0 | NS: 4)
#>   Thresholds : FC >= 1.5000 (up) | <= 0.5000 (down)
#>   Plots      : Fold Change vs. Mean Abundance
plot_fcmetric(fc, metric = "sd", base_size = 12)
#> -- plot_fcmetric results -----------------------------------
#>   Metric     : sd
#>   Features   : 4  (Up: 0 | Down: 0 | NS: 4)
#>   Thresholds : FC >= 1.5000 (up) | <= 0.5000 (down)
#>   Plots      : Fold Change vs. Standard Deviation

# \donttest{
dif <- run_diff(x, meta, outcome = names(x), group = "Species",
                summary_table = TRUE, verbose = FALSE)
plot_fcmetric(list(fc = fc, diff = dif), metric = "effect_size",
              x_cut = c(0.2, 0.5, 0.8))
#> Warning: Effect sizes use different metrics (Rank-biserial correlation, Cohen's d); they are not directly comparable.
#> -- plot_fcmetric results -----------------------------------
#>   Metric     : effect_size
#>   Features   : 4  (Up: 0 | Down: 0 | NS: 4)
#>   Thresholds : FC >= 1.5000 (up) | <= 0.5000 (down)
#>   Plots      : Fold Change vs. Effect Size
# }
```
