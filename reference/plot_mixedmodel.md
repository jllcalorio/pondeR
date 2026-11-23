# Plot Results from a Linear Mixed-Effects Model Analysis

Produces a volcano plot from the `combined_fixed_effects` component of a
[`run_mixedmodel`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md)
object, with log2 fold-change (or raw estimates) on the x-axis and
\\-\log\_{10}(p\_{\text{adj}})\\ on the y-axis. When the model formula
contained more than one fixed effect term, panels are automatically
faceted by `term` so each predictor is shown separately. Points are
colored by significance and direction, and labels can be added to the
top-\\n\\ features per panel or to a user-specified set of features.

## Usage

``` r
plot_mixedmodel(
  x,
  type = "volcano",
  term = NULL,
  use_padj = TRUE,
  p_threshold = 0.05,
  estimate_threshold = 0,
  label_features = NULL,
  minmax = NULL,
  iqr_k = 1.5,
  label_size = NULL,
  label_max_overlap = 20L,
  label_box = FALSE,
  label_color = "match",
  point_size = 2,
  point_size_sig = 2.5,
  alpha = 0.6,
  label_sig_only = TRUE,
  color_up = "#D85A30",
  color_down = "#0072B2",
  color_ns = "#888780",
  x_label = "Estimate (fixed effect)",
  y_label = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  facet_scales = "fixed",
  facet_ncol = NULL,
  theme = "nature",
  base_size = 11,
  font_family = ""
)
```

## Arguments

- x:

  A
  [`run_mixedmodel`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md)
  object returned by
  [`run_mixedmodel`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md).

- type:

  Character string. Currently `"volcano"` is supported. Reserved for
  future types (`"effect"`, `"heatmap"`). Default: `"volcano"`.

- term:

  A character vector of fixed-effect term name(s) to include. When
  `NULL` (default), all non-intercept terms are included. Useful for
  focusing on a single predictor when the formula has many.

- use_padj:

  Logical. If `TRUE` (default), the y-axis reflects `p_adj`
  (Benjamini-Hochberg or whichever method was passed to
  [`run_mixedmodel`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md)).
  If `FALSE`, raw `p_value` is used.

- p_threshold:

  Numeric. Horizontal dashed significance line drawn at
  \\-\log\_{10}(\text{p\\threshold})\\. Default: `0.05`.

- estimate_threshold:

  Numeric. Vertical dashed lines drawn at \\\pm\\`estimate_threshold`.
  Set to `0` to suppress. Default: `0` (no vertical lines).

- label_features:

  A character vector of specific feature names to label, regardless of
  significance rank. Combined with `label_top` (union). Default: `NULL`.

- minmax:

  A numeric vector of length 2, e.g. `c(-0.02, 0.02)`. Features with
  estimate \\\le\\ `minmax[1]` or \\\ge\\ `minmax[2]` are annotated.
  Overrides automatic IQR-based annotation. Default: `NULL`.

- iqr_k:

  Numeric. IQR fence multiplier used for automatic annotation when both
  `minmax` and `label_features` are `NULL`. Features beyond \\Q1 - k
  \times IQR\\ or \\Q3 + k \times IQR\\ are labelled. Set to `Inf` to
  disable automatic annotation entirely. Default: `1.5`.

- label_size:

  Numeric. Font size for point labels, in points. Default: `NULL`,
  derived from `base_size`.

- label_max_overlap:

  Integer. Passed to
  [`ggrepel::geom_text_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)'s
  `max.overlaps`. Default: `20`.

- label_box:

  Logical. If `TRUE`, draws a border box around labels via
  [`ggrepel::geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Default: `FALSE`.

- label_color:

  Character. Color of label text. `"match"` inherits the point's
  significance color. Any valid color string overrides. Default:
  `"match"`.

- point_size:

  Numeric. Base point size for non-significant features. Default: `2`.

- point_size_sig:

  Numeric. Point size for significant features. Default: `2.5`.

- alpha:

  Numeric in (0, 1\]. Transparency for non-significant points. Default:
  `0.6`.

- label_sig_only:

  Logical. If `TRUE` (default), only features passing `p_threshold` are
  eligible for `minmax` or IQR-based auto-annotation. Features in
  `label_features` are always annotated regardless of this setting.
  Default: `TRUE`.

- color_up:

  Character. Color for significantly up-regulated (positive estimate +
  significant) features. Default: `"#D85A30"` (Okabe-Ito vermillion).

- color_down:

  Character. Color for significantly down-regulated (negative estimate +
  significant) features. Default: `"#0072B2"` (Okabe-Ito blue).

- color_ns:

  Character. Color for non-significant features. Default: `"#888780"`.

- x_label:

  Character. x-axis label. Default: `"Estimate (fixed effect)"`.

- y_label:

  Character. y-axis label. Default: `expression(-log[10](p[adj]))` or
  `expression(-log[10](p))` depending on `use_padj`.

- title:

  Character. Plot title. Default: `NULL` (no title).

- subtitle:

  Character. Plot subtitle. Default: `NULL`.

- caption:

  Character. Plot caption. `NULL` produces an automatic caption
  reporting the p-value threshold and adjustment method. `""` suppresses
  the caption entirely. Default: `NULL`.

- facet_scales:

  Character. Passed to
  [`ggplot2::facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)'s
  `scales` argument when multiple terms are plotted. One of `"fixed"`,
  `"free_x"`, `"free_y"`, `"free"`. Default: `"fixed"`.

- facet_ncol:

  Integer or `NULL`. Number of columns in the facet grid. `NULL` lets
  `ggplot2` decide. Default: `NULL`.

- theme:

  Character string. One of `"nature"` (default, built on `theme_bw`),
  `"minimal"`, `"classic"`, `"bw"`.

- base_size:

  Base font size in points; every text element scales from it. Default:
  `11`, which suits a default 7 x 5 inch Quarto figure.

- font_family:

  Character. Font family. Default: `""` (system default).

## Value

A `ggplot` object (single-term) or a patchwork-composed plot (multi-term
with `facet_wrap`). The object can be passed to
[`save_plots()`](https://jllcalorio.github.io/pondeR/reference/save_plots.md)
or modified with standard `ggplot2` additions.

## Details

**Significance classes**: Each feature is assigned one of three classes
per term panel:

- **Up** - estimate \> `estimate_threshold` AND p (adjusted or raw) \<
  `p_threshold`.

- **Down** - estimate \< \\-\\`estimate_threshold` AND p \<
  `p_threshold`.

- **NS** - does not meet either criterion.

When `estimate_threshold = 0` (default), the up/down split is purely on
the sign of the estimate - useful for a single continuous predictor
(e.g., age, time as numeric) where there is no conventional fold-change
cutoff.

**ggrepel**: If ggrepel is installed (strongly recommended for dense
plots), it is used automatically. Install with
`install.packages("ggrepel")`. If absent, `geom_text` is used with a
slight upward nudge.

## References

Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
67(1), 1-48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *Journal of the
Royal Statistical Society: Series B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Okabe M, Ito K (2008). *Color Universal Design (CUD): How to make
figures and presentations that are friendly to colorblind people*.
<https://jfly.uni-koeln.de/color/>

## Author

John Lennon L. Calorio

## Examples

``` r
set.seed(42)
meta <- data.frame(subject_id = factor(rep(1:15, each = 2)),
                   time       = rep(c("pre", "post"), 15))
x <- data.frame(
  glucose = rnorm(30, ifelse(meta$time == "post", 5.6, 5), 0.4),
  insulin = rnorm(30, 80, 10),
  lactate = rnorm(30, 1.2, 0.3)
)
fit <- run_mixedmodel(x, meta, formula = ".feature ~ time + (1 | subject_id)",
                      verbose = FALSE)

plot_mixedmodel(fit)

plot_mixedmodel(fit, label_features = "glucose", label_box = TRUE)
```
