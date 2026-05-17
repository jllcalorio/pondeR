# Plot Results from a Linear Mixed-Effects Model Analysis

Produces a volcano plot from the `combined_fixed_effects` component of a
`run_mixedmodel` object, with log2 fold-change (or raw estimates) on the
x-axis and \\-\log\_{10}(p\_{\text{adj}})\\ on the y-axis. When the
model formula contained more than one fixed effect term, panels are
automatically faceted by `term` so each predictor is shown separately.
Points are colored by significance and direction, and labels can be
added to the top-\\n\\ features per panel or to a user-specified set of
features.

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
  label_size = 3,
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
  base_size = 12,
  font_family = "",
  zoom = 1
)
```

## Arguments

- x:

  A `run_mixedmodel` object returned by
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
  `run_mixedmodel`). If `FALSE`, raw `p_value` is used.

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

  Numeric. Font size for point labels (in ggplot2 pt units). Default:
  `3`.

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

  Numeric. Base font size passed to the theme. Default: `12`.

- font_family:

  Character. Font family. Default: `""` (system default).

- zoom:

  Numeric scalar. Multiplies all point sizes, label sizes, line widths,
  and base font size proportionally for export scaling. Default: `1`.

## Value

A `ggplot` object (single-term) or a patchwork-composed plot (multi-term
with `facet_wrap`). The object can be passed to
[`save_plots()`](https://jllcalorio.github.io/pondeR/reference/save_plots.md)
or modified with standard `ggplot2` additions.

## Details

**Significance classes**: Each feature is assigned one of three classes
per term panel:

- **Up** — estimate \> `estimate_threshold` AND p (adjusted or raw) \<
  `p_threshold`.

- **Down** — estimate \< \\-\\`estimate_threshold` AND p \<
  `p_threshold`.

- **NS** — does not meet either criterion.

When `estimate_threshold = 0` (default), the up/down split is purely on
the sign of the estimate — useful for a single continuous predictor
(e.g., age, time as numeric) where there is no conventional fold-change
cutoff.

**ggrepel**: If ggrepel is installed (strongly recommended for dense
plots), it is used automatically. Install with
`install.packages("ggrepel")`. If absent, `geom_text` is used with a
slight upward nudge.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_mixedmodel(
  x        = x,
  metadata = metadata,
  formula  = ".feature ~ time + age + (1 | subject_id)"
)

# Default: all terms, auto-label top 10 significant features
plot_mixedmodel(result)

# Focus on one term, label specific features
plot_mixedmodel(
  result,
  term           = "timepost",
  label_features = c("glucose", "insulin"),
  color_up       = "#CC79A7",
  alpha          = 0.3,
  alpha_sig      = 0.9,
  title          = "Time effect (pre → post)"
)

# No auto-labels, only highlight manually specified features
plot_mixedmodel(
  result,
  label_features = c("glucose", "insulin", "PA O-28:1"),
  label_box      = TRUE,
  label_color    = "black"
)
} # }
```
