# Plot Multiple Variables on Primary and Secondary Y Axes

Creates a dual-axis ggplot with a shared x-axis, allowing simultaneous
display of multiple variables on independent primary and secondary y
scales. Supports four geometry types (line, area, histogram, dot) on
each axis, optional faceting, colorblind-safe palettes, and full
publication-ready typography controls.

## Usage

``` r
plot_multi_y(
  x,
  x_axis,
  main_as_line = NULL,
  main_as_area = NULL,
  main_as_hist = NULL,
  main_as_dot = NULL,
  secondary_as_line = NULL,
  secondary_as_area = NULL,
  secondary_as_hist = NULL,
  secondary_as_dot = NULL,
  order_main = c("hist", "area", "line", "dot"),
  order_secondary = c("hist", "area", "line", "dot"),
  color_main_line = NULL,
  color_main_area = NULL,
  color_main_hist = NULL,
  color_main_dot = NULL,
  color_secondary_line = NULL,
  color_secondary_area = NULL,
  color_secondary_hist = NULL,
  color_secondary_dot = NULL,
  linewidth_main_line = 0.8,
  linewidth_secondary_line = 0.8,
  size_main_dot = 2,
  size_secondary_dot = 2,
  area_alpha_main = 0.35,
  area_alpha_secondary = 0.35,
  facet = NULL,
  facet_ncol = NULL,
  legend_pos = "top",
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  x_lab = NULL,
  y_lab_main = NULL,
  y_lab_secondary = "Secondary axis",
  global_font_size = 20,
  title_font_size = NULL,
  subtitle_font_size = NULL,
  caption_font_size = NULL,
  axis_font_size = NULL,
  legend_font_size = NULL,
  facet_font_size = NULL,
  theme = "nature",
  zoom = 1
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, or `matrix` (converted internally). Column
  names may contain special characters.

- x_axis:

  A single character string giving the column of `x` to use on the
  x-axis. May be numeric, Date/POSIXct, or categorical.

- main_as_line:

  Character vector of column name(s) to plot as line geometry on the
  primary y-axis. `NULL` omits this geometry.

- main_as_area:

  Character vector of column name(s) to plot as area geometry on the
  primary y-axis. `NULL` omits this geometry.

- main_as_hist:

  Character vector of column name(s) to plot as histogram (bar) geometry
  on the primary y-axis. `NULL` omits this geometry.

- main_as_dot:

  Character vector of column name(s) to plot as point geometry on the
  primary y-axis. `NULL` omits this geometry.

- secondary_as_line:

  Character vector of column name(s) to plot as line geometry on the
  secondary y-axis. `NULL` omits this geometry.

- secondary_as_area:

  Character vector of column name(s) to plot as area geometry on the
  secondary y-axis. `NULL` omits this geometry.

- secondary_as_hist:

  Character vector of column name(s) to plot as histogram (bar) geometry
  on the secondary y-axis. `NULL` omits this geometry.

- secondary_as_dot:

  Character vector of column name(s) to plot as point geometry on the
  secondary y-axis. `NULL` omits this geometry.

- order_main:

  Character vector specifying the back-to-front layer order of primary
  geometries. Default `c("hist", "area", "line", "dot")`.

- order_secondary:

  Same as `order_main` but for secondary geometries.

- color_main_line:

  Character vector of R color names or hex codes for primary line
  geometry. Recycled if fewer colors than series. Defaults to the
  Okabe-Ito colorblind-safe palette.

- color_main_area:

  Character vector of colors for primary area geometry.

- color_main_hist:

  Character vector of colors for primary histogram geometry.

- color_main_dot:

  Character vector of colors for primary dot geometry.

- color_secondary_line:

  Character vector of colors for secondary line geometry.

- color_secondary_area:

  Character vector of colors for secondary area geometry.

- color_secondary_hist:

  Character vector of colors for secondary histogram geometry.

- color_secondary_dot:

  Character vector of colors for secondary dot geometry.

- linewidth_main_line:

  Numeric. Line width for primary line geometry. Default `0.8`.

- linewidth_secondary_line:

  Numeric. Line width for secondary line geometry. Default `0.8`.

- size_main_dot:

  Numeric. Point size for primary dot geometry. Default `2`.

- size_secondary_dot:

  Numeric. Point size for secondary dot geometry. Default `2`.

- area_alpha_main:

  Numeric in \\\[0, 1\]\\. Fill transparency for primary area and
  histogram geometries. Default `0.35`.

- area_alpha_secondary:

  Numeric in \\\[0, 1\]\\. Fill transparency for secondary area and
  histogram geometries. Default `0.35`.

- facet:

  A single character string naming a categorical column of `x` to use
  for
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
  `NULL` (default) disables faceting.

- facet_ncol:

  Integer. Number of columns when `facet` is set. `NULL` lets ggplot2
  choose automatically.

- legend_pos:

  Character. Legend position: `"top"` (default), `"bottom"`, `"left"`,
  or `"right"`.

- title:

  Character. Plot title. `NULL` omits the title.

- subtitle:

  Character. Plot subtitle. `NULL` omits the subtitle.

- caption:

  Character. Plot caption. `NULL` omits the caption.

- x_lab:

  Character. x-axis label. Defaults to the value of `x_axis`.

- y_lab_main:

  Character. Primary y-axis label. Defaults to `"Primary axis"` when
  secondary variables are supplied, otherwise an empty string.

- y_lab_secondary:

  Character. Secondary y-axis label. Defaults to `"Secondary axis"`.

- global_font_size:

  Numeric. Base font size in points. Default `20`.

- title_font_size:

  Numeric. Override font size for the plot title. Defaults to
  `1.20 * global_font_size`.

- subtitle_font_size:

  Numeric. Override font size for the subtitle. Defaults to
  `0.90 * global_font_size`.

- caption_font_size:

  Numeric. Override font size for the caption. Defaults to
  `0.75 * global_font_size`.

- axis_font_size:

  Numeric. Override font size for axis titles and text. Defaults to
  `1.00 * global_font_size`.

- legend_font_size:

  Numeric. Override font size for legend text. Defaults to
  `0.90 * global_font_size`.

- facet_font_size:

  Numeric. Override font size for facet strip labels. Only relevant when
  `facet` is non-`NULL`. Defaults to `0.85 * global_font_size`.

- theme:

  Character. Visual theme. One of `"nature"` (default), `"minimal"`, or
  `"classic"`.

- zoom:

  Numeric. Multiplicative scaling factor applied to all resolved sizes
  at construction time. Default `1`.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

**Dual-axis implementation.** ggplot2 implements secondary axes as
linear rescalings of the primary axis
([`sec_axis`](https://ggplot2.tidyverse.org/reference/sec_axis.html)).
`plot_multi_y` computes the linear transform mapping the secondary
variables' observed range into the primary range with a 5\\

**NA handling.** Rows where `x_axis` is `NA` are always dropped before
plotting. Rows where a y-variable is `NA` are retained so that other
series at the same x position are still drawn; ggplot2 renders gaps in
those series naturally. A warning is issued listing which columns
contain `NA` values so you remain aware of the gaps.

**Long-format pivoting.** Each geometry group is pivoted to long format
internally, with a `.series` factor carrying the original column name.
This is what feeds the colour/fill aesthetics and drives the legend,
ensuring the manual colour scale always has matching levels.

**Layer ordering.** `order_main` and `order_secondary` control the
painter's order. The default `c("hist", "area", "line", "dot")` places
histograms at the back, then areas, then lines, with dots in the
foreground.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Basic: area on primary, two lines on primary, faceted by Year ---
set.seed(42)
df <- data.frame(
  Year         = factor(rep(c("2021", "2022"), each = 20)),
  `Week No.`   = rep(1:20, 2),
  `Weekly cases` = c(rpois(20, 60), rpois(20, 80)),
  `RdRP GC/L`  = c(rep(NA, 5), runif(15, 1e4, 1e5),
                   rep(NA, 3), runif(17, 2e4, 2e5)),
  `N GC/L`     = c(rep(NA, 8), runif(12, 5e3, 8e4),
                   rep(NA, 4), runif(16, 1e4, 1e5)),
  check.names  = FALSE
)

plot_multi_y(
  x             = df,
  x_axis        = "Week No.",
  main_as_area  = "Weekly cases",
  main_as_line  = c("RdRP GC/L", "N GC/L"),
  facet         = "Year",
  title         = NULL
)

# --- Dual axis: cases on primary, viral load on secondary ---
plot_multi_y(
  x                  = df,
  x_axis             = "Week No.",
  main_as_area       = "Weekly cases",
  secondary_as_line  = c("RdRP GC/L", "N GC/L"),
  y_lab_main         = "Weekly cases",
  y_lab_secondary    = "GC/L",
  facet              = "Year"
)
} # }
```
