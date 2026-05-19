# Plot Pairwise Before-and-After Comparison for Selected Features

Generates a multi-panel plot comparing the signal of selected features
before and after any data transformation, correction, or manipulation.
Each panel displays one feature across the injection sequence, colored
by sample group and shaped by batch, with an optional LOESS trend line
fitted to a user-defined reference group (e.g., QC samples).

## Usage

``` r
plot_beforeafter(
  x,
  y,
  metadata,
  plot_what = NULL,
  col_injection = "InjectionSequence",
  col_batch = "Batch",
  col_group = "Group",
  col_qc_label = "QC",
  x_label = "Before",
  y_label = "After",
  log_transform = TRUE,
  point_size = 2,
  point_alpha = 0.7,
  plot_cols = NULL,
  theme = "nature",
  base_size = 11,
  font_family = "sans",
  seed = 123,
  global_font_size = NULL,
  title_size = NULL,
  subtitle_size = NULL,
  xlab_size = NULL,
  ylab_size = NULL
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, or `matrix` representing the data *before*
  the transformation. Rows are samples; columns are features. Column
  names may contain special characters.

- y:

  A `data.frame`, `tibble`, or `matrix` representing the data *after*
  the transformation. Must have the same dimensions and column names as
  `x`.

- metadata:

  A `data.frame` containing sample-level metadata shared by both `x` and
  `y`. Must have the same number of rows as `x` and `y`.

- plot_what:

  Character vector or `NULL`. Column names from `x` to be plotted. If
  `NULL` (default), the function will randomly select up to 6 feature
  names from `x`. If `x` has fewer than 6 features, all available
  features will be plotted. This argument is required.

- col_injection:

  Character. Name of the column in `metadata` containing the injection
  sequence (numeric order of runs). Default is `"InjectionSequence"`.

- col_batch:

  Character. Name of the column in `metadata` containing batch
  identifiers. Default is `"Batch"`.

- col_group:

  Character. Name of the column in `metadata` containing sample group
  labels (e.g., `"QC"`, `"Sample"`). Default is `"Group"`.

- col_qc_label:

  Character. The value in the `col_group` column that identifies the
  reference group for LOESS trend fitting (typically QC samples). If
  `NULL`, no LOESS line is drawn. Default is `"QC"`.

- x_label:

  Character. Label for the before-correction facet panel. Default is
  `"Before"`.

- y_label:

  Character. Label for the after-correction facet panel. Default is
  `"After"`.

- log_transform:

  Logical. If `TRUE`, applies a `log10` transformation (with a small
  floor of `1e-10`) to feature values before plotting. Default is
  `TRUE`.

- point_size:

  Numeric. Size of individual data points. Default is `2`.

- point_alpha:

  Numeric in `[0, 1]`. Transparency of data points. Default is `0.7`.

- plot_cols:

  Integer or `NULL`. Number of columns in the combined plot grid. If
  `NULL`, determined automatically. Default is `NULL`.

- theme:

  Character. ggplot2 theme to apply to each panel. Options: `"nature"`,
  `"minimal"`, `"classic"`, `"bw"`, `"light"`, `"dark"`. Default is
  `"nature"` (a clean, publication-ready style based on `theme_bw()`).

- base_size:

  Numeric. Base font size for the theme (pts). Default is `11`.

- font_family:

  Character. Font family for all text elements. Default is `"sans"`.

- seed:

  Numeric or `NULL`. Random seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before any internal
  sampling operations. Default is `123`.

- global_font_size:

  Numeric or `NULL`. Base font size for proportional scaling of all text
  elements. Default is `NULL`.

- title_size:

  Numeric or `NULL`. Feature panel title font size. Default is `NULL`.

- subtitle_size:

  Numeric or `NULL`. Main plot title (top grob) font size. Default is
  `NULL`.

- xlab_size:

  Numeric or `NULL`. X-axis label font size. Default is `NULL`.

- ylab_size:

  Numeric or `NULL`. Y-axis label font size. Default is `NULL`.

## Value

Invisibly returns a `gtable` grob object (from
[`gridExtra::arrangeGrob()`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html))
containing all feature panels. The plot is also drawn to the active
graphics device.

## Details

Plot Pairwise Before-and-After Comparison for Selected Features

`x` and `y` must share identical row counts, column counts, and column
names. Row order is assumed to correspond between `x`, `y`, and
`metadata`.

When `log_transform = TRUE`, values are floored at `1e-10` before log
transformation to avoid `-Inf` results from zero or negative values.

Font sizes are resolved with the following priority (highest to lowest):
explicit size argument (e.g., `title_size`) \> `global_font_size`
scaling \> `base_size`-derived default. When `global_font_size` is
provided and a specific size is not, text elements scale proportionally
relative to a base size of 11.

The LOESS trend line is fitted only to the group identified by
`col_qc_label`. If fewer than 3 observations of that group are present,
a warning is issued and the line is suppressed for that feature.

The combined plot is rendered via
[`gridExtra::arrangeGrob()`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
and displayed with
[`grid::grid.draw()`](https://rdrr.io/r/grid/grid.draw.html). The
returned object is an invisible `gtable` grob that can be redrawn or
exported with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
n_samples  <- 30
n_features <- 4

x <- as.data.frame(
  matrix(abs(rnorm(n_samples * n_features, mean = 1000, sd = 200)),
         nrow = n_samples,
         dimnames = list(NULL, paste0("Feature_", seq_len(n_features))))
)
y <- as.data.frame(x + rnorm(n_samples * n_features, sd = 50))

meta <- data.frame(
  InjectionSequence = seq_len(n_samples),
  Batch             = rep(c("B1", "B2"), each = n_samples / 2),
  Group             = rep(c("Sample", "QC"), times = c(n_samples - 6, 6))
)

plot_beforeafter(
  x             = x,
  y             = y,
  metadata      = meta,
  plot_what     = c("Feature_1", "Feature_2"),
  col_injection = "InjectionSequence",
  col_batch     = "Batch",
  col_group     = "Group",
  col_qc_label  = "QC",
  theme         = "nature"
)
} # }
```
