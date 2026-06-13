# Plot Distribution Comparison Before and After Data Transformation

Generates density and box plots comparing data distributions before and
after any transformation, normalization, scaling, or correction step.
Plots can be oriented by sample or by feature. Returns a named list of
individual `ggplot2` objects and a combined `gtable` grob.

## Usage

``` r
plot_dist_beforeafter(
  x,
  y,
  group_by = "sample",
  plot_what = NULL,
  n_random = 30L,
  x_label = "Before",
  y_label = "After",
  plot_title = NULL,
  x_fill = "#56B4E9",
  y_fill = "#E69F00",
  point_alpha = 0.6,
  theme = "nature",
  base_size = 15,
  font_family = "sans",
  plot_cols = NULL,
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
  transformation. Rows are samples; columns are features. Column names
  may contain special characters.

- y:

  A `data.frame`, `tibble`, or `matrix` representing the data *after*
  transformation. Must have the same dimensions and column names as `x`.

- group_by:

  Character. Perspective for box plots. `"sample"` shows one box per
  sample; `"feature"` shows one box per feature. Default is `"sample"`.

- plot_what:

  Character vector or `NULL`. Subset of column names (features) from `x`
  and `y` to include in the plots. If `NULL`, all features are used.
  Default is `NULL`.

- n_random:

  Integer or `NULL`. Maximum number of samples (when
  `group_by = "sample"`) or features (when `group_by = "feature"`) to
  display in box plots. If `NULL` or if the data has fewer items than
  `n_random`, all are shown. Default is `30`.

- x_label:

  Character. Label describing the before state, used in plot titles.
  Default is `"Before"`.

- y_label:

  Character. Label describing the after state, used in plot titles.
  Default is `"After"`.

- plot_title:

  Character or `NULL`. Main title for the combined plot. If `NULL`
  (default), a title is automatically generated.

- x_fill:

  Character. Fill color for before density and box plots. Default is
  `"#56B4E9"` (Okabe-Ito sky blue).

- y_fill:

  Character. Fill color for after density and box plots. Default is
  `"#E69F00"` (Okabe-Ito orange).

- point_alpha:

  Numeric in `[0, 1]`. Transparency for density fill and box fill.
  Default is `0.6`.

- theme:

  Character. ggplot2 theme to apply. Options: `"nature"`, `"minimal"`,
  `"classic"`, `"bw"`, `"light"`, `"dark"`. Default is `"nature"`.

- base_size:

  Numeric. Base font size for the theme (pts). Default is `15`.

- font_family:

  Character. Font family for all text elements. Default is `"sans"`.

- plot_cols:

  Integer or `NULL`. Number of columns in the combined output grid. If
  `NULL`, defaults to `2` (density row + boxplot row). Default is
  `NULL`.

- seed:

  Numeric or `NULL`. Random seed for reproducible sampling when
  `n_random` is used. Default is `123`.

- global_font_size:

  Numeric or `NULL`. Base font size for proportional scaling of all text
  elements. Default is `NULL`.

- title_size:

  Numeric or `NULL`. Plot panel title font size. Default is `NULL`.

- subtitle_size:

  Numeric or `NULL`. Combined plot main title font size. Default is
  `NULL`.

- xlab_size:

  Numeric or `NULL`. X-axis label font size. Default is `NULL`.

- ylab_size:

  Numeric or `NULL`. Y-axis label font size. Default is `NULL`.

## Value

A named list containing:

- `plot_density_before`:

  Density plot of `x`.

- `plot_density_after`:

  Density plot of `y`.

- `plot_box_before`:

  Box plot of `x`.

- `plot_box_after`:

  Box plot of `y`.

- `plot_combined`:

  A `gtable` grob of all four panels, drawn to the active device.

## Details

Plot Distribution Comparison Before and After Data Transformation

`x` and `y` must share identical row counts, column counts, and column
names. Row order is assumed to correspond between `x` and `y`.

When `plot_what` is specified, only those features are included in both
density and box plots. When `n_random` is also specified, random
sub-sampling is applied on top of the `plot_what` subset for box plots
only; density plots always use the full `plot_what` subset.

Box plot items are ordered by median value for readability.

Font sizes follow the same resolution priority as `plot_beforeafter`:
explicit argument \> `global_font_size` scaling \> `base_size` derived
default.

The function renders the combined plot to the active graphics device and
also returns a named list so individual panels can be extracted,
re-styled, or exported independently.

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/pondeR/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
n_samples  <- 40
n_features <- 20

x <- as.data.frame(
  matrix(abs(rnorm(n_samples * n_features, mean = 1000, sd = 300)),
         nrow = n_samples,
         dimnames = list(
           paste0("S", seq_len(n_samples)),
           paste0("Feature_", seq_len(n_features))
         ))
)
y <- as.data.frame(scale(x))

# Sample-wise comparison, default settings
result <- plot_dist_beforeafter(
  x       = x,
  y       = y,
  group_by = "sample"
)

# Feature-wise, subset of features, classic theme
result <- plot_dist_beforeafter(
  x         = x,
  y         = y,
  group_by  = "feature",
  plot_what = paste0("Feature_", 1:10),
  theme     = "classic"
)

# Extract individual panels
result$plot_density_before
result$plot_box_after
} # }
```
