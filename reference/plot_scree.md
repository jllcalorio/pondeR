# Plot PCA Scree Plot

Creates a scree plot showing variance explained or eigenvalues per
principal component. Helps determine the number of meaningful PCs to
retain.

## Usage

``` r
plot_scree(
  pca_result,
  title = "Scree Plot",
  subtitle = NULL,
  caption = NULL,
  position = "center",
  what = "variance",
  type = "line",
  max_pc = 15,
  show_cumulative = TRUE,
  show_values = TRUE,
  bar_color = "#0072B2",
  line_color = "#0072B2",
  cumulative_color = "#D55E00",
  theme = "nature",
  base_size = 11,
  font_family = "sans",
  axis_title_size = NULL,
  axis_text_size = NULL,
  plot_title_size = NULL,
  zoom = 1,
  verbose = TRUE
)
```

## Arguments

- pca_result:

  List. Output from
  [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md).

- title:

  Character. Main plot title. Default: `"Scree Plot"`.

- subtitle:

  Character. Plot subtitle. Default: `NULL` (no subtitle).

- caption:

  Character. Plot caption (bottom-right). Default: `NULL` (no caption).

- position:

  Character. Horizontal alignment of title, subtitle, and caption.
  Options: `"left"`, `"center"`, `"right"`. Default: `"center"`.

- what:

  Character. What to plot: `"variance"` (variance explained %) or
  `"eigen"` (eigenvalues). Default: `"variance"`.

- type:

  Character. Plot geometry: `"bar"`, `"line"`, or `"both"`. Default:
  `"line"`.

- max_pc:

  Integer. Maximum number of PCs to display. Default: `15`.

- show_cumulative:

  Logical. Show cumulative variance explained as a dashed line with a
  secondary y-axis. Only available when `what = "variance"`. Default:
  `TRUE`.

- show_values:

  Logical. Display value labels on the plot using `ggrepel` to avoid
  overlap. Falls back to `geom_text` if `ggrepel` is not installed.
  Default: `TRUE`.

- bar_color:

  Character. Fill color for bars (when `type = "bar"` or `"both"`).
  Default: `"#0072B2"` (Okabe-Ito blue).

- line_color:

  Character. Color for the primary line and points (when `type = "line"`
  or `"both"`). Default: `"#0072B2"` (Okabe-Ito blue).

- cumulative_color:

  Character. Color for the cumulative variance line, points, and labels.
  Default: `"#D55E00"` (Okabe-Ito vermillion).

- theme:

  Character. ggplot2 theme to apply. Options: `"nature"`, `"minimal"`,
  `"classic"`, `"bw"`, `"light"`, `"dark"`. Default: `"nature"` (a
  clean, publication-ready style based on `theme_bw()` with no
  gridlines).

- base_size:

  Numeric. Base font size for the theme (pts). Default: `11`.

- font_family:

  Character. Font family for all text elements. Default: `"sans"`.

- axis_title_size:

  Numeric. Font size for axis titles. Default: `base_size + 1`.

- axis_text_size:

  Numeric. Font size for axis tick labels. Default: `base_size - 1`.

- plot_title_size:

  Numeric. Font size for the plot title. Default: `base_size + 3`.

- zoom:

  Numeric. Uniform scaling factor applied to all text sizes and
  point/line sizes. Values \> 1 enlarge, values \< 1 shrink. Default:
  `1`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A `ggplot2` object.

## Details

**Scree Plot Interpretation:**

The scree plot helps determine how many principal components to retain:

- **Elbow Method**: Look for an "elbow" where the curve flattens. PCs
  before the elbow capture meaningful variance; those after capture
  mostly noise.

- **Kaiser Criterion**: Retain PCs with eigenvalues \> 1 (when
  `what = "eigen"`). These PCs explain more variance than a single
  original variable.

- **Cumulative Variance**: Retain enough PCs to explain 70–90% of total
  variance (shown when `show_cumulative = TRUE`).

**Plot Types:**

- `"bar"`: Traditional bar chart; good for discrete comparison.

- `"line"`: Line plot; better for visualising trends and the elbow.

- `"both"`: Combines bars and a line for a comprehensive view.

**Cumulative Variance Line:**

When `show_cumulative = TRUE`, a dashed line showing cumulative variance
is overlaid on a shared y-axis scale. A secondary axis on the right is
labelled accordingly. Only available when `what = "variance"`.

**When to Use Each Metric:**

- **Variance explained**: More intuitive (percentages). Recommended for
  most users.

- **Eigenvalues**: Useful for the Kaiser criterion and advanced
  analysis.

**Themes:**

The `"nature"` theme (default) applies `theme_bw()` with all gridlines
removed — suitable for journal figures. Other options map directly to
their ggplot2 equivalents, also with gridlines suppressed.

## References

Cattell, R.B. (1966). The scree test for the number of factors.
*Multivariate Behavioral Research*, 1(2), 245–276.
[doi:10.1207/s15327906mbr0102_10](https://doi.org/10.1207/s15327906mbr0102_10)

Kaiser, H.F. (1960). The application of electronic computers to factor
analysis. *Educational and Psychological Measurement*, 20(1), 141–151.
[doi:10.1177/001316446002000116](https://doi.org/10.1177/001316446002000116)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Run PCA first
pca_result <- run_pca(scaled_data, metadata, scale_method = "auto")

# Basic scree plot (nature theme, Okabe-Ito colors, defaults)
plot_scree(pca_result)

# Eigenvalue plot, bar + line geometry
plot_scree(pca_result,
           what     = "eigen",
           type     = "both",
           title    = "PCA Eigenvalues",
           position = "left")

# Variance plot without cumulative line, classic theme
plot_scree(pca_result,
           show_cumulative = FALSE,
           max_pc          = 10,
           theme           = "classic")

# Larger font, serif family, custom zoom
plot_scree(pca_result,
           base_size   = 13,
           font_family = "serif",
           zoom        = 1.2)
} # }
```
