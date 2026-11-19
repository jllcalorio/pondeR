# Scores Plot for Multivariate Ordination Results

Produces a publication-ready scores plot from the output of multivariate
ordination and dimensionality-reduction methods, including Principal
Component Analysis
([`run_pca`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)),
Partial Least Squares and PLS-DA
([`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)),
and Principal Coordinate Analysis
([`run_pcoa`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md)).
Supports discrete group coloring with confidence ellipses, continuous
gradient coloring, outlier detection and labeling, and a range of
publication themes.

## Usage

``` r
plot_score(res, ...)

# S3 method for class 'run_pca'
plot_score(
  res,
  pc = c(1, 2),
  color_by,
  points_from,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  position = "center",
  arrange_levels = NULL,
  ellipse = TRUE,
  ellipse_type = "t",
  ellipse_level = 0.95,
  legend = NULL,
  legend_position = "bottom",
  show_outliers = FALSE,
  label_outliers = "all",
  outlier_margin = 1,
  outlier_label_size = NULL,
  show_groups = NULL,
  colors = NULL,
  point_size = 10,
  point_alpha = 0.8,
  theme = "nature",
  base_size = 30,
  font_family = "sans",
  axis_title_size = NULL,
  axis_text_size = NULL,
  plot_title_size = NULL,
  legend_title_size = NULL,
  legend_text_size = NULL,
  zoom = 1,
  verbose = TRUE,
  ...
)

# S3 method for class 'run_pls'
plot_score(res, pc = c(1, 2), title = NULL, subtitle = NULL, ...)

# S3 method for class 'run_pcoa'
plot_score(
  res,
  pc = c(1, 2),
  color_by,
  points_from,
  metadata = NULL,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  position = "center",
  arrange_levels = NULL,
  ellipse = TRUE,
  ellipse_type = "t",
  ellipse_level = 0.95,
  legend = NULL,
  legend_position = "bottom",
  show_outliers = FALSE,
  label_outliers = "all",
  colors = NULL,
  point_size = 3,
  point_alpha = 0.8,
  theme = "nature",
  base_size = 11,
  font_family = "sans",
  axis_title_size = NULL,
  axis_text_size = NULL,
  plot_title_size = NULL,
  legend_title_size = NULL,
  legend_text_size = NULL,
  zoom = 1,
  verbose = TRUE,
  ...
)

# Default S3 method
plot_score(res, ...)
```

## Arguments

- res:

  List. Output from
  [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md),
  [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md),
  or
  [`run_pcoa()`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md).
  The appropriate plot method is selected automatically via S3 dispatch
  based on the class of `res`.

- ...:

  Additional arguments passed to methods.

- pc:

  Character vector of length 2. Principal components to plot, e.g.,
  `c("PC1", "PC2")`. Can also be numeric, e.g., `c(1, 2)`. Default:
  `c(1, 2)`.

- color_by:

  Character. Name of column in `metadata` to use for coloring points. If
  the column contains categorical data, points are colored by discrete
  groups with optional confidence ellipses. If the column contains
  numeric data, points are colored using a continuous gradient scale.
  Required parameter with no default.

- points_from:

  Character. Name of column in `metadata` to use for labeling points.
  This determines which values are used for sample labels in the plot
  and for outlier labeling when enabled. If the specified column is not
  found in metadata, sequential labels will be used instead. Required
  parameter with no default.

- title:

  Character. Main plot title. Default: auto-generated based on PCs.

- subtitle:

  Character. Plot subtitle. Default: `NULL` (no subtitle).

- caption:

  Character. Plot caption (bottom-right). Default: `NULL` (no caption).

- position:

  Character. Horizontal alignment of title, subtitle, and caption.
  Options: `"left"`, `"center"`, `"right"`. Default: `"center"`.

- arrange_levels:

  Character vector. Custom order for group levels in legend (only
  applies to categorical `color_by`). Default: `NULL` (alphabetical
  order).

- ellipse:

  Logical. Draw confidence ellipses around groups. Only applies when
  `color_by` is categorical. Requires at least 3 samples per group.
  Default: `TRUE`.

- ellipse_type:

  Character. Type of ellipse: `"t"` (t-distribution, default) or
  `"norm"` (normal distribution). Default: `"t"`.

- ellipse_level:

  Numeric. Confidence level for ellipses (0-1). Default: `0.95`.

- legend:

  Character or `NULL`. Custom legend title. Default: `NULL` (uses
  `color_by` column name).

- legend_position:

  Character. Position of the legend. Options: `"bottom"`, `"top"`,
  `"left"`, `"right"`, `"none"`. Default: `"bottom"`.

- show_outliers:

  Logical. Label samples outside the confidence ellipse defined by
  `ellipse_level`. Only works when ellipses are drawn (categorical
  coloring). Default: `FALSE`.

- label_outliers:

  Character. Which groups to check for outliers. Use `"all"` (default)
  or a character vector of specific group names from the `color_by`
  variable.

- outlier_margin:

  Numeric. Multiplier applied to the Mahalanobis squared-distance
  threshold that exactly matches the drawn `stat_ellipse` boundary (same
  covariance estimator and critical value formula). Default `1.0` labels
  any sample that falls outside that boundary. Values \> 1 require the
  sample to be proportionally farther beyond it: `1.0`–`1.5` labels most
  ellipse-breaching samples; `1.5`–`3.0` labels clearly separated
  samples; `> 3.0` labels only extreme outliers.

- outlier_label_size:

  Numeric or `NULL`. Font size (ggplot2 text units) for outlier labels.
  Default: `NULL` (resolves to `base_size / 3.5`, consistent with
  `base_size`).

- show_groups:

  Character vector or `NULL`. Subset of group values from `color_by` to
  include in the plot. Samples not in `show_groups` are hidden. Default:
  `NULL` (all groups plotted). Only applies when `color_by` is
  categorical. When supplied, axis labels are removed because variance
  explained is computed on the full dataset, not the shown subset.

- colors:

  Character vector. Custom color palette for discrete groups. If `NULL`
  (default), the Okabe-Ito colorblind-friendly palette is used. Ignored
  when `color_by` is numeric (continuous scale is always used for
  numeric variables).

- point_size:

  Numeric. Size of points. Default: `10`.

- point_alpha:

  Numeric. Transparency of points (0-1). Default: `0.8`.

- theme:

  Character. ggplot2 theme to apply. Options: `"nature"`, `"minimal"`,
  `"classic"`, `"bw"`, `"light"`, `"dark"`. Default: `"nature"` (a
  clean, publication-ready style based on `theme_bw()`).

- base_size:

  Numeric. Base font size for the theme (pts). Default: `30`.

- font_family:

  Character. Font family for all text elements. Default: `"sans"`.

- axis_title_size:

  Numeric. Font size for axis titles. Default: `base_size + 1`.

- axis_text_size:

  Numeric. Font size for axis tick labels. Default: `base_size - 1`.

- plot_title_size:

  Numeric. Font size for the plot title. Default: `base_size + 3`.

- legend_title_size:

  Numeric. Font size for the legend title. Default: `base_size`.

- legend_text_size:

  Numeric. Font size for legend text. Default: `base_size - 1`.

- zoom:

  Numeric. Scaling factor applied uniformly to all text and point size
  elements. Values \> 1 enlarge, values \< 1 shrink. Default: `1`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

- metadata:

  A `data.frame` with one row per observation (sample). Required for
  `plot_score.run_pcoa` because
  [`run_pcoa()`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md)
  does not bundle metadata inside its result (unlike
  [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md)).
  Row order must match the row order of the data matrix originally
  passed to
  [`run_pcoa()`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md).

## Value

A `ggplot2` object.

## Details

**Scores Plot Interpretation:**

- **Clustering**: Samples that cluster together have similar metabolic
  profiles.

- **Separation**: Distance between groups indicates metabolic
  differences.

- **Outliers**: Samples far from their group may represent biological
  outliers, technical issues, or mislabeled samples.

**Coloring Options:**

The `color_by` parameter determines how points are colored:

1.  **Categorical column**: Colors by discrete groups (e.g., Treatment,
    Genotype).

    - Points are both colored and shaped by this variable.

    - Confidence ellipses can be drawn with `ellipse = TRUE`.

    - Outliers can be detected and labeled.

    - Colors default to the Okabe-Ito colorblind-safe palette; override
      with `colors`.

2.  **Numeric column**: Gradient coloring by continuous variable (e.g.,
    Time, Dose, Age).

    - Uses the viridis "plasma" color scale.

    - Points are shaped by the same variable (binned into quartiles).

    - No ellipses or outlier detection available.

**Color Defaults:**

By default, discrete groups use the Okabe-Ito palette, which is
perceptually distinct and safe for the most common forms of color vision
deficiency (deuteranopia, protanopia, tritanopia). The palette supports
up to 8 groups. For more than 8 groups, supply a custom `colors` vector.

**Point Labeling:**

The `points_from` parameter controls how samples are labeled:

- If the specified column exists in metadata, those values are used as
  labels.

- If the column doesn't exist, sequential labels (`Sample_1`,
  `Sample_2`, ...) are used.

- Labels are shown when `show_outliers = TRUE` for outlier samples only.

- Common choices: `"Sample"`, `"SampleID"`, `"SubjectID"`.

**Outlier Detection:**

Outliers are identified using Mahalanobis distance from the group
centroid. A sample is flagged as an outlier if it falls outside the
confidence region defined by `ellipse_level`. The default (`0.95`)
identifies samples outside the 95% confidence region. This exactly
mirrors
[`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html)
behavior.

Outlier detection requires:

- Categorical `color_by`

- `ellipse = TRUE`

- At least 3 samples per group

**Sample Exclusion:**

To exclude samples (e.g., QC samples) from the plot, exclude them during
PCA computation using `run_pca(exclude = ...)` rather than filtering at
the plot stage.

**Themes:**

The `"nature"` theme (default) is a clean, publication-ready style with
a white background, minimal gridlines, and no top/right panel border -
suitable for journal figures. Other options map directly to their
ggplot2 equivalents.

## References

Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). *Multivariate Analysis*.
London: Academic Press.

Brereton, R.G., & Lloyd, G.R. (2014). Partial least squares discriminant
analysis: taking the magic away. *Journal of Chemometrics*, 28(4),
213-225. [doi:10.1002/cem.2609](https://doi.org/10.1002/cem.2609)

Okabe, M., & Ito, K. (2002). *Color Universal Design (CUD): How to make
figures and presentations that are friendly to colorblind people*.
<https://jfly.uni-koeln.de/color/>

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html),
[`run_pca`](https://jllcalorio.github.io/pondeR/reference/run_pca.md),
[`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)

## Author

John Lennon L. Calorio

## Examples

``` r
meta <- data.frame(Sample = rownames(iris), Species = iris$Species)
pca  <- run_pca(scale(iris[, 1:4]), meta, group = "Species", verbose = FALSE)

plot_score(pca, color_by = "Species", points_from = "Sample")
#> Using discrete coloring by 'Species'
#> Using labels from column 'Sample'
#> Scores plot created successfully.


# PC2 vs PC3, classic theme, outliers flagged
plot_score(pca, pc = c(2, 3), color_by = "Species", points_from = "Sample",
           show_outliers = TRUE, theme = "classic")
#> Using discrete coloring by 'Species'
#> Using labels from column 'Sample'
#> Computing outliers (ellipse level: 95.0%)...
#> Outliers detected: 12 (8.0%)
#> Scores plot created successfully.

```
