# Plot PCA Scores Plot

Creates a scores plot showing sample positions in principal component
space. Supports discrete group coloring, continuous variable coloring,
confidence ellipses, and outlier detection.

## Usage

``` r
plot_score(
  pca_result,
  pc = c(1, 2),
  color_by,
  points_from,
  title = NULL,
  subtitle = NULL,
  position = "center",
  arrange_levels = NULL,
  ellipse = TRUE,
  ellipse_type = "t",
  ellipse_level = 0.95,
  legend = NULL,
  show_outliers = FALSE,
  label_outliers = "all",
  point_size = 3,
  point_alpha = 0.8,
  theme_base_size = 11,
  verbose = TRUE
)
```

## Arguments

- pca_result:

  List. Output from
  [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md).

- pc:

  Character vector of length 2. Principal components to plot, e.g.,
  c("PC1", "PC2"). Can also be numeric, e.g., c(1, 2). Default: c(1, 2).

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

  Character. Plot subtitle. Default: NULL (no subtitle).

- position:

  Character. Position of title and subtitle. Options: "left", "center",
  "right". Default: "center".

- arrange_levels:

  Character vector. Custom order for group levels in legend (only
  applies to categorical `color_by`). Default: NULL (alphabetical
  order).

- ellipse:

  Logical. Draw confidence ellipses around groups. Only applies when
  `color_by` is categorical. Requires at least 3 samples per group.
  Default: TRUE.

- ellipse_type:

  Character. Type of ellipse: "t" (t-distribution, default) or "norm"
  (normal distribution). Default: "t".

- ellipse_level:

  Numeric. Confidence level for ellipses (0-1). Default: 0.95 (95%
  confidence).

- legend:

  Character or NULL. Custom legend title. Default: NULL (uses `color_by`
  column name).

- show_outliers:

  Logical. Label samples outside the confidence ellipse defined by
  `ellipse_level`. Only works when ellipses are drawn (categorical
  coloring). Default: FALSE.

- label_outliers:

  Character. Which groups to check for outliers. Options: "all"
  (default) or specific group names from the `color_by` variable. Can be
  a vector.

- point_size:

  Numeric. Size of points. Default: 3.

- point_alpha:

  Numeric. Transparency of points (0-1). Default: 0.8.

- theme_base_size:

  Numeric. Base font size for theme. Default: 11.

- verbose:

  Logical. Print messages. Default: TRUE.

## Value

A ggplot2 object.

## Details

**Scores Plot Interpretation:**

- **Clustering**: Samples that cluster together have similar metabolic
  profiles

- **Separation**: Distance between groups indicates metabolic
  differences

- **Outliers**: Samples far from their group may represent:

  - Biological outliers (interesting phenotypes)

  - Technical issues (sample prep, analysis)

  - Mislabeled samples

**Coloring Options:**

The `color_by` parameter determines how points are colored:

1.  **Categorical column**: Colors by discrete groups (e.g., Treatment,
    Genotype)

    - Points are both colored and shaped by this variable

    - Confidence ellipses can be drawn with `ellipse = TRUE`

    - Outliers can be detected and labeled

2.  **Numeric column**: Gradient coloring by continuous variable (e.g.,
    Time, Dose, Age)

    - Uses viridis color scale

    - Points are shaped by the same variable (binned into quartiles)

    - No ellipses or outlier detection available

**Point Labeling:**

The `points_from` parameter controls how samples are labeled:

- If the specified column exists in metadata, those values are used as
  labels

- If the column doesn't exist, sequential labels (Sample_1, Sample_2,
  ...) are used

- Labels are shown when `show_outliers = TRUE` for outlier samples

- Common choices: "Sample", "SampleID", "SubjectID"

**Outlier Detection:**

Outliers are identified using Mahalanobis distance from the group
centroid. A sample is flagged as an outlier if it falls outside the
confidence ellipse defined by `ellipse_level`. The default (0.95)
identifies samples outside the 95% confidence region. This exactly
mirrors
[`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html)
behavior.

**When Outlier Detection Works:**

- Only with categorical `color_by`

- When `ellipse = TRUE`

- When groups have at least 3 samples

**Sample Exclusion:**

To exclude samples (e.g., QC samples) from the plot, exclude them during
PCA computation using `run_pca(exclude = ...)` rather than filtering in
the plot.

## References

Mardia, K.V., Kent, J.T., & Bibby, J.M. (1979). Multivariate Analysis.
London: Academic Press.

Brereton, R.G., & Lloyd, G.R. (2014). Partial least squares discriminant
analysis: taking the magic away. Journal of Chemometrics, 28(4),
213-225. [doi:10.1002/cem.2609](https://doi.org/10.1002/cem.2609)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(rnorm(100 * 50, mean = 100, sd = 20),
            nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)
rownames(x) <- paste0("Sample", 1:100)

metadata <- data.frame(
  Sample = paste0("Sample", 1:100),
  Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20)),
  Batch = rep(1:4, each = 25),
  Time = rep(c(0, 6, 12, 24), 25)
)

# Run PCA (excluding QC samples)
scaled_result <- run_scale(x, method = "auto")
pca_result <- run_pca(scaled_result$data, metadata, 
                      group = "Group",
                      exclude = "QC")

# Basic scores plot (colored by Group, labeled by Sample)
plot_score(pca_result, 
           color_by = "Group",
           points_from = "Sample")

# PC2 vs PC3 with custom title
plot_score(pca_result, 
           pc = c(2, 3),
           color_by = "Group",
           points_from = "Sample",
           title = "PCA: PC2 vs PC3",
           position = "left")

# Color by batch
plot_score(pca_result,
           color_by = "Batch",
           points_from = "Sample")

# Continuous coloring by time with custom legend
plot_score(pca_result,
           color_by = "Time",
           points_from = "Sample",
           legend = "Collection Time (h)")

# Show outliers in specific groups
plot_score(pca_result,
           color_by = "Group",
           points_from = "Sample",
           show_outliers = TRUE,
           label_outliers = c("Control", "Treatment"),
           ellipse_level = 0.95)

# Stricter outlier detection (99% confidence)
plot_score(pca_result,
           color_by = "Group",
           points_from = "Sample",
           show_outliers = TRUE,
           ellipse_level = 0.99)
} # }
```
