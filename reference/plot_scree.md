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
  position = "center",
  what = "variance",
  type = "line",
  max_pc = 15,
  show_cumulative = TRUE,
  show_values = TRUE,
  theme_base_size = 11,
  verbose = TRUE
)
```

## Arguments

- pca_result:

  List. Output from
  [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md).

- title:

  Character. Main plot title. Default: "Scree Plot".

- subtitle:

  Character. Plot subtitle. Default: NULL (no subtitle).

- position:

  Character. Position of title and subtitle. Options: "left", "center",
  "right". Default: "center".

- what:

  Character. What to plot: "variance" (variance explained %) or "eigen"
  (eigenvalues). Default: "variance".

- type:

  Character. Plot type: "bar", "line", or "both". Default: "line".

- max_pc:

  Integer. Maximum number of PCs to display. Default: 15.

- show_cumulative:

  Logical. Show cumulative variance explained as a line on secondary
  axis. Only works when `what = "variance"`. Default: TRUE.

- show_values:

  Logical. Display values on plot using ggrepel to avoid overlap.
  Default: TRUE.

- theme_base_size:

  Numeric. Base font size for theme. Default: 11.

- verbose:

  Logical. Print messages. Default: TRUE.

## Value

A ggplot2 object.

## Details

**Scree Plot Interpretation:**

The scree plot helps determine how many principal components to retain:

- **Elbow Method**: Look for an "elbow" where the curve flattens. PCs
  before the elbow capture meaningful variance; those after capture
  mostly noise.

- **Kaiser Criterion**: Retain PCs with eigenvalues \> 1 (when
  `what = "eigen"`). These PCs explain more variance than a single
  original variable.

- **Cumulative Variance**: Retain enough PCs to explain 70-90% of total
  variance (shown when `show_cumulative = TRUE`).

**Plot Types:**

- **bar**: Traditional bar chart, good for discrete comparison

- **line**: Line plot, better for seeing trends

- **both**: Combines bars and line for comprehensive view

**Cumulative Variance Line:**

When `show_cumulative = TRUE`, a line showing cumulative variance is
added to the plot on a secondary y-axis at the top. This helps visualize
how many PCs are needed to explain a target percentage of variance
(e.g., 80%).

**When to Use Each Metric:**

- **Variance explained**: More intuitive (percentages). Recommended for
  most users.

- **Eigenvalues**: Useful for Kaiser criterion and advanced analysis.

## References

Cattell, R.B. (1966). The scree test for the number of factors.
Multivariate Behavioral Research, 1(2), 245-276.
[doi:10.1207/s15327906mbr0102_10](https://doi.org/10.1207/s15327906mbr0102_10)

Kaiser, H.F. (1960). The application of electronic computers to factor
analysis. Educational and Psychological Measurement, 20(1), 141-151.
[doi:10.1177/001316446002000116](https://doi.org/10.1177/001316446002000116)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Run PCA first
pca_result <- run_pca(scaled_data, metadata, scale_method = "auto")

# Basic scree plot
plot_scree(pca_result)

# Eigenvalue plot with custom styling
plot_scree(pca_result, 
          what = "eigen",
          type = "both",
          title = "PCA Eigenvalues",
          position = "left")

# Variance plot without cumulative line
plot_scree(pca_result,
          show_cumulative = FALSE,
          max_pc = 10)
} # }
```
