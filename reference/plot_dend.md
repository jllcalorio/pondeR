# Plot a Hierarchical Clustering Dendrogram

Produces a publication-ready dendrogram from a data frame, tibble, or
matrix using hierarchical clustering. Supports row-wise and column-wise
clustering with configurable distance metrics and agglomeration methods.
Returns a ggplot2 object via ggdendro.

## Usage

``` r
plot_dend(
  x,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  agglomeration_method = "ward.D",
  orientation = "rows",
  theme = "nature",
  title = NULL,
  subtitle = NULL,
  base_size = 11,
  title_size = NULL,
  subtitle_size = NULL,
  label_size = NULL,
  leaf_label_angle = NULL,
  leaf_labels = NULL,
  line_width = 0.5,
  line_color = "black",
  ...
)
```

## Arguments

- x:

  A data frame, tibble, or matrix of numeric values to cluster. Column
  names may contain special characters.

- clustering_distance_rows:

  A single character string specifying the distance metric for row
  clustering. One of `"euclidean"` (default), `"maximum"`,
  `"manhattan"`, `"canberra"`, `"binary"`, or `"minkowski"`. See
  [`dist`](https://rdrr.io/r/stats/dist.html).

- clustering_distance_cols:

  A single character string specifying the distance metric for column
  clustering. Same options as `clustering_distance_rows`. Default
  `"euclidean"`.

- agglomeration_method:

  A single character string specifying the agglomeration method passed
  to [`hclust`](https://rdrr.io/r/stats/hclust.html). One of `"ward.D"`
  (default), `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
  `"mcquitty"`, `"median"`, or `"centroid"`.

- orientation:

  A single character string. One of `"rows"` (default; clusters rows),
  `"cols"` (clusters columns), or `"both"` (returns a named list of two
  `ggplot2` objects).

- theme:

  A single character string specifying the plot theme. One of `"nature"`
  (default), `"classic"`, `"plos"`, `"minimal"`, `"bw"`, `"light"`,
  `"dark"`, or `"gray"`.

- title:

  A single character string for the plot title. Default `NULL` (no
  title).

- subtitle:

  A single character string for the plot subtitle. Default `NULL` (no
  subtitle).

- base_size:

  Base font size in points; every text element scales from it. Default
  `11`, which suits a default 7 x 5 inch Quarto figure.

- title_size, subtitle_size, label_size:

  Font sizes in points for the title, subtitle and leaf labels. Default
  `NULL` derives each from `base_size`.

- leaf_label_angle:

  A single numeric specifying the rotation angle (in degrees) for leaf
  label text. Default `NULL`. When `NULL`, the angle defaults to 0 if
  `orientation = "rows"` (horizontal labels) or 90 if
  `orientation = "cols"` (vertical labels).

- leaf_labels:

  A character vector of custom leaf names. Default `NULL`. Must have one
  name per row for `orientation = "rows"` and one per column for
  `orientation = "cols"`; not used with `orientation = "both"`.

- line_width:

  A single positive numeric for dendrogram branch width. Default `0.5`.

- line_color:

  A single character string for branch colour. Default `"black"`.

- ...:

  Additional arguments passed to
  [`dist`](https://rdrr.io/r/stats/dist.html) or
  [`hclust`](https://rdrr.io/r/stats/hclust.html).

## Value

A `ggplot2` object, or a named list of two `ggplot2` objects (`rows` and
`cols`) when `orientation = "both"`.

## Details

Numeric columns are extracted and scaled (zero mean, unit variance)
before computing distances. Rows or columns containing all `NA` values
are silently dropped. Remaining `NA`s are imputed with column means
before scaling.

When `orientation = "both"`, the returned object is a named list with
elements `rows` and `cols`.

## References

Ward JH (1963). Hierarchical grouping to optimize an objective function.
*Journal of the American Statistical Association*, 58(301), 236-244.
[doi:10.1080/01621459.1963.10500845](https://doi.org/10.1080/01621459.1963.10500845)

Murtagh F, Contreras P (2012). Algorithms for hierarchical clustering:
an overview. *WIREs Data Mining and Knowledge Discovery*, 2(1), 86-97.
[doi:10.1002/widm.53](https://doi.org/10.1002/widm.53)

## Author

John Lennon L. Calorio

## Examples

``` r
## Row dendrogram (default)
plot_dend(mtcars)


## Column dendrogram with Manhattan distance
plot_dend(mtcars, orientation = "cols",
          clustering_distance_cols = "manhattan")


## Both axes
dends <- plot_dend(mtcars, orientation = "both")
dends$rows

dends$cols


## Custom styling
plot_dend(mtcars, theme = "minimal", line_color = "#2166AC",
          line_width = 0.8, base_size = 10,
          title = "Sample Clustering")
```
