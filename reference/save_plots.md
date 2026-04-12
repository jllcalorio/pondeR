# Save Plots to Image or Vector Files

Exports a single plot or a named list of plots to individual files.
Supported plot types include ggplot2 objects, plotly widgets, grid
grobs/gtables, `recordedplot` objects, zero-argument base-R drawing
functions, and pondeR result objects produced by functions such as
[`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md),
[`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md),
[`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md),
[`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md),
and
[`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md).

When `x` is a single plot, `filename` is used as the output file name.
When `x` is a list of plots, `filename` is ignored and each element's
list name is used instead. If `x` contains no direct plot objects at the
top level but holds a list of plots one level deep (e.g.
`result$plots`), that nested list is transparently unpacked.

## Usage

``` r
save_plots(
  x,
  filename = "plot",
  format = NULL,
  width = 15,
  height = 12,
  dpi = 300,
  folder = NULL,
  transparent = FALSE,
  overwrite = FALSE,
  fd_date_stamp = FALSE,
  fd_time_stamp = FALSE,
  fn_date_stamp = FALSE,
  fn_time_stamp = FALSE
)
```

## Arguments

- x:

  A single plot object or a named list of plot objects. Accepted types:
  `ggplot`, `plotly`, `grob`/`gtable`, `recordedplot`, a zero-argument
  base-R drawing function, or a pondeR result object. One level of list
  nesting is automatically unpacked when no direct plot objects are
  found at the top level.

- filename:

  A single character string for the output file name. Used when `x` is a
  single plot. When `x` is a multi-plot list this argument is ignored
  and each element's name is used as the file name. The extension in
  `filename` sets the format unless `format` is also supplied. Defaults
  to `"plot"` when omitted.

- format:

  Optional single character string overriding the file extension for all
  exported files, e.g. `"pdf"`, `"svg"`, `"png"`. Case-insensitive;
  accepts with or without a leading dot. When supplied, the extension in
  `filename` is ignored.

- width:

  Numeric. Output width in inches. Default `15`.

- height:

  Numeric. Output height in inches. Default `12`.

- dpi:

  Numeric. Resolution in dots per inch for raster formats. Default
  `300`. Ignored for vector formats.

- folder:

  A single character string for the target directory. Supports nested
  folders (e.g., `"Main folder/sub folder"`). Created recursively if
  absent. Defaults to `NULL` (current working directory).

- transparent:

  Logical. Use a transparent background. Defaults to `FALSE`. A warning
  is issued and white is used as fallback for JPEG, which does not
  support transparency.

- overwrite:

  Logical. Overwrite existing files when `TRUE`. When `FALSE` (default)
  a Windows-style counter suffix is appended, e.g. `"plot (1).png"`.

- fd_date_stamp:

  Logical. Appends `"MmmDDYYYY"` date to `folder`. Default `FALSE`.

- fd_time_stamp:

  Logical. Appends `"HHMM"` time to `folder`. Default `FALSE`.

- fn_date_stamp:

  Logical. Appends date to each file name. Default `FALSE`.

- fn_time_stamp:

  Logical. Appends time to each file name. Default `FALSE`.

## Value

Called for its side effect. Returns the full normalised file path(s) of
the saved plot(s) invisibly as a character vector.

## Details

**Single vs. multi-plot behaviour.** A single plot object uses
`filename` exactly (after sanitisation). A list with multiple elements
ignores `filename`; each element's list name becomes its file name.
Illegal characters are stripped, and duplicate names are resolved with
`(1)`, `(2)`, etc.

**Unnamed list elements.** Elements without a name receive the
placeholder name `"plot_N"` where `N` is their position.

**pondeR result objects.** Objects whose class starts with `"plot_"` or
matches a known pondeR class are unpacked automatically. The first
extractable plot is exported; a note is printed when multiple plots are
present.

**plotly HTML export.** Requires htmlwidgets. Static raster export of
plotly objects requires plotly with a working reticulate/kaleido
installation: `reticulate::py_install("kaleido")`.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

# Single ggplot — filename is respected
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
save_plots(p, filename = "AUROC All (Neg)")

# List of ggplots — filename is ignored; list names are used
p2 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
        geom_point()
save_plots(list(scatter = p, iris_plot = p2))

# Accessing a list slot with multiple plots (auto-unpacked)
res <- run_auc(data = mydata, response = "group", predictors = c("A","B"))
save_plots(res$plots, folder = "figures")

# Override format for all plots
save_plots(list(scatter = p, iris_plot = p2), format = "pdf")

# Base-R plot via a zero-argument function
save_plots(list(hist_mpg = function() hist(mtcars$mpg, main = "MPG")))
} # }
```
