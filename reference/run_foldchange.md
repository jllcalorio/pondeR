# Fold Change Analysis Across Groups

Computes pairwise fold changes (and optionally log2 fold changes) for
each feature (column) in `x` across ordered group levels. Group means
are used as the basis of comparison. When `log2 = TRUE` the data are
shifted so that all values are strictly positive before
log-transformation. The function is designed to integrate seamlessly
with
[`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
and
[`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md).

## Usage

``` r
run_foldchange(
  x,
  metadata,
  group,
  arrange = NULL,
  filter = NULL,
  select = NULL,
  sort = TRUE,
  log2 = TRUE,
  eps = 1e-08
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, or named `matrix` whose columns are numeric
  features (variables) to be analysed. Column names may contain special
  characters.

- metadata:

  A `data.frame`, `tibble`, or named `matrix` with the same number of
  rows as `x` and containing sample annotations. Must include the column
  specified by `group`.

- group:

  A single character string naming the column in `metadata` that
  contains the grouping factor. Must have at least two distinct non-`NA`
  levels after optional filtering via `filter`.

- arrange:

  An optional character vector listing *all* levels of `group` (after
  filtering) in the desired comparison order. Fold changes are computed
  as consecutive ratios: `arrange[1] / arrange[2]`,
  `arrange[2] / arrange[3]`, etc. If `NULL` (default), levels are sorted
  alphanumerically.

- filter:

  An optional character vector of group levels to *exclude* before
  analysis. Rows whose `group` value appears in `filter` are dropped.
  Default `NULL` (no filtering).

- select:

  An optional character vector of column names in `x` to include in the
  analysis. Default `NULL` uses all columns of `x`.

- sort:

  Logical. If `TRUE` (default), each fold-change column in the output is
  sorted from highest to lowest absolute fold change.

- log2:

  Logical. If `TRUE` (default), log2 fold changes are also computed.
  Data are shifted prior to log-transformation so that the minimum value
  equals `eps` (see below), guaranteeing all values are strictly
  positive.

- eps:

  A small positive numeric value. When `log2 = TRUE`, data are shifted
  so that the global minimum of `x` (after filtering and selection)
  becomes `eps`. Default `1e-8`.

## Value

A named list of class `"run_foldchange"` containing:

- `fc_table`:

  A `data.frame` of pairwise fold changes. Rows are features; columns
  are labelled `"<level_i>_vs_<level_{i+1}>"`. When `sort = TRUE` each
  column is independently sorted (descending).

- `log2fc_table`:

  A `data.frame` of log2 fold changes in the same structure as
  `fc_table`. `NULL` when `log2 = FALSE`.

- `shifted_data`:

  The (possibly shifted) numeric data matrix used for log2 computation.
  `NULL` when `log2 = FALSE`.

- `min_value`:

  The global minimum of the analysis data *before* shifting (numeric
  scalar).

- `shift`:

  The numeric shift \\\delta\\ applied. `0` when no shift was needed or
  `log2 = FALSE`.

- `group_means`:

  A `data.frame` of per-group, per-feature means (unshifted) with one
  row per group level.

- `comparisons`:

  A character vector of comparison labels, e.g. `c("A_vs_B", "B_vs_C")`.

- `params`:

  A list of all resolved parameter values for reproducibility.

## Details

**Fold change computation**  
Group means are computed per feature. For \\k\\ ordered levels \\g_1,
g_2, \ldots, g_k\\, the fold change for comparison \\i\\ is
\\\bar{x}\_{g_i} / \bar{x}\_{g\_{i+1}}\\, yielding \\k - 1\\ pairwise
comparisons.

**Log2 shifting**  
When `log2 = TRUE`, let \\m\\ be the global minimum of the (filtered,
selected) data matrix. A shift of \\\delta = \text{eps} - m\\ is added
to every value when \\m \le 0\\; otherwise no shift is applied (\\\delta
= 0\\). Log2 fold changes are then computed from the shifted group
means.

**Integration with
[`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
and
[`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)**  
The returned list exposes `fc_table`, `log2fc_table`, and `shifted_data`
in a consistent format that
[`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md)
can consume directly.

## Author

John Lennon L. Calorio

## Examples

``` r
## -----------------------------------------------------------------------
## Example 1 — basic usage with iris
## -----------------------------------------------------------------------
data(iris)
x_iris  <- iris[, 1:4]
meta_iris <- data.frame(Species = iris$Species)

res <- run_foldchange(
  x        = x_iris,
  metadata = meta_iris,
  group    = "Species"
)
print(res$fc_table)
#>        feature setosa_vs_versicolor versicolor_vs_virginica
#> 1  Sepal.Width            1.2375451               0.9314055
#> 2 Sepal.Length            0.8433288               0.9010322
#> 3 Petal.Length            0.3431925               0.7672911
#> 4  Petal.Width            0.1855204               0.6544916
print(res$log2fc_table)
#>        feature setosa_vs_versicolor versicolor_vs_virginica
#> 1  Sepal.Width            0.3074811              -0.1025187
#> 2 Sepal.Length           -0.2458328              -0.1503495
#> 3 Petal.Length           -1.5429101              -0.3821541
#> 4  Petal.Width           -2.4303506              -0.6115534

## -----------------------------------------------------------------------
## Example 2 — custom order, filter one level, select two features
## -----------------------------------------------------------------------
res2 <- run_foldchange(
  x        = x_iris,
  metadata = meta_iris,
  group    = "Species",
  arrange  = c("virginica", "versicolor", "setosa"),
  filter   = NULL,
  select   = c("Sepal.Length", "Petal.Length")
)
print(res2$comparisons)
#> [1] "virginica_vs_versicolor" "versicolor_vs_setosa"   
print(res2$fc_table)
#>        feature virginica_vs_versicolor versicolor_vs_setosa
#> 1 Petal.Length                1.303286             2.913817
#> 2 Sepal.Length                1.109838             1.185777

## -----------------------------------------------------------------------
## Example 3 — data with special characters in column names
## -----------------------------------------------------------------------
set.seed(1)
n   <- 60
xsc <- data.frame(
  `LPC 18:2`  = rnorm(n, 5, 1),
  `PA O-28:1` = rnorm(n, 3, 0.5),
  check.names = FALSE
)
meta_sc <- data.frame(group = rep(c("Healthy", "Disease"), each = 30))

res3 <- run_foldchange(x = xsc, metadata = meta_sc, group = "group")
print(res3$fc_table)
#>     feature Disease_vs_Healthy
#> 1  LPC 18:2             1.0099
#> 2 PA O-28:1             1.0005
```
