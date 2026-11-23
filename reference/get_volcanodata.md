# Get Volcano Plot Data

Extracts and combines results from
[`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
and
[`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
into a single data frame. This function automatically aligns
features/outcomes, retrieves fold-change and significance thresholds
from the provided objects, and (optionally) keeps only "Up" or "Down"
regulated features.

## Usage

``` r
get_volcanodata(
  ...,
  up = NULL,
  down = NULL,
  pval = NULL,
  filter = TRUE,
  sort_es = c("desc", "asc", "none")
)
```

## Arguments

- ...:

  Objects of class
  [`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  and/or
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md).

- up, down, pval:

  Optional numeric overrides for thresholds. If `NULL` (default),
  thresholds are retrieved from the objects' parameters.

- filter:

  Logical. If `TRUE` (default), the data frame is filtered to include
  only "Up" or "Down" regulated features.

- sort_es:

  Row order by effect size: `"desc"` (default) puts the largest effects
  first, `"asc"` the smallest first, and `"none"` keeps the original
  order. Features are ordered by effect size magnitude (very large,
  large, medium, small, very small) and then by the absolute effect size
  estimate. Ignored when no
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  effect sizes are available. Note that estimates from different effect
  size metrics (e.g. Cohen's d vs. rank-biserial correlation) are not
  strictly comparable.

## Value

A `data.frame` containing merged data from both objects: fold change
columns, effect size (95\\ p-value column (`Adjusted p-value` when
present, otherwise `p-value`). Column names are ready for reports. When
the
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
summary table is simplified, the p-value column keeps its superscript
test markers (legend stored in the `"pvalue_legend"` attribute and shown
by [`as_gt`](https://jllcalorio.github.io/pondeR/reference/as_gt.md)).
Up/Down/NS classification is used only for `filter` and is not returned.

## Details

The function automatically detects the thresholds:

- `up` and `down` are retrieved from the
  [`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  parameters.

- Significance threshold (`pval`) is retrieved from
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  parameters (defaulting to 0.05 if not found).

- When
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  result has multiple outcomes, an `Adjusted p-value` column is present;
  this function prefers it over `p-value`.

## References

Cohen J (1988). *Statistical Power Analysis for the Behavioral Sciences*
(2nd ed.). Lawrence Erlbaum Associates.

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *Journal of the
Royal Statistical Society: Series B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## See also

[`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md),
[`plot_volcano`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md),
[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html)

## Author

John Lennon L. Calorio

## Examples

``` r
x    <- iris[iris$Species != "setosa", 1:4]
meta <- droplevels(iris[iris$Species != "setosa", "Species", drop = FALSE])

fc  <- run_foldchange(x, metadata = meta, group = "Species")
dif <- run_diff(x, meta, outcome = names(x), group = "Species",
                summary_table = TRUE, verbose = FALSE)

get_volcanodata(fc, dif, filter = FALSE)
#>                Comparison      Feature Mean (group 1) Mean (group 2)
#> 4 versicolor vs virginica  Petal.Width          1.326          2.026
#> 3 versicolor vs virginica Petal.Length          4.260          5.552
#> 2 versicolor vs virginica Sepal.Length          5.936          6.588
#> 1 versicolor vs virginica  Sepal.Width          2.770          2.974
#>   Fold change    log2 FC          ES (95% CI)  Magnitude Adjusted p-value
#> 4   0.6544916 -0.6115534 -0.96 (-0.97, -0.94) very large           <.001ᵇ
#> 3   0.7672911 -0.3821541 -2.52 (-3.04, -1.99)      large           <.001ᵃ
#> 2   0.9010322 -0.1503495 -1.13 (-1.55, -0.70)      large           <.001ᵃ
#> 1   0.9314055 -0.1025187 -0.33 (-0.51, -0.11)      large            .005ᵇ
```
