# pondeR

pondeR turns a data frame into publication-ready statistical tables and
figures with very little code. It picks a suitable test, reports effect
sizes and post hoc comparisons, and returns tables whose column names
can go straight into a report.

## Installation

``` r
# install.packages("pak")
pak::pak("jllcalorio/pondeR")
```

A few functions use optional packages. If one is missing, the function
tells you which package to install.

## Quick start

``` r
library(pondeR)

# Descriptive table with p-values
run_summarytable(iris, split_by = "Species", add_inferential_pvalues = TRUE)

# Automatic test selection, effect sizes and post hoc comparisons
res <- run_diff(iris, outcome = c("Sepal.Length", "Sepal.Width"),
                group = "Species", summary_table = TRUE)
res$summary_table

# A figure from the same result
plot_diff(res, outcome = "Sepal.Length")$plots[[1]]

# Correlations with confidence intervals
run_correl(mtcars[c("mpg", "wt", "hp")], summary_table_only = TRUE)
```

## Working with dimsprepr

pondeR accepts a [dimsprepr](https://github.com/jllcalorio/dimsprepr)
[`run_DIpreprocess()`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html)
result directly. The data, metadata, group column and QC samples are
picked up for you, and QC samples are left out of group comparisons.

``` r
library(dimsprepr)
library(pondeR)

di <- run_DIpreprocess(dims_data, dims_metadata, verbose = FALSE)

# Differential analysis
fc  <- run_foldchange(di)
dif <- run_diff(di, summary_table = TRUE)
plot_volcano(list(fc, dif))$plots[[1]]

# Multivariate analysis
pca <- run_pca(di)
plot_score(pca)
plot_scree(pca)
```

If you preprocessed with several methods at once, choose one result
first, for example `run_pca(di[["median_log10"]])`.

## Figures and tables in reports

Every plot defaults to `base_size = 11`, which suits the 7 x 5 inch
figures of a default Quarto document. Changing `base_size` rescales all
the text in a figure. Arguments such as `title_size` and `xlab_size`
override a single element.

Table columns already have report-ready names, such as
`Adjusted p-value` and `95% CI lower`, so tables can be printed without
renaming anything:

``` r
knitr::kable(res$Sepal.Length$data_summary, digits = 2)
as_gt(res$summary_table)
```

## What’s included

| Task | Functions |
|:---|:---|
| Summary tables | [`run_summarytable()`](https://jllcalorio.github.io/pondeR/reference/run_summarytable.md), [`run_anthroindex()`](https://jllcalorio.github.io/pondeR/reference/run_anthroindex.md) |
| Group comparisons | [`run_diff()`](https://jllcalorio.github.io/pondeR/reference/run_diff.md), [`plot_diff()`](https://jllcalorio.github.io/pondeR/reference/plot_diff.md), [`run_foldchange()`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md), [`plot_volcano()`](https://jllcalorio.github.io/pondeR/reference/plot_volcano.md), [`plot_fcmetric()`](https://jllcalorio.github.io/pondeR/reference/plot_fcmetric.md), [`run_lefse()`](https://jllcalorio.github.io/pondeR/reference/run_lefse.md) |
| Associations and correlations | [`run_assoc()`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md), [`run_single_assoc()`](https://jllcalorio.github.io/pondeR/reference/run_single_assoc.md), [`run_correl()`](https://jllcalorio.github.io/pondeR/reference/run_correl.md), [`plot_heatmap()`](https://jllcalorio.github.io/pondeR/reference/plot_heatmap.md), [`plot_dend()`](https://jllcalorio.github.io/pondeR/reference/plot_dend.md) |
| Modelling | [`run_logreg()`](https://jllcalorio.github.io/pondeR/reference/run_logreg.md), [`run_regreg()`](https://jllcalorio.github.io/pondeR/reference/run_regreg.md), [`run_mixedmodel()`](https://jllcalorio.github.io/pondeR/reference/run_mixedmodel.md), [`plot_mixedmodel()`](https://jllcalorio.github.io/pondeR/reference/plot_mixedmodel.md), [`plot_oddsratio()`](https://jllcalorio.github.io/pondeR/reference/plot_oddsratio.md), [`run_auc()`](https://jllcalorio.github.io/pondeR/reference/run_auc.md) |
| Multivariate analysis | [`run_pca()`](https://jllcalorio.github.io/pondeR/reference/run_pca.md), [`run_pls()`](https://jllcalorio.github.io/pondeR/reference/run_pls.md), [`run_pcoa()`](https://jllcalorio.github.io/pondeR/reference/run_pcoa.md), [`plot_score()`](https://jllcalorio.github.io/pondeR/reference/plot_score.md), [`plot_scree()`](https://jllcalorio.github.io/pondeR/reference/plot_scree.md) |
| Other figures | [`plot_bar()`](https://jllcalorio.github.io/pondeR/reference/plot_bar.md), [`plot_meanmap()`](https://jllcalorio.github.io/pondeR/reference/plot_meanmap.md), [`plot_multi_y()`](https://jllcalorio.github.io/pondeR/reference/plot_multi_y.md), [`plot_relabund()`](https://jllcalorio.github.io/pondeR/reference/plot_relabund.md) |
| Content validity | [`run_contentvalidity()`](https://jllcalorio.github.io/pondeR/reference/run_contentvalidity.md), [`run_contentvaliditycvr()`](https://jllcalorio.github.io/pondeR/reference/run_contentvaliditycvr.md) |
| Export and reports | [`save_excel()`](https://jllcalorio.github.io/pondeR/reference/save_excel.md), [`save_plots()`](https://jllcalorio.github.io/pondeR/reference/save_plots.md), [`render_metabolomics()`](https://jllcalorio.github.io/pondeR/reference/render_metabolomics.md) |
| Data helpers | [`run_relabel()`](https://jllcalorio.github.io/pondeR/reference/run_relabel.md), [`run_randomize()`](https://jllcalorio.github.io/pondeR/reference/run_randomize.md), [`run_subset()`](https://jllcalorio.github.io/pondeR/reference/run_subset.md), [`run_reduce()`](https://jllcalorio.github.io/pondeR/reference/run_reduce.md), [`run_vector2df()`](https://jllcalorio.github.io/pondeR/reference/run_vector2df.md), [`run_move_names()`](https://jllcalorio.github.io/pondeR/reference/run_move_names.md) |

See the [function
reference](https://jllcalorio.github.io/pondeR/reference/) for details.

## Contributing

Bug reports and feature requests are welcome on the [issue
tracker](https://github.com/jllcalorio/pondeR/issues). Please note that
this project has a [Code of
Conduct](https://github.com/jllcalorio/pondeR/blob/main/CODE_OF_CONDUCT.md).

## License

MIT
