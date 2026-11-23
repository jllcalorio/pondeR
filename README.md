# pondeR

pondeR turns a data frame into publication-ready statistical tables and figures
with very little code. It picks a suitable test, reports effect sizes and post
hoc comparisons, and returns tables whose column names can go straight into a
report.

## Installation

```r
# install.packages("pak")
pak::pak("jllcalorio/pondeR")
```

A few functions use optional packages. If one is missing, the function tells you
which package to install.

## Quick start

```r
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
`run_DIpreprocess()` result directly. The data, metadata, group column and QC
samples are picked up for you, and QC samples are left out of group comparisons.

```r
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

If you preprocessed with several methods at once, choose one result first, for
example `run_pca(di[["median_log10"]])`.

## Figures and tables in reports

Every plot defaults to `base_size = 11`, which suits the 7 x 5 inch figures of a
default Quarto document. Changing `base_size` rescales all the text in a figure.
Arguments such as `title_size` and `xlab_size` override a single element.

Table columns already have report-ready names, such as `Adjusted p-value` and
`95% CI lower`, so tables can be printed without renaming anything:

```r
knitr::kable(res$Sepal.Length$data_summary, digits = 2)
as_gt(res$summary_table)
```

## What's included

| Task | Functions |
|:--|:--|
| Summary tables | `run_summarytable()`, `run_anthroindex()` |
| Group comparisons | `run_diff()`, `plot_diff()`, `run_foldchange()`, `plot_volcano()`, `plot_fcmetric()`, `run_lefse()` |
| Associations and correlations | `run_assoc()`, `run_single_assoc()`, `run_correl()`, `plot_heatmap()`, `plot_dend()` |
| Modelling | `run_logreg()`, `run_regreg()`, `run_mixedmodel()`, `plot_mixedmodel()`, `plot_oddsratio()`, `run_auc()` |
| Multivariate analysis | `run_pca()`, `run_pls()`, `run_pcoa()`, `plot_score()`, `plot_scree()` |
| Other figures | `plot_bar()`, `plot_meanmap()`, `plot_multi_y()`, `plot_relabund()` |
| Content validity | `run_contentvalidity()`, `run_contentvaliditycvr()` |
| Export and reports | `save_excel()`, `save_plots()`, `render_metabolomics()` |
| Data helpers | `run_relabel()`, `run_randomize()`, `run_subset()`, `run_reduce()`, `run_vector2df()`, `run_move_names()` |

See the [function reference](https://jllcalorio.github.io/pondeR/reference/) for
details.

## Contributing

Bug reports and feature requests are welcome on the
[issue tracker](https://github.com/jllcalorio/pondeR/issues). Please note that
this project has a
[Code of Conduct](https://github.com/jllcalorio/pondeR/blob/main/CODE_OF_CONDUCT.md).

## License

MIT
