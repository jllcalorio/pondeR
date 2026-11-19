# Comprehensive Area Under the ROC Curve (AUROC) Analysis

Performs a comprehensive Receiver Operating Characteristic (ROC)
analysis for one or more continuous predictors against a binary response
variable. Wraps [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html)
to compute AUC and confidence intervals, and optionally produces
publication-ready ROC plots via ggplot2.

## Usage

``` r
run_auc(
  x,
  y = NULL,
  remove = NULL,
  metadata = NULL,
  set_control = NULL,
  include_cols = NULL,
  include_rows = NULL,
  direction = "auto",
  compute_auc = TRUE,
  compute_ci = TRUE,
  ci_level = 0.95,
  compute_performance_metrics = TRUE,
  plot = "none",
  theme = "nature",
  ncol = NULL,
  nrow = NULL,
  plot_title = NULL,
  plot_subtitle = NULL,
  plot_legend = "AUC (95% CI)",
  xlab = NULL,
  ylab = NULL,
  linewidth = 3,
  global_font_size = 25,
  title_font_size = NULL,
  subtitle_font_size = NULL,
  xlab_font_size = NULL,
  ylab_font_size = NULL,
  inplot_font_size = NULL,
  min_vip = 1,
  up = 1.5,
  down = 0.6,
  p_value = 0.05,
  ...
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, named `matrix`, `run_DIpreprocess` object,
  or a `list` of objects from the `run_*` ecosystem. If a list is
  provided, it MUST contain a `run_DIpreprocess` object. Optional
  inclusion of
  [`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md),
  [`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md),
  and/or
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  objects will trigger automatic feature filtering based on `min_vip`,
  `up`, `down`, and `p_value` thresholds.

- y:

  A single character string naming the column in `x` that contains the
  binary grouping variable (cases vs. controls). The column may be
  categorical (`character` / `factor`) with exactly two levels after
  optional removal via `remove`, or numeric containing only `0`
  (control) and `1` (case) after optional removal. If `x` is a
  `run_DIpreprocess` object or a list containing one, `y` will default
  to `x$parameters$group_col` if not explicitly provided.

- remove:

  Optional. A vector of values to remove from `y` before analysis, used
  to reduce a multi-level variable down to exactly two. If `x` is a
  `run_DIpreprocess` object or a list containing one, `remove` will
  default to `x$parameters$qc_types` if not explicitly provided. May be
  a character vector (for categorical `y`), groups. May be a character
  vector (for categorical `y`), a numeric vector (for numeric `y`), or a
  mixed list. Ignored when `y` already contains exactly two categories
  or only `0`/`1` values. Default `NULL`.

- metadata:

  An optional object of the same class and row count as `x` carrying
  additional sample annotations. Currently reserved for future use and
  is ignored by the analysis. Default `NULL`.

- set_control:

  A single string specifying which level of a categorical `y` represents
  the control group. Ignored when `y` is numeric. Default `NULL` (the
  first level alphabetically or by factor order).

- include_cols:

  A character vector of column names, or an integer vector of column
  positions, selecting the predictor columns of `x` to be analysed.
  Default `NULL` includes all non-`y` columns.

- include_rows:

  A character vector of row names, or an integer vector of row
  positions, selecting the rows of `x` to be analysed. Row names are
  only considered when `rownames(x)` are not the default integer
  sequence. Default `NULL` includes all rows.

- direction:

  A single string controlling the direction of comparison passed to
  [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html). Choices:
  `"auto"` (default), `">"` (controls \> cases), or `"<"` (controls \<
  cases). See [`?pROC::roc`](https://rdrr.io/pkg/pROC/man/roc.html) for
  details.

- compute_auc:

  Logical. Whether to compute the AUC. Default `TRUE`.

- compute_ci:

  Logical. Whether to compute bootstrap or DeLong confidence intervals
  for the AUC. Default `TRUE`.

- ci_level:

  A numeric scalar in \\(0, 1)\\ specifying the confidence level.
  Default `0.95`.

- compute_performance_metrics:

  Logical. Whether to compute additional performance metrics
  (sensitivity, specificity, NPV, PPV, accuracy, Youden index, J-stat)
  at the optimal threshold determined by Youden's index. Default `TRUE`.

- plot:

  Controls ROC plot generation. Accepted values:

  `"none"`

  :   No plots are produced.

  `"all"`

  :   All selected predictors are overlaid on one plot.

  `"separate"`

  :   One plot per predictor, returned as a named list.

  Integer \\\ge 1\\

  :   Plot the top *n* predictors (by AUC) on one combined chart.

  Numeric in \\\[0, 1)\\

  :   Plot only predictors with \\AUC \ge\\ this threshold.

  `c("all", n)` or `c("separate", n)`

  :   Combine a string mode with an integer or numeric filter as
      described above.

  Default `"none"`.

- theme:

  A single string specifying the ggplot2 theme for plots. Supported
  values: `"nature"` (default), `"minimal"`, `"bw"`, `"classic"`,
  `"plos"`, `"gray"`, `"light"`, `"dark"`. Colour-blind-friendly
  palettes are applied automatically.

- ncol:

  Integer. Number of legend columns when `plot = "all"` or
  `plot = c("all", ...)`. Default `NULL` (auto).

- nrow:

  Integer. Number of legend rows when `plot = "all"` or
  `plot = c("all", ...)`. Default `NULL` (auto).

- plot_title:

  A character vector of plot titles. `NULL` (default) uses `"ROC Curve"`
  for combined plots and `"ROC Curve for <predictor>"` for separate
  plots. When `plot = "separate"`, the vector must match the number of
  predictors being plotted.

- plot_subtitle:

  A single string used as the plot subtitle. Default `NULL` (no
  subtitle).

- plot_legend:

  A single string for the legend title. Default `"AUC (95% CI)"`.

- xlab:

  A single string for the x-axis label. Default `NULL` renders
  `"1 - Specificity (False Positive Rate)"`.

- ylab:

  A single string for the y-axis label. Default `NULL` renders
  `"Sensitivity (True Positive Rate)"`.

- linewidth:

  A positive numeric value controlling the ROC curve line width. Default
  `3`.

- global_font_size:

  A positive numeric value setting the base font size for all plot text
  elements. Default `25`. Derived sizes: title \\= 1.20 \times\\,
  subtitle \\= 0.90 \times\\, axis labels / ticks \\= 1.00 \times\\.

- title_font_size:

  Override for the title font size. Default `NULL` (derived from
  `global_font_size`).

- subtitle_font_size:

  Override for the subtitle font size. Default `NULL` (derived from
  `global_font_size`).

- xlab_font_size:

  Override for the x-axis label and tick font size. Default `NULL`
  (derived from `global_font_size`).

- ylab_font_size:

  Override for the y-axis label and tick font size. Default `NULL`
  (derived from `global_font_size`).

- inplot_font_size:

  Font size of the in-plot legend text (e.g., AUC labels in combined
  plots). Default `NULL` (derived from `global_font_size`).

- min_vip:

  A single numeric value specifying the minimum VIP score for features
  to be considered. Only applied if a
  [`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md)
  object is provided in the `x` list. Default `1.0`.

- up:

  A single numeric value specifying the minimum fold change for
  up-regulated features. Only applied if a
  [`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  object is provided in the `x` list. Default is in its `params$up`.
  When `x` is `data.frame`, `tibble`, named `matrix`, the default is
  `1.5`.

- down:

  A single numeric value specifying the maximum fold change for
  down-regulated features. Only applied if a
  [`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md)
  object is provided in the `x` list. Default is in its `params$down`.
  When `x` is `data.frame`, `tibble`, named `matrix`, the default is
  `0.6`.

- p_value:

  A single numeric value specifying the maximum acceptable p-value for
  features. Only applied if a
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  object is provided in the `x` list. Default `0.05`.

- ...:

  Additional arguments forwarded to
  [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html),
  [`pROC::auc()`](https://rdrr.io/pkg/pROC/man/auc.html), or
  [`pROC::ci()`](https://rdrr.io/pkg/pROC/man/ci.html).

## Value

A named list with the following elements:

- `roc_list`:

  Named list of [`pROC::roc`](https://rdrr.io/pkg/pROC/man/roc.html)
  objects, one per analysed predictor.

- `auc_table`:

  A `data.frame` with columns `predictor`, `auc`, `ci_lower`,
  `ci_upper`, `ci_level`, `direction`, `n_cases`, and `n_controls`, plus
  `auc_interpretation` (conventional label; see Details) right after
  `ci_upper`. When a
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  object with effect sizes is supplied in the `x` list, `ES (95% CI)`
  and `ES_magnitude` columns are added after `n_controls`. When
  `compute_performance_metrics = TRUE`, additional columns are included:
  `sensitivity`, `specificity`, `ppv` (positive predictive value), `npv`
  (negative predictive value), `accuracy`, `youden_index`, and
  `optimal_threshold` (all computed at the threshold that maximizes
  Youden's index: sensitivity + specificity - 1).

- `plots`:

  A named list of `ggplot` objects (or `NULL` when `plot = "none"`). The
  element `"combined"` is present for combined-mode plots; individual
  predictor names are used as keys for separate-mode plots.

- `params`:

  A list echoing the validated call parameters for reproducibility.

## Details

**Workflow overview**  

1.  `x` is coerced to a `data.frame`; columns and rows are subset
    according to `include_cols` / `include_rows`.

2.  `y` values listed in `remove` are dropped, leaving exactly two
    groups.

3.  When `y` is categorical, the two remaining levels are mapped to `0`
    (control, as specified by `set_control`) and `1` (case). If `y`
    already encodes `0`/`1` as a character or factor, it is silently
    converted to numeric.

4.  [`pROC::roc()`](https://rdrr.io/pkg/pROC/man/roc.html) is called for
    each predictor column; AUC and, optionally, DeLong confidence
    intervals ([`pROC::ci()`](https://rdrr.io/pkg/pROC/man/ci.html)) are
    extracted.

5.  Plots are assembled with ggplot2 using a colour-blind-friendly
    palette (Okabe-Ito). AUC labels are rounded to two decimal places;
    values that round to `1.00` are progressively displayed to more
    decimal places (up to five), after which scientific notation with
    two significant figures is used.

**AUC interpretation**  
The `auc_interpretation` column of `auc_table` labels each AUC using the
conventional thresholds below. These are only a convention (a rule of
thumb), not a rule; acceptable discrimination depends on the clinical or
scientific context, and the confidence interval should be considered
alongside the point estimate.

|               |                            |
|---------------|----------------------------|
| **AUC**       | **Conventional threshold** |
| 0.5           | No discrimination          |
| 0.5 to \< 0.7 | Poor                       |
| 0.7 to \< 0.8 | Acceptable                 |
| 0.8 to \< 0.9 | Excellent                  |
| 0.9 or above  | Outstanding                |

AUC values below 0.5 (possible when `direction` is forced) are labelled
`"Worse than chance"`.

**Reference line**  
A grey dashed diagonal line representing an AUC of 0.5 (random
classifier) is always drawn.

**Legend placement**  
For combined plots (`plot = "all"`) the legend is positioned at the
bottom-right inside the plot area by default.

**pROC dependency**  
This function requires pROC (\\\ge\\ 1.18) and ggplot2 (\\\ge\\ 3.4).
Both are listed as `Imports` in the package `DESCRIPTION` file.

## References

Robin X, Turck N, Hainard A, et al. (2011). pROC: an open-source package
for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*,
12, 77.
[doi:10.1186/1471-2105-12-77](https://doi.org/10.1186/1471-2105-12-77)

DeLong ER, DeLong DM, Clarke-Pearson DL (1988). Comparing the areas
under two or more correlated receiver operating characteristic curves: a
nonparametric approach. *Biometrics*, 44(3), 837-845.
[doi:10.2307/2531595](https://doi.org/10.2307/2531595)

Youden WJ (1950). Index for rating diagnostic tests. *Cancer*, 3(1),
32-35.

Okabe M, Ito K (2008). *Color Universal Design (CUD): How to make
figures and presentations that are friendly to colorblind people*.
<https://jfly.uni-koeln.de/color/>

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.html),
[`run_pls`](https://jllcalorio.github.io/pondeR/reference/run_pls.md),
[`run_foldchange`](https://jllcalorio.github.io/pondeR/reference/run_foldchange.md),
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)

## Author

John Lennon L. Calorio

## Examples

``` r
## -----------------------------------------------------------------------
## Example 1 - numeric binary response (0/1), all predictors
## -----------------------------------------------------------------------
set.seed(42)
n  <- 120
df <- data.frame(
  group    = sample(0:1, n, replace = TRUE),
  marker_A = rnorm(n, mean = ifelse(sample(0:1, n, replace = TRUE) == 1, 1, 0)),
  marker_B = rnorm(n),
  marker_C = runif(n)
)

res <- run_auc(x = df, y = "group", plot = "none")
print(res$auc_table)
#>   predictor       auc  ci_lower  ci_upper auc_interpretation ci_level direction
#> 1  marker_A 0.5133765 0.4085016 0.6182514               Poor     0.95         >
#> 2  marker_B 0.5243593 0.4185301 0.6301886               Poor     0.95         <
#> 3  marker_C 0.5350605 0.4300641 0.6400569               Poor     0.95         >
#>   n_cases n_controls sensitivity specificity       ppv       npv  accuracy
#> 1      67         53   0.3731343   0.7924528 0.6944444 0.5000000 0.5583333
#> 2      67         53   0.7761194   0.3207547 0.5909091 0.5312500 0.5750000
#> 3      67         53   0.4029851   0.6981132 0.6279070 0.4805195 0.5333333
#>   youden_index optimal_threshold
#> 1   0.16558716        -0.2014022
#> 2   0.09687412        -0.7084071
#> 3   0.10109828         0.3504492

## -----------------------------------------------------------------------
## Example 2 - categorical response, remove one group, combined plot
## -----------------------------------------------------------------------
df2 <- data.frame(
  status   = sample(c("control", "case", "unknown"), n, replace = TRUE),
  score_1  = rnorm(n),
  score_2  = rnorm(n)
)

res2 <- run_auc(
  x           = df2,
  y           = "status",
  remove      = "unknown",
  set_control = "control",
  plot        = "none",
  compute_ci  = FALSE
)
#> Moving non-numeric columns to metadata (predictors must be numeric): status
print(res2$auc_table)
#>   predictor       auc ci_lower ci_upper auc_interpretation ci_level direction
#> 1   score_1 0.5257453       NA       NA               Poor     0.95         >
#> 2   score_2 0.5372629       NA       NA               Poor     0.95         <
#>   n_cases n_controls sensitivity specificity       ppv       npv  accuracy
#> 1      36         41   0.6388889   0.4878049 0.5227273 0.6060606 0.5584416
#> 2      36         41   0.2500000   0.8780488 0.6428571 0.5714286 0.5844156
#>   youden_index optimal_threshold
#> 1    0.1266938         0.3111415
#> 2    0.1280488         1.0426505

## -----------------------------------------------------------------------
## Example 3 - access performance metrics at optimal Youden threshold
## -----------------------------------------------------------------------
res3 <- run_auc(
  x                     = df,
  y                     = "group",
  plot                  = "none",
  compute_performance_metrics = TRUE
)
# The auc_table now includes sensitivity, specificity, PPV, NPV, accuracy,
# Youden index, and optimal threshold values
head(res3$auc_table[, c("predictor", "auc", "sensitivity", "specificity", 
                          "ppv", "npv", "youden_index")])
#>   predictor       auc sensitivity specificity       ppv       npv youden_index
#> 1  marker_A 0.5133765   0.3731343   0.7924528 0.6944444 0.5000000   0.16558716
#> 2  marker_B 0.5243593   0.7761194   0.3207547 0.5909091 0.5312500   0.09687412
#> 3  marker_C 0.5350605   0.4029851   0.6981132 0.6279070 0.4805195   0.10109828

## -----------------------------------------------------------------------
## Example 4 - column names with spaces, select top-2 by AUC
## -----------------------------------------------------------------------
df4 <- data.frame(
  `my group`  = sample(0:1, n, replace = TRUE),
  `feature 1` = rnorm(n, mean = 0.5),
  `feature 2` = rnorm(n),
  check.names = FALSE
)

res4 <- run_auc(
  x    = df4,
  y    = "my group",
  plot = "none"
)
print(res4$auc_table)
#>   predictor       auc  ci_lower  ci_upper auc_interpretation ci_level direction
#> 1 feature 1 0.5167558 0.4120382 0.6214735               Poor     0.95         >
#> 2 feature 2 0.5516756 0.4461648 0.6571863               Poor     0.95         >
#>   n_cases n_controls sensitivity specificity       ppv       npv  accuracy
#> 1      67         53   0.1343284   0.9811321 0.9000000 0.4727273 0.5083333
#> 2      67         53   0.7761194   0.3773585 0.6117647 0.5714286 0.6000000
#>   youden_index optimal_threshold
#> 1    0.1154604        -1.1017566
#> 2    0.1534779         0.6624894
```
