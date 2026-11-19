# Performs Descriptive Statistics on a Dataframe

Generates a comprehensive descriptive statistics table from a given
dataframe, optionally stratified by a categorical variable. Provides
customisable summaries for both continuous and categorical variables,
including control over rounding, missing data reporting, labelling,
formatting, and table aesthetics. The output is a publication-ready
`gtsummary` / `gt` table.

When `add_inferential_pvalues = TRUE`, p-values are computed *outside*
`gtsummary` using
[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
(for continuous variables) and
[`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
(for categorical variables), then injected into the table. Each p-value
cell carries a letter superscript identifying the exact test used; a
legend is appended as table footnotes. Variables for which a p-value
could not be computed are shown as `---†` with a corresponding footnote.

## Usage

``` r
run_summarytable(
  x,
  summarize_what = NULL,
  split_by = NULL,
  filter = NULL,
  split_by_header = NULL,
  strata_by = NULL,
  subgroup = NULL,
  rename_variables = NULL,
  combine_categories = NULL,
  continuous_statistics = "meanSD",
  categorical_statistics = "n_percent",
  force_continuous = NULL,
  force_categorical = NULL,
  n_digits_continuous = c(2, 2),
  n_digits_categorical = c(0, 2),
  zero_as_exp = TRUE,
  display_missing = "ifany",
  missing_text = "No data/missing",
  missing_stat = "n_percent",
  include_missing_in_splits = FALSE,
  sort_categorical_variables_by = "alphanumeric",
  calc_percent_by = "column",
  calc_col_percent_using = "n_valid_in_column",
  include_missing_in_denom = FALSE,
  add_inferential_pvalues = FALSE,
  test_type_continuous = c("auto", "parametric", "nonparametric"),
  test_type_categorical = c("auto", "chisq", "fisher"),
  paired = FALSE,
  n_digits_pvalues = 3,
  bold_significant_pvalues = TRUE,
  bold_significant_pvalues_at = 0.05,
  mark_posthoc_significance = TRUE,
  p_adjust_method = "bonferroni",
  header = "Variable",
  show_n_header = TRUE,
  bold_labels = TRUE,
  italicize_levels = TRUE,
  clean_table = TRUE,
  table_name = "auto",
  .in_subgroup_recursion = FALSE
)
```

## Arguments

- x:

  Dataframe. The input dataframe.

- summarize_what:

  Character vector. Column names to include in the analysis. When `NULL`
  (default), all columns except those named in `split_by` and
  `strata_by` are included.

- split_by:

  String. A column name in `x` (must be categorical with at least 2
  unique levels) by which to stratify the table. Required when
  `add_inferential_pvalues = TRUE`. Default is `NULL`.

- filter:

  Optional character vector of levels in `split_by` to exclude from the
  analysis. If `x` is a `run_DIpreprocess` object, defaults to the
  identified QC types.

- split_by_header:

  String. Display header for the `split_by` variable. Defaults to the
  value of `split_by`.

- strata_by:

  String. A column name in `x` (must be categorical with at least 2
  unique levels) by which to stratify *before* `split_by`. Cannot be
  combined with `add_inferential_pvalues = TRUE`. Default is `NULL`.

- subgroup:

  String or `NULL`. A column name in `x` (must be categorical) for
  subgroup analysis. When provided, the function computes separate
  summary tables for each unique value in the column and returns a list
  of tables. Cannot be the same column as `split_by` or `strata_by`.
  Default is `NULL`.

- rename_variables:

  List. Formulas of the form `list("original" ~ "new", ...)` used to
  relabel variables in the table.

- combine_categories:

  Named list of named lists. Used to combine levels of categorical
  variables. The top-level names correspond to column names in `x`. Each
  column's value is a named list where the names are the new category
  labels and the values are character vectors of original levels to
  combine. If a column is not categorical by default, it must be listed
  in `force_categorical` first. Defaults to `NULL`. See examples for
  usage.

- continuous_statistics:

  String. Summary statistic(s) for continuous variables. One of:
  `"meanSD"`, `"meanSD2"`, `"medianIQR"`, `"mean"`, `"sd"`, `"median"`,
  `"p0"`, `"p25"`, `"p50"`, `"p75"`, `"p100"`, `"IQR"`. Default is
  `"meanSD"`.

- categorical_statistics:

  String. Summary statistic(s) for categorical variables. One of:
  `"n_percent"`, `"n"`, `"percent"`. Default is `"n_percent"`.

- force_continuous:

  Character vector. Column names to treat as continuous regardless of
  their storage type. The same coercion is applied before passing data
  to
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md).

- force_categorical:

  Character vector. Column names to treat as categorical regardless of
  their storage type. The same coercion is applied before passing data
  to
  [`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md).

- n_digits_continuous:

  Numeric vector `c(d1, d2)` or `"dynamic"`. Number of decimal places
  for the first and subsequent numeric tokens in continuous summary
  cells (e.g., mean and SD). Default is `c(2, 2)`.

- n_digits_categorical:

  Numeric vector `c(d1, d2)` or `"dynamic"`. Decimal places for count
  and percentage in categorical summary cells. Default is `c(0, 2)`.

- zero_as_exp:

  Logical. When `TRUE` (default), continuous summary tokens that round
  to `0` at the requested precision but are genuinely non-zero are
  displayed in scientific notation (e.g., `1.00e-3`).

- display_missing:

  String. How to report missing values: `"ifany"` (default), `"no"`, or
  `"always"`.

- missing_text:

  String. Label for missing-value rows. Default is `"No data/missing"`.

- missing_stat:

  String. Format for missing-value counts: one of `"n_percent"`
  (default), `"n"`, `"percent"`.

- include_missing_in_splits:

  Logical. If `TRUE`, missing values in `split_by` / `strata_by` columns
  are kept as an explicit level. Defaults to `FALSE`.

- sort_categorical_variables_by:

  String. Level ordering for categorical variables: `"alphanumeric"`
  (default) or `"frequency"`.

- calc_percent_by:

  String. Percentage base for categorical variables: `"column"`
  (default), `"row"`, `"cell"`, or `"total"`. `"total"` calculates
  percentages based on the total number of samples, study participants,
  or simply number of rows in `x`.

- calc_col_percent_using:

  String. Denominator when `calc_percent_by = "column"`:
  `"n_valid_in_column"` (default) or `"n_in_column"`.

- include_missing_in_denom:

  Logical. Only applies when `calc_percent_by = "column"`. Controls
  whether missing observations are included in the percentage
  denominator for categorical levels. When `FALSE` (default),
  non-missing level percentages are calculated using only the valid
  (non-missing) observations as the denominator, while the "No
  data/missing" row uses the total group N. When `TRUE`, all level
  percentages - including the missing row - are calculated using the
  total group N (i.e., missing observations count toward the
  denominator). For example, with 2 Normocardia, 2 Tachycardia, and 41
  missing (N = 45): `FALSE` gives 50\\, 50\\, 91.11\\; `TRUE` gives
  4.44\\, 4.44\\, 91.11\\. Defaults to `FALSE`.

- add_inferential_pvalues:

  Logical. If `TRUE`, p-values are computed via
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  (continuous variables) and
  [`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  (categorical variables) and injected into the table. Requires
  `split_by` to be non-`NULL` and `strata_by` to be `NULL`. Default is
  `FALSE`.

- test_type_continuous:

  String. Test strategy forwarded to
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  when `add_inferential_pvalues = TRUE`. One of `"auto"` (default),
  `"parametric"`, or `"nonparametric"`. See
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  for details.

- test_type_categorical:

  String. Test strategy forwarded to
  [`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  when `add_inferential_pvalues = TRUE`. One of `"auto"` (default),
  `"chisq"`, or `"fisher"`. See
  [`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  for details.

- paired:

  Logical. If `TRUE`, paired tests are used when computing p-values.
  Forwarded to both
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  and
  [`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md).
  Note that the descriptive statistics shown in the table body are
  always unpaired summaries. Default is `FALSE`.

- n_digits_pvalues:

  Numeric or `NULL`. Decimal places for displayed p-values. Default is
  `3`. Set to `NULL` for automatic formatting (up to 3 d.p. when \\p \<
  0.10\\, otherwise 1 d.p.).

- bold_significant_pvalues:

  Logical. If `TRUE` (default), p-values below
  `bold_significant_pvalues_at` are bolded.

- bold_significant_pvalues_at:

  Numeric. Significance threshold for bolding. Default is `0.05`.

- mark_posthoc_significance:

  Logical. Only takes effect when `add_inferential_pvalues = TRUE`. When
  `TRUE`, performs post-hoc testing for variables whose omnibus p-value
  is below `bold_significant_pvalues_at` and marks the significant
  summary cells with numeric superscripts (numbering starts at 2 for
  each variable, as 1 is the summary-statistic footnote). For continuous
  variables with \\\>2\\ groups (ANOVA / Welch's ANOVA /
  Kruskal-Wallis), pairwise post-hoc comparisons (Tukey HSD /
  Games-Howell / Dunn's test, as already computed by
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md))
  mark the two group columns involved in each significant pair. Not
  applicable to two-group comparisons (t-test / Mann-Whitney), which
  have no post-hoc step. For categorical variables (Chi-square /
  Fisher's exact / McNemar's), adjusted standardized residuals (Z-tests
  for proportions) are computed per cell of the level-by-group table;
  cells with an adjusted p-value below `bold_significant_pvalues_at` are
  individually marked. This includes binary variables that `gtsummary`
  collapses to a single dichotomous row (e.g. Yes/No): the displayed
  level is identified from its printed count, and only the ambiguous
  case of tied counts across levels is left unmarked. A footnote block
  lists each numbered comparison/cell per variable. For two-group
  categorical variables (2 x 2 tables) cell-level post-hoc is skipped,
  since the cells mirror the omnibus test; a one-line direction note
  (which level is more common in the second group, with column
  percentages) is added to the footnote instead. For other tables with a
  two-level `split_by`, the cell p-values are adjusted over one column
  only (the two columns mirror each other, so there are only as many
  independent tests as levels). When the omnibus test is significant but
  no pair/cell is significant after adjustment, a footnote says so. Has
  no effect when `add_inferential_pvalues = FALSE`. Default is `TRUE`.

- p_adjust_method:

  String. P-value adjustment method used for post-hoc comparisons when
  `mark_posthoc_significance = TRUE`. Forwarded to
  [`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
  and
  [`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
  (affects Dunn's test pairwise p-values for Kruskal-Wallis follow-up),
  and used directly to adjust the per-cell
  adjusted-standardized-residual p-values for categorical variables. One
  of [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html)
  (`"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"` (default), `"BH"`,
  `"BY"`, `"fdr"`, `"none"`).

- header:

  String. Header for the variable-label column. Default is `"Variable"`.

- show_n_header:

  Logical. When `TRUE` (default), appends the total sample size to the
  `header` (e.g., "Variable (N = 100)").

- bold_labels:

  Logical. If `TRUE` (default), variable labels are bolded.

- italicize_levels:

  Logical. If `TRUE` (default), categorical level labels are italicised.

- clean_table:

  Logical. If `TRUE` (default), cells containing a zero count (e.g.,
  `0 (0.00%)`) are replaced with empty strings.

- table_name:

  String or `NULL`. Table title. `"auto"` (default) generates a title
  automatically. `NULL` suppresses the title.

- .in_subgroup_recursion:

  Internal use only. Logical flag set to `TRUE` when
  `run_summarytable()` calls itself recursively to build subgroup
  tables; not intended to be set by users.

## Value

A `gt_tbl` object (via
[`gtsummary::as_gt()`](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html))
representing the descriptive (and optionally inferential) statistics
table.

## References

Sjoberg DD, Whiting K, Curry M, Lavery JA, Larmarange J (2021).
Reproducible summary tables with the gtsummary package. *The R Journal*,
13(1), 570-580.
[doi:10.32614/RJ-2021-053](https://doi.org/10.32614/RJ-2021-053)

Welch BL (1947). The generalization of 'Student's' problem when several
different population variances are involved. *Biometrika*, 34(1-2),
28-35.
[doi:10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28)

Mann HB, Whitney DR (1947). On a test of whether one of two random
variables is stochastically larger than the other. *The Annals of
Mathematical Statistics*, 18(1), 50-60.
[doi:10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)

Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance
analysis. *Journal of the American Statistical Association*, 47(260),
583-621.
[doi:10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)

Tukey JW (1949). Comparing individual means in the analysis of variance.
*Biometrics*, 5(2), 99-114.
[doi:10.2307/3001913](https://doi.org/10.2307/3001913)

Games PA, Howell JF (1976). Pairwise multiple comparison procedures with
unequal n's and/or variances: a Monte Carlo study. *Journal of
Educational Statistics*, 1(2), 113-125.
[doi:10.3102/10769986001002113](https://doi.org/10.3102/10769986001002113)

Dunn OJ (1964). Multiple comparisons using rank sums. *Technometrics*,
6(3), 241-252.
[doi:10.1080/00401706.1964.10490181](https://doi.org/10.1080/00401706.1964.10490181)

Pearson K (1900). On the criterion that a given system of deviations
from the probable in the case of a correlated system of variables is
such that it can be reasonably supposed to have arisen from random
sampling. *Philosophical Magazine*, 50(302), 157-175.
[doi:10.1080/14786440009463897](https://doi.org/10.1080/14786440009463897)

Fisher RA (1922). On the interpretation of chi-squared from contingency
tables, and the calculation of P. *Journal of the Royal Statistical
Society*, 85(1), 87-94.
[doi:10.2307/2340521](https://doi.org/10.2307/2340521)

McNemar Q (1947). Note on the sampling error of the difference between
correlated proportions or percentages. *Psychometrika*, 12(2), 153-157.
[doi:10.1007/BF02295996](https://doi.org/10.1007/BF02295996)

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *Journal of the
Royal Statistical Society: Series B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## See also

[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
for the continuous-variable test engine.
[`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
for the categorical-variable test engine.

## Author

John Lennon L. Calorio

## Examples

``` r
# Quick example
run_summarytable(iris, split_by = "Species")


  


Descriptive Statistics Stratified by Species Column
```
