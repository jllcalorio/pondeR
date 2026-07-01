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
  add_inferential_pvalues = FALSE,
  test_type_continuous = c("auto", "parametric", "nonparametric"),
  test_type_categorical = c("auto", "chisq", "fisher"),
  paired = FALSE,
  n_digits_pvalues = 3,
  bold_significant_pvalues = TRUE,
  bold_significant_pvalues_at = 0.05,
  header = "Variable",
  show_n_header = TRUE,
  bold_labels = TRUE,
  italicize_levels = TRUE,
  clean_table = TRUE,
  table_name = "auto"
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

- rename_variables:

  List. Formulas of the form `list("original" ~ "new", ...)` used to
  relabel variables in the table.

- combine_categories:

  Named list of named lists. Used to combine levels of categorical
  variables. The top-level names correspond to column names in `x`. Each
  column's value is a named list where the names are the new category
  labels and the values are character vectors of original levels to
  combine. If a column is not categorical by default, it must be listed
  in `force_categorical` first. Defaults to `NULL`.

  Example format:
  ` combine_categories = list( Education = list( "Higher Ed" = c("Bachelors", "Masters"), "Schooling" = c("Elementary", "High School") ), Gender = list( "Non-Male" = c("Female", "Other") ) ) `

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

## Value

A `gt_tbl` object (via
[`gtsummary::as_gt()`](https://www.danieldsjoberg.com/gtsummary/reference/as_gt.html))
representing the descriptive (and optionally inferential) statistics
table.

## See also

[`run_diff`](https://jllcalorio.github.io/pondeR/reference/run_diff.md)
for the continuous-variable test engine.
[`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
for the categorical-variable test engine.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
library(gtsummary)
library(dplyr)

set.seed(123)
sample_df <- data.frame(
  Age       = sample(c(18:65, NA), 100, replace = TRUE,
                     prob = c(rep(1, 48), 5)),
  Gender    = sample(c("Male", "Female", NA), 100, replace = TRUE,
                     prob = c(0.48, 0.48, 0.04)),
  Education = sample(c("Elementary", "High School",
                        "Bachelors", "Masters", NA), 100,
                     replace = TRUE,
                     prob = c(0.40, 0.30, 0.15, 0.10, 0.05)),
  Income    = sample(c(seq(5000, 40000, by = 5000), NA), 100,
                     replace = TRUE, prob = c(rep(1, 8), 2)),
  Smoker    = sample(c(1, 0, NA), 100, replace = TRUE,
                     prob = c(0.25, 0.65, 0.10)),
  Region    = sample(c("North", "South", "East", "West", NA), 100,
                     replace = TRUE,
                     prob = c(0.24, 0.24, 0.24, 0.24, .04))
)

# Basic descriptive statistics (no grouping)
sample_df |> run_summarytable()

# Descriptive statistics split by Gender
sample_df |>
  run_summarytable(
    split_by        = "Gender",
    split_by_header = "Participant Gender",
    rename_variables = list(Age    ~ "Age (Years)",
                            Income ~ "Annual Income",
                            Smoker ~ "Smoking Status")
  )

# Median and IQR for continuous variables
sample_df |>
  run_summarytable(
    continuous_statistics  = "medianIQR",
    categorical_statistics = "n",
    n_digits_continuous    = c(0, 0),
    n_digits_categorical   = c(0, 0)
  )

# Inferential p-values using run_diff (continuous) and run_assoc
# (categorical), auto test selection
sample_df |>
  run_summarytable(
    split_by                = "Gender",
    add_inferential_pvalues = TRUE
  )

# Force non-parametric tests for continuous variables
sample_df |>
  run_summarytable(
    split_by                = "Gender",
    add_inferential_pvalues = TRUE,
    test_type_continuous    = "nonparametric",
    test_type_categorical   = "fisher"
  )

# Paired inferential p-values (descriptive body remains unpaired)
sample_df |>
  run_summarytable(
    split_by                = "Gender",
    add_inferential_pvalues = TRUE,
    paired                  = TRUE
  )
} # }
```
