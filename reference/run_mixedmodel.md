# Perform Linear Mixed-Effects Model Analysis

Fits Linear Mixed-Effects Models (LMEMs) across all features (columns)
in a data matrix or data frame, using metadata for model covariates and
grouping structure. Results are structured to support downstream
diagnostic and effect visualizations including residuals vs. fitted, Q-Q
plots, observed vs. predicted, effect plots, and random effects
caterpillar plots.

## Usage

``` r
run_mixedmodel(
  x,
  metadata,
  formula,
  subject_col = NULL,
  group = NULL,
  p_adjust_method = "BH",
  exclude = NULL,
  verbose = FALSE
)
```

## Arguments

- x:

  A data frame, tibble, or matrix of numeric features to be modeled
  (e.g., metabolites, genes). Rows must correspond to observations;
  columns to features. Special characters in column names are supported.

- metadata:

  A data frame or tibble with the same number of rows as `x`, containing
  covariates, grouping variables, and subject identifiers used in the
  model formula.

- formula:

  A character string specifying the LMEM formula. The response variable
  should be written as `.feature` as a placeholder, which will be
  substituted for each column in `x` during iteration (e.g.,
  `".feature ~ time + (1 | subject_id)"`). All variables referenced in
  the formula (other than `.feature`) must be present in `metadata`.

- subject_col:

  A single character string specifying the column name in `metadata` for
  subject identifiers used as random effects. Defaults to `NULL`, which
  triggers automatic detection by searching for columns named `"ID"`,
  `"id"`, `"SubjectID"`, or `"subject_id"` in that order. A warning is
  issued if no match is found and `subject_col` remains `NULL`.

- group:

  A single character string specifying a column name in `metadata` that
  encodes group membership (e.g., treatment arm, timepoint). Used for
  structuring summaries and downstream visualizations. Defaults to
  `NULL`.

- p_adjust_method:

  A character string specifying the p-value adjustment method for the
  combined fixed effects summary. Passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). One of `"holm"`,
  `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, or
  `"none"`. Default is `"BH"`.

- exclude:

  A character vector of column names present in `x`, `metadata`, or
  both, to be removed before any analysis. Defaults to `NULL` (no
  exclusions).

- verbose:

  Logical. If `TRUE`, emits a message for each feature as it is
  processed. Default is `FALSE`.

## Value

A named list of class `c("run_mixedmodel", "list")` with the following
top-level components:

- `<feature_name>`:

  One entry per feature in `x`. Each is a list with:

  `status`

  :   Character string - `"success"` or failure reason.

  `model`

  :   Fitted `lmerModLmerTest` object or `NULL`.

  `fixed_effects`

  :   Data frame of fixed effects.

  `random_effects`

  :   From [`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html).

  `n_subjects`

  :   Integer.

  `n_obs`

  :   Integer.

  `convergence_warnings`

  :   Character vector or `NULL`.

  `plot_data`

  :   A list of plotting data:

      `resid_vs_fitted`

      :   Columns: `fitted`, `residuals`.

      `qq`

      :   Columns: `theoretical`, `sample`.

      `obs_vs_pred`

      :   Columns: `observed`, `predicted`.

      `effects`

      :   Effect estimates with CIs.

      `ranef_caterpillar`

      :   Random effects summaries.

- `combined_fixed_effects`:

  Data frame combining fixed effects from all successful models.

- `parameters`:

  List of input parameters for reproducibility.

- `processing_summary`:

  Summary with `n_features`, `n_success`, `n_failed`, and
  `success_rate`.

## Details

The function iterates over all columns in `x` (after exclusions),
binding each as the response variable into a combined analysis data
frame alongside `metadata`. The model formula is constructed by
substituting the literal token `.feature` with the backtick-quoted
column name of the current feature. This allows formulas with arbitrary
fixed and random effects structures to be reused across all features
without modification.

Models are fit using
[`lmerTest::lmer`](https://rdrr.io/pkg/lmerTest/man/lmer.html) with REML
estimation, which provides Satterthwaite-approximated degrees of freedom
and p-values for fixed effects.

For each successfully fitted model, the following diagnostic data are
extracted and stored to support downstream visualization:

- **Residuals vs. Fitted**: raw residuals and fitted values.

- **Q-Q plot**: standardized residuals and theoretical quantiles.

- **Observed vs. Predicted**: response vector alongside fitted values.

- **Effect plots**: fixed effect estimates with 95\\

- **Random effects (caterpillar)**: conditional modes and conditional
  standard deviations per random effect grouping factor.

Adjusted p-values in `Combined_Fixed_Effects_Summary` are computed
separately within each unique fixed effect term across all features,
using the method specified in `p_adjust_method`.

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
n_subj  <- 20
n_times <- 2
n       <- n_subj * n_times

metadata <- data.frame(
  subject_id = factor(rep(seq_len(n_subj), each = n_times)),
  time       = rep(c("pre", "post"), n_subj),
  age        = rep(round(runif(n_subj, 25, 60)), each = n_times),
  stringsAsFactors = FALSE
)

x <- data.frame(
  glucose   = rnorm(n, mean = ifelse(metadata$time == "post", 5.5, 5.0), sd = 0.5),
  insulin   = rnorm(n, mean = ifelse(metadata$time == "post", 90,  80 ), sd = 10),
  `PA O-28:1` = rnorm(n, mean = 1.2, sd = 0.3),
  check.names = FALSE
)

result <- run_mixedmodel(
  x        = x,
  metadata = metadata,
  formula  = ".feature ~ time + age + (1 | subject_id)",
  group    = "time"
)

print(result)
summary(result)
head(result$combined_fixed_effects)
} # }
```
