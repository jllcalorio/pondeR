# Perform Binary Logistic Regression (Standard or Firth-Corrected) with Iterative Feature Combinations

A robust wrapper for [`glm`](https://rdrr.io/r/stats/glm.html) (standard
logistic regression) and
[`logistf`](https://rdrr.io/pkg/logistf/man/logistf.html) (Firth's
bias-reduced penalized-likelihood logistic regression) that performs
binary logistic regression. It automatically handles column names with
spaces, computes advanced model fit metrics, maps specific reference
categories, removes specified categories from predictor variables, and
optionally iterates through all possible combinations of independent
variables. When `apply_firth = TRUE`, Firth's correction is applied
selectively — only to models where at least one predictor-outcome
cross-tabulation contains a zero cell, which is the primary condition
warranting penalized likelihood. A column in the output indicates
whether Firth's correction was applied to each model. Output is
formatted as a clean data frame.

## Usage

``` r
run_logreg(
  x,
  y,
  confounders = NULL,
  indep = NULL,
  ref,
  ref_levels = NULL,
  remove = NULL,
  exclude = NULL,
  iterate = FALSE,
  add_vif = TRUE,
  apply_firth = TRUE,
  force_apply_firth = FALSE
)
```

## Arguments

- x:

  A data frame, matrix, or tibble containing the dataset.

- y:

  A string specifying the valid, categorical column name in `x` to be
  used as the dependent variable.

- confounders:

  A character vector of valid column names in `x` to be used as
  confounding factors. Default is `NULL`.

- indep:

  A character vector of valid column names in `x` to be used as
  independent variables. Default is `NULL`.

- ref:

  A string specifying the valid category of `y` to be used as the
  reference category.

- ref_levels:

  A list of two-sided formulas mapping categorical predictor column
  names to their desired reference categories. The left-hand side must
  be a valid column name of `x` (quoted as a string), and the right-hand
  side must be a valid category of that column (e.g.,
  `list("Sex" ~ "Female", "District" ~ "North")`). Default is `NULL`.

- remove:

  A list of two-sided formulas specifying categories to drop from
  predictor variables before analysis. The left-hand side must be a
  valid column name of `x` (quoted as a string), and the right-hand side
  must be a valid category of that column to be removed. Multiple
  mappings may target the same column to remove multiple categories
  (e.g., `list("Sex" ~ "Unknown", "District" ~ "Other")`). Default is
  `NULL`.

- exclude:

  A character vector of categories in `y` to drop before analysis.
  Default is `NULL`.

- iterate:

  Logical; if `FALSE` (default), combines `confounders` and `indep` into
  a single model. If `TRUE`, fits separate models for all possible
  combinations of `indep` added to `confounders`.

- add_vif:

  Logical; if `TRUE` (default), adds Variance Inflation Factor (VIF)
  values to the output. VIFs are computed via
  [`car::vif()`](https://rdrr.io/pkg/car/man/vif.html), which returns
  GVIF for categorical predictors and standard VIF for continuous
  predictors, consistent with jamovi's output. Note that VIF is computed
  from a standard [`glm()`](https://rdrr.io/r/stats/glm.html) fit even
  when Firth's correction is applied, as
  [`car::vif()`](https://rdrr.io/pkg/car/man/vif.html) does not support
  `logistf` objects.

- apply_firth:

  Logical; if `TRUE` (default), Firth's bias-reduced
  penalized-likelihood logistic regression
  ([`logistf::logistf()`](https://rdrr.io/pkg/logistf/man/logistf.html))
  is applied *selectively* — only to models where at least one predictor
  variable has a zero cell count in its cross-tabulation with the
  outcome. Models with no zero cells are fitted using standard
  [`glm()`](https://rdrr.io/r/stats/glm.html) regardless of this
  setting. A logical column `Firth Corrected` is added to the `Metrics`
  table (and attached metrics for `iterate = FALSE`) indicating whether
  Firth's correction was applied to each model.

- force_apply_firth:

  Logical; if `TRUE` (default `FALSE`), overrides the dynamic zero-cell
  detection of `apply_firth` and applies Firth's correction to *all*
  models regardless of whether zero cells are present. Ignored if
  `apply_firth = FALSE`.

## Value

If `iterate = FALSE`, returns a data frame containing Predictor, Log
Odds, Std Error, p-value, Significance, OR, 95\\ `Firth Corrected`) are
attached as `attr(result, "model_metrics")`.

If `iterate = TRUE`, returns a named list containing:

- Tables:

  A named list of result data frames for each model combination.

- Metrics:

  A data frame comparing model metrics (N, Deviance, AIC, BIC, McFadden
  \\R^2\\, Cox-Snell \\R^2\\, Nagelkerke \\R^2\\, Tjur \\R^2\\, Has
  Significant Predictor, and Firth Corrected) across all combinations.

- filtered_Tables:

  A subset of `Tables` excluding only models where *all* non-intercept
  `95% CI Low` values or *all* non-intercept `95% CI High` values are
  `Inf`.

- p_filtered_Tables:

  A subset of `filtered_Tables` retaining only models where at least one
  non-intercept `p-value` is less than .05.

- vif_filtered_Tables:

  A subset of `p_filtered_Tables` retaining only models where all
  non-`NA` VIF values are less than 5, ensuring no violation of
  multicollinearity.

## Details

When `iterate = TRUE`, the function generates all possible combinations
(the power set) of the `indep` vector. For example, if
`indep = c("A", "B", "C")`, the combinations evaluated alongside
`confounders` will be: `A`, `B`, `C`, `A+B`, `A+C`, `B+C`, and `A+B+C`.

**Firth's correction (`apply_firth`):** Firth's penalized-likelihood
method
([`logistf::logistf()`](https://rdrr.io/pkg/logistf/man/logistf.html))
addresses the monotone likelihood problem caused by complete or
quasi-complete separation, which arises when a predictor level perfectly
predicts one outcome group — typically evidenced by a zero cell in the
cross-tabulation. When `apply_firth = TRUE`, the function checks each
model's predictors against the outcome for zero cells. If any are found,
`logistf()` is used in place of
[`glm()`](https://rdrr.io/r/stats/glm.html) for that model. If
`force_apply_firth = TRUE`, this check is skipped and all models are
fitted with `logistf()`. The `Firth Corrected` column in the output
records which models received the correction. Note that Firth-corrected
models report profile likelihood confidence intervals rather than Wald
intervals, consistent with `logistf`'s default behavior.

**Category removal (`remove`):** Rows where the specified predictor
variable equals the given category are dropped from `x` before any model
is fitted. This is applied after `exclude` but before factor releveling.
Multiple formulas targeting the same column are each applied in
sequence, allowing several categories to be removed from one variable.

**Aliased coefficients:** When standard
[`glm()`](https://rdrr.io/r/stats/glm.html) cannot estimate a
coefficient due to perfect or quasi-complete separation, it sets that
coefficient to `NA`. [`summary()`](https://rdrr.io/r/base/summary.html)
silently omits these aliased rows while
[`confint.default()`](https://rdrr.io/r/stats/confint.html) retains
them, causing a row-count mismatch. `run_logreg()` detects this
automatically, aligns both objects to their shared row names, and emits
a warning naming every dropped aliased term. This situation is largely
avoided when `apply_firth = TRUE` since such models will be routed to
`logistf()` instead.

**VIF computation:** This function uses
[`car::vif()`](https://rdrr.io/pkg/car/man/vif.html) to match jamovi's
output. For categorical predictors with more than one degree of freedom,
[`car::vif()`](https://rdrr.io/pkg/car/man/vif.html) returns the
Generalized VIF (GVIF) and its degree-freedom-adjusted form
\\GVIF^{1/(2\cdot df)}\\. The value stored in the `VIF` column is the
raw GVIF (or standard VIF for continuous/binary predictors). Because
[`car::vif()`](https://rdrr.io/pkg/car/man/vif.html) does not support
`logistf` objects, VIF for Firth-corrected models is computed from a
parallel standard [`glm()`](https://rdrr.io/r/stats/glm.html) fit on the
same formula and data.

**Significance codes:** The `Significance` column uses the conventional
scheme: `***` if \\p \< .001\\, `**` if \\p \< .01\\, `*` if \\p \<
.05\\, and blank otherwise.

## References

Fox J, Weisberg S (2019). *An R Companion to Applied Regression*, 3rd
ed. Sage.

Firth D (1993). "Bias reduction of maximum likelihood estimates."
*Biometrika*, 80(1), 27–38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27)

Heinze G, Schemper M (2002). "A solution to the problem of separation in
logistic regression." *Statistics in Medicine*, 21(16), 2409–2419.
[doi:10.1002/sim.1047](https://doi.org/10.1002/sim.1047)

Heinze G, Ploner M, Jiricka L, Steuer A (2024). *logistf: Firth's
Bias-Reduced Logistic Regression*. R package.
<https://CRAN.R-project.org/package=logistf>

Lüdecke D, Ben-Shachar MS, Patil I, Waggoner P, Makowski D (2021).
"performance: An R Package for Assessment, Comparison and Testing of
Statistical Models." *Journal of Open Source Software*, 6(60), 3139.
[doi:10.21105/joss.03139](https://doi.org/10.21105/joss.03139)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Example 1: Single model with Firth's correction applied dynamically ---
# Firth's correction is applied automatically only if zero cells are detected.
result <- run_logreg(
  x          = my_data,
  y          = "Outcome",
  confounders = c("Age", "Sex"),
  indep      = c("Biomarker_A", "Biomarker_B"),
  ref        = "Control",
  ref_levels = list("Sex" ~ "Female"),
  apply_firth = TRUE
)
result
attr(result, "model_metrics")  # includes Firth Corrected column

# --- Example 2: Iterative models with forced Firth and VIF filtering -------
# All models receive Firth's correction regardless of zero-cell status.
results <- run_logreg(
  x                = my_data,
  y                = "Outcome",
  confounders      = c("Age", "Sex", "Site"),
  indep            = c("Biomarker_A", "Biomarker_B", "Biomarker_C"),
  ref              = "Control",
  ref_levels       = list("Sex" ~ "Female", "Site" ~ "Main"),
  exclude          = "Indeterminate",
  iterate          = TRUE,
  apply_firth      = TRUE,
  force_apply_firth = TRUE
)
results$Metrics
results$vif_filtered_Tables  # models with p < .05 and all VIF < 5
} # }
```
