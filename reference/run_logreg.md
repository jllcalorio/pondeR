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

When `crude = TRUE`, a separate simple logistic regression is fitted for
each variable in `indep` against `y`, producing unadjusted (crude) odds
ratios. Results are stacked into a single tidy data frame. The
`confounders` and `iterate` arguments are ignored in this mode.

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
  crude = FALSE,
  add_vif = TRUE,
  apply_firth = TRUE,
  force_apply_firth = FALSE,
  firth_control = NULL,
  firth_plcontrol = NULL,
  stepwise = FALSE,
  step_n = 1L,
  verbose = FALSE
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
  confounding factors. Ignored when `crude = TRUE` (a warning is
  emitted). Default is `NULL`.

- indep:

  A character vector of valid column names in `x` to be used as
  independent variables. When `crude = TRUE`, each variable is modelled
  individually against `y`. Default is `NULL`.

- ref:

  A string specifying the valid category of `y` to be used as the
  reference category.

- ref_levels:

  A list of two-sided formulas mapping categorical predictor column
  names to their desired reference categories. The left-hand side must
  be a valid column name of `x` (quoted as a string), and the right-hand
  side must be a valid category of that column (e.g.,
  `list("Sex" ~ "Female", "District" ~ "North")`). Applied in both
  standard and crude modes. Default is `NULL`.

- remove:

  A list of two-sided formulas specifying categories to drop from
  predictor variables before analysis. The left-hand side must be a
  valid column name of `x` (quoted as a string), and the right-hand side
  must be a valid category of that column to be removed. Multiple
  mappings may target the same column to remove multiple categories
  (e.g., `list("Sex" ~ "Unknown", "District" ~ "Other")`). Applied in
  both standard and crude modes. Default is `NULL`.

- exclude:

  A character vector of categories in `y` to drop before analysis.
  Applied in both standard and crude modes. Default is `NULL`.

- iterate:

  Logical; if `FALSE` (default), combines `confounders` and `indep` into
  a single model. If `TRUE`, fits separate models for all possible
  combinations of `indep` added to `confounders`. Ignored when
  `crude = TRUE`.

- crude:

  Logical; if `TRUE`, fits a separate simple (unadjusted) logistic
  regression for each variable in `indep` against `y`, producing crude
  odds ratios. Results are stacked into a single tidy data frame with
  the same column structure as the standard output. `confounders` and
  `iterate` are ignored in this mode (a warning is emitted if
  `confounders` is non-`NULL`). All other parameters (`ref`,
  `ref_levels`, `remove`, `exclude`, `add_vif`, `apply_firth`,
  `force_apply_firth`) remain applicable. Default is `FALSE`.

- add_vif:

  Logical; if `TRUE` (default), adds Variance Inflation Factor (VIF)
  values to the output. VIFs are computed via
  [`car::vif()`](https://rdrr.io/pkg/car/man/vif.html), which returns
  GVIF for categorical predictors and standard VIF for continuous
  predictors, consistent with jamovi's output. Note that VIF is computed
  from a standard [`glm()`](https://rdrr.io/r/stats/glm.html) fit even
  when Firth's correction is applied, as
  [`car::vif()`](https://rdrr.io/pkg/car/man/vif.html) does not support
  `logistf` objects. In crude mode, VIF is not meaningful for simple
  (one-predictor) models and will be `NA` throughout.

- apply_firth:

  Logical; if `TRUE` (default), Firth's bias-reduced
  penalized-likelihood logistic regression
  ([`logistf::logistf()`](https://rdrr.io/pkg/logistf/man/logistf.html))
  is applied *selectively* — only to models where at least one predictor
  variable has a zero cell count in its cross-tabulation with the
  outcome. Models with no zero cells are fitted using standard
  [`glm()`](https://rdrr.io/r/stats/glm.html) regardless of this
  setting. A logical column `Firth Corrected` is added to the output
  indicating whether Firth's correction was applied to each model.
  Applicable in both standard and crude modes.

- force_apply_firth:

  Logical; if `TRUE` (default `FALSE`), overrides the dynamic zero-cell
  detection of `apply_firth` and applies Firth's correction to *all*
  models regardless of whether zero cells are present. Ignored if
  `apply_firth = FALSE`. Applicable in both standard and crude modes.

- firth_control:

  Optional
  [`logistf.control`](https://rdrr.io/pkg/logistf/man/logistf.control.html)
  object passed through to `logistf()`'s `control` argument, governing
  the iterative fitting algorithm (e.g.
  `logistf::logistf.control(maxit = 200)`). Use this when a
  Firth-corrected model reports "Maximum number of iterations ...
  exceeded" for the main fit. Ignored for
  [`glm()`](https://rdrr.io/r/stats/glm.html) models and when no model
  in the run requires Firth's correction. Default is `NULL` (logistf's
  own default control).

- firth_plcontrol:

  Optional
  [`logistf.control`](https://rdrr.io/pkg/logistf/man/logistf.control.html)
  object (despite the name, built the same way as `firth_control`)
  passed through to `logistf()`'s `plcontrol` argument, governing
  convergence of the profile-likelihood confidence intervals. Use this
  when a Firth-corrected model reports "Nonconverged PL confidence
  limits". Default is `NULL` (logistf's own default control).

- stepwise:

  Logical; if `TRUE` (default `FALSE`), performs VIF-based backward
  elimination on `indep` before the final model is fitted. At each step,
  a VIF is computed per variable (via a parallel
  [`glm()`](https://rdrr.io/r/stats/glm.html) fit, one value per
  variable rather than per dummy level of a factor) for the current
  `confounders + indep` set; if any `indep` variable has `VIF >= 5`, the
  `step_n` variable(s) with the highest VIF are dropped and the process
  repeats. It stops once no remaining `indep` variable has `VIF >= 5`,
  or when only one `indep` variable remains (VIF is undefined for a
  single predictor). `confounders` are never removed. Mutually exclusive
  with `iterate` (an error is raised if both are `TRUE`) and ignored
  (with a warning) when `crude = TRUE`. Forces `add_vif = TRUE` if it
  was set to `FALSE`.

- step_n:

  A positive integer (default `1`), at most
  `length(indep) + length(confounders)`, giving how many high-VIF
  `indep` variables are removed per step when `stepwise = TRUE`. Ignored
  otherwise.

- verbose:

  Logical; if `TRUE` (default `FALSE`), prints progress messages via
  [`message()`](https://rdrr.io/r/base/message.html). When
  `iterate = TRUE`, a running `X/Y` counter (completed combinations out
  of the total) is printed as each combination finishes fitting. When
  `iterate = FALSE` (and `crude = FALSE`), a single message reports the
  number of observations used and whether Firth's correction was applied
  to the base model.

## Value

If `crude = TRUE`, returns a single tidy data frame with rows for all
`indep` variables stacked sequentially (separated by blank spacer rows).
Columns: Predictor, Log Odds, Std Error, p-value, Significance, OR, 95\\
Deviance, AIC, BIC, McFadden R2, CoxSnell R2, Nagelkerke R2, Tjur R2,
Firth Corrected. Inline model metrics are placed in the first
(intercept) row of each variable's block.

If `crude = FALSE` and `iterate = FALSE`, returns a data frame
containing Predictor, Log Odds, Std Error, p-value, Significance, OR,
95\\ (Deviance, AIC, BIC, McFadden R2, CoxSnell R2, Nagelkerke R2, Tjur
R2, Firth Corrected) in the first row. Model metrics are also attached
as `attr(result, "model_metrics")`. When `stepwise = TRUE`, the
variables removed at each step (and their triggering VIF) are
additionally attached as `attr(result, "stepwise_removed")`.

If `crude = FALSE` and `iterate = TRUE`, returns a named list
containing:

- Tables:

  A named list of result data frames, each with inline model metrics in
  row 1.

- Metrics:

  A data frame comparing model metrics across all combinations.

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

**Crude odds ratios (`crude = TRUE`):** Each variable in `indep` is
fitted in a separate simple logistic regression model with `y` as the
only outcome and no covariates, yielding unadjusted (crude) odds ratios.
All preprocessing steps — category exclusion (`exclude`), predictor
category removal (`remove`), dependent variable releveling (`ref`), and
predictor releveling (`ref_levels`) — are applied before each model is
fitted. The results are stacked row-wise into a single tidy data frame.
Each variable's rows are separated by a blank spacer row to visually
group the predictors. Inline model-fit metrics (Deviance, AIC, BIC, and
R2 measures) are placed in the first row of each variable's block
(aligned with its intercept row), and `NA` elsewhere. Because each model
contains only one predictor, VIF is not meaningful and the `VIF` column
will be `NA` for all rows. `confounders` and `iterate` are silently
ignored after emitting a warning.

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

**Inline model metrics:** Model-fit statistics (Deviance, AIC, BIC,
McFadden R2, CoxSnell R2, Nagelkerke R2, Tjur R2, Firth Corrected) are
appended as additional columns to the right of the `VIF` column. Values
are placed only in the first row of each model block (aligned with the
intercept); all subsequent rows contain `NA` for these columns.

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
d <- mtcars
d$am <- factor(d$am, labels = c("auto", "manual"))
d$vs <- factor(d$vs, labels = c("V", "S"))

# Adjusted model (Firth correction is applied automatically when needed)
fit <- run_logreg(d, y = "am", confounders = "wt", indep = "hp", ref = "auto")
fit
#>     Predictor   Log Odds  Std Error     p-value Significance           OR
#> 1 (Intercept) 18.8662987 7.44355806 0.011258199            * 1.561455e+08
#> 2          wt -8.0834752 3.06867511 0.008433813           ** 3.085967e-04
#> 3          hp  0.0362556 0.01773415 0.040914646            * 1.036921e+00
#>     95% CI Low  95% CI High      VIF Deviance      AIC      BIC McFadden R2
#> 1 7.203795e+01 3.384524e+14       NA 10.05911 16.05911 20.45632   0.7673104
#> 2 7.539288e-07 1.263142e-01 2.444297       NA       NA       NA          NA
#> 3 1.001498e+00 1.073596e+00 2.444297       NA       NA       NA          NA
#>   CoxSnell R2 Nagelkerke R2   Tjur R2 Firth Corrected
#> 1   0.6453351      0.870897 0.8112073           FALSE
#> 2          NA            NA        NA              NA
#> 3          NA            NA        NA              NA
attr(fit, "model_metrics")
#>        Model  N Deviance      AIC      BIC McFadden R2 CoxSnell R2
#> 1 Base Model 32 10.05911 16.05911 20.45632   0.7673104   0.6453351
#>   Nagelkerke R2   Tjur R2 Firth Corrected
#> 1      0.870897 0.8112073           FALSE

# Crude (unadjusted) odds ratios, one model per variable
run_logreg(d, y = "am", indep = c("wt", "hp", "vs"), ref = "auto",
           ref_levels = list("vs" ~ "V"), crude = TRUE)
#>     Predictor     Log Odds   Std Error     p-value Significance           OR
#> 1 (Intercept) 12.040369659 4.509706105 0.007587858           ** 1.694596e+05
#> 2          wt -4.023969940 1.436415784 0.005088198           ** 1.788183e-02
#> 3                       NA          NA          NA         <NA>           NA
#> 4 (Intercept)  0.776613552 0.915428632 0.396237197              2.174097e+00
#> 5          hp -0.008117121 0.006074337 0.181452089              9.919157e-01
#> 6                       NA          NA          NA         <NA>           NA
#> 7 (Intercept) -0.693147181 0.499999908 0.165656960              5.000000e-01
#> 8         vsS  0.693147181 0.731924992 0.343628884              2.000000e+00
#>    95% CI Low  95% CI High VIF Deviance      AIC      BIC McFadden R2
#> 1 24.56955593 1.168786e+09  NA 19.17608 23.17608 26.10756  0.55641445
#> 2  0.00107086 2.986011e-01  NA       NA       NA       NA          NA
#> 3          NA           NA  NA       NA       NA       NA          NA
#> 4  0.36146372 1.307655e+01  NA 41.22757 45.22757 48.15904  0.04631460
#> 5  0.98017652 1.003796e+00  NA       NA       NA       NA          NA
#> 6          NA           NA  NA       NA       NA       NA          NA
#> 7  0.18765896 1.332204e+00  NA 42.32263 46.32263 49.25410  0.02098329
#> 8  0.47644669 8.395483e+00  NA       NA       NA       NA          NA
#>   CoxSnell R2 Nagelkerke R2    Tjur R2 Firth Corrected
#> 1  0.52842471    0.71312328 0.62400276           FALSE
#> 2          NA            NA         NA              NA
#> 3          NA            NA         NA              NA
#> 4  0.06065057    0.08184957 0.07202501           FALSE
#> 5          NA            NA         NA              NA
#> 6          NA            NA         NA              NA
#> 7  0.02794894    0.03771784 0.02834008           FALSE
#> 8          NA            NA         NA              NA
```
