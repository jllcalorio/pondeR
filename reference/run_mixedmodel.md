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

## References

Bates D, Machler M, Bolker B, Walker S (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
67(1), 1-48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)

Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package:
tests in linear mixed effects models. *Journal of Statistical Software*,
82(13), 1-26.
[doi:10.18637/jss.v082.i13](https://doi.org/10.18637/jss.v082.i13)

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *Journal of the
Royal Statistical Society: Series B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## Author

John Lennon L. Calorio

## Examples

``` r
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
#> === run_mixedmodel ===
#>   Features analysed : 3
#>   Successful models : 3
#>   Failed / skipped  : 0
#>   Success rate      : 100.0%
summary(result)
#> ---------------------------------------
#> run_mixedmodel - Analysis Summary
#> ---------------------------------------
#>   Formula       : .feature ~ time + age + (1 | subject_id)
#>   Subject col   : subject_id
#>   P-adj method  : BH
#>   Features in   : 3
#>   Converged     : 3
#>   Failed        : 0
#> 
#> -- Top Features by Adjusted P-value --
#>    Feature        Term      Estimate      p-value Adjusted p-value Significance
#>    glucose (Intercept)  5.417529e+00 4.463334e-13     1.339000e-12          ***
#>    insulin (Intercept)  8.477767e+01 1.029527e-09     1.544290e-09          ***
#>  PA O-28:1 (Intercept)  1.314088e+00 1.724036e-04     1.724036e-04          ***
#>    insulin     timepre -1.249302e+01 1.701740e-04     5.105219e-04          ***
#>    glucose     timepre -5.272202e-01 9.172948e-03     1.375942e-02            *
#>  PA O-28:1     timepre -1.146388e-01 1.904673e-01     1.904673e-01             
#>    insulin         age  1.708731e-01 2.997114e-01     8.991343e-01             
#>    glucose         age  1.155782e-04 9.911166e-01     9.911166e-01             
#>  PA O-28:1         age -2.594189e-03 6.639419e-01     9.911166e-01             
head(result$combined_fixed_effects)
#>   Feature        Term      Estimate Std. error       df           t
#> 1 glucose (Intercept)  5.417529e+00 0.49825353 37.00000 10.87303635
#> 2 glucose         age  1.155782e-04 0.01031082 37.00000  0.01120942
#> 3 glucose     timepre -5.272202e-01 0.19173681 37.00000 -2.74970806
#> 4 insulin (Intercept)  8.477767e+01 7.70562141 19.12266 11.00205465
#> 5 insulin         age  1.708731e-01 0.16001951 18.00000  1.06782656
#> 6 insulin     timepre -1.249302e+01 2.68017781 19.00000 -4.66126710
#>        p-value 95% CI lower 95% CI upper Adjusted p-value Significance
#> 1 4.463334e-13   4.44096979    6.3940878     1.339000e-12          ***
#> 2 9.911166e-01  -0.02009325    0.0203244     9.911166e-01             
#> 3 9.172948e-03  -0.90301748   -0.1514230     1.375942e-02            *
#> 4 1.029527e-09  69.67492741   99.8804083     1.544290e-09          ***
#> 5 2.997114e-01  -0.14275940    0.4845056     8.991343e-01             
#> 6 1.701740e-04 -17.74607665   -7.2399727     5.105219e-04          ***
```
