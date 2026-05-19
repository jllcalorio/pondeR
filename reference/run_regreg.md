# Regularized Regression with Cross-Validation (Ridge, Elastic Net, LASSO)

Fits regularized regression models — Ridge (alpha = 0), LASSO (alpha =
1), or Elastic Net (0 \< alpha \< 1) — with k-fold cross-validation via
[`cv.glmnet`](https://rdrr.io/pkg/glmnet/man/cv.glmnet.html). Supports
both categorical (binary and multinomial classification) and continuous
(Gaussian regression) dependent variables. When multiple alpha values
are supplied, all are evaluated on identical cross-validation folds and
the best-performing model is returned alongside a comparison table.

Unpenalized covariates (e.g., batch indicators, demographic confounders)
can be included through the `not_penalized` argument; their coefficients
are estimated freely while all other predictors remain subject to
regularization.

## Usage

``` r
run_regreg(
  x,
  metadata,
  pred,
  ref = NULL,
  not_penalized = NULL,
  is_numeric = NULL,
  train_percent = 0.8,
  alpha = 0.5,
  lambda = "1se",
  cv_folds = 10L,
  type_measure = NULL,
  parallel = FALSE,
  standardize = FALSE,
  maxit = 100000L,
  seed = 123
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, or `matrix` of predictor variables. Must
  have column names (special characters are supported). All columns are
  treated as features to be penalized unless overridden by
  `not_penalized`. Rows must correspond to the same samples as
  `metadata`.

- metadata:

  A `data.frame` or `tibble` containing sample-level metadata. Must have
  the same number of rows as `x` and in the same order. Must contain the
  column specified by `pred`.

- pred:

  A single character string naming a column in `metadata` to use as the
  dependent variable. If the column is numeric, Gaussian regression is
  performed. If it is character or factor, classification is performed.

- ref:

  A single character string specifying the reference category for
  categorical `pred`. Required for classification; ignored for numeric
  `pred`.

- not_penalized:

  Optional character vector of column names from `metadata` and/or `x`
  to include in the model without penalization (penalty factor = 0).
  Useful for known confounders or covariates that must be controlled.
  The column named in `pred` is automatically excluded with a warning if
  accidentally included (e.g. when passing `names(metadata)`). Default:
  `NULL`.

- is_numeric:

  Optional character vector of column names that are listed in
  `not_penalized` and should be treated as already numeric, bypassing
  the automatic non-numeric detection and integer-code conversion. Use
  this when a column is genuinely numeric but is stored as `character`
  or `factor` (e.g. due to upstream data import). All names must exist
  in `metadata` or `x` and must also appear in `not_penalized`. Default:
  `NULL`.

- train_percent:

  A numeric scalar in `(0, 1)` specifying the proportion of samples used
  for model training. The remaining `1 - train_percent` is held out for
  evaluation. Default: `0.8`.

- alpha:

  A numeric scalar or vector with values in `[0, 1]`:

  0

  :   Ridge regression (L2 penalty; retains all predictors)

  1

  :   LASSO regression (L1 penalty; performs variable selection)

  (0, 1)

  :   Elastic Net (blend of L1 and L2)

  When a vector is supplied all values are evaluated and the best model
  is selected. Default: `0.5`.

- lambda:

  A single string controlling lambda selection after cross-validation:

  `"1se"`

  :   Lambda within one standard error of the minimum CV error (more
      conservative; fewer selected predictors). Default.

  `"min"`

  :   Lambda minimizing CV error (less conservative; more selected
      predictors).

- cv_folds:

  A positive integer between 3 and 20 specifying the number of
  cross-validation folds. Default: `10`.

- type_measure:

  A character string specifying the CV loss function. `NULL` (default)
  selects automatically based on `pred` type:

  Categorical

  :   `"class"` (misclassification error). Also accepts `"auc"` (binary
      only; falls back to `"class"` for multinomial) or `"deviance"`.

  Numeric

  :   `"mse"` (mean squared error). Also accepts `"mae"` or
      `"deviance"`.

- parallel:

  Logical. If `TRUE`, cross-validation is parallelized using a
  registered foreach backend (e.g., via
  [`doParallel::registerDoParallel()`](https://rdrr.io/pkg/doParallel/man/registerDoParallel.html)).
  A warning is issued if `TRUE` but no backend is registered. Default:
  `FALSE`.

- standardize:

  Logical. Whether to standardize predictors inside `glmnet` prior to
  fitting. Set to `FALSE` (default) when data have already been
  preprocessed/scaled. Default: `FALSE`.

- maxit:

  A positive integer specifying the maximum number of iterations for
  coordinate descent convergence. Default: `100000`.

- seed:

  A single numeric value passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducibility
  of the train/test split and cross-validation fold assignment. Set to
  `NULL` to disable. Default: `123`.

## Value

A named list of class `"run_regreg"` with the following elements:

- `modelsummary`:

  A `data.frame` summarising performance metrics (and a rank column) for
  every alpha tested.

- `sampledistribution`:

  A `data.frame` showing the number of training and testing samples per
  group (classification) or overall (regression), with a `Total` row
  appended.

- `BestModel`:

  A named list for the best-performing alpha:

  `Model`

  :   The fitted `cv.glmnet` object.

  `Performance`

  :   A one-row `data.frame` of metrics (Accuracy + Kappa for
      classification; MSE, RMSE, MAE, R2 for regression).

  `Coefficients`

  :   A `data.frame` of non-zero coefficients (intercept included if
      non-zero). Columns: `Feature`, `Coefficient`, and `OddsRatio`
      (classification) or just `Feature` and `Coefficient` (regression).
      Multinomial models include a leading `Comparison` column.

  `N_Coefficients`

  :   Integer; total non-zero coefficients including intercept(s).

  `ConfusionMatrix`

  :   A `confusionMatrix` object (caret); `NULL` for regression.

  `Predictions`

  :   A `data.frame` with columns `Actual` and `Predicted`.

  `Lambda`

  :   Numeric; the selected lambda value.

  `Alpha`

  :   Numeric; the alpha value of this model.

- `AllModels`:

  A named list (`"alpha_<value>"`) with the same structure as
  `BestModel` for every alpha tested.

- `AlphaComparison`:

  A `data.frame` comparing all tested alpha values by performance metric
  and number of non-zero coefficients. `NULL` when only one alpha is
  tested.

## Details

Perform Regularized Regression Analysis

**Data assembly:** The feature matrix is taken from `x`. If
`not_penalized` column names belonging to `metadata` are provided, those
columns are appended to the feature matrix. A `penalty.factor` vector of
length `ncol(x)` is set to `1` for all penalized predictors and `0` for
unpenalized ones, and is passed to
[`glmnet::cv.glmnet`](https://rdrr.io/pkg/glmnet/man/cv.glmnet.html).

**Cross-validation consistency:** Fold IDs are generated once (using
`seed`) and reused across all tested alpha values, ensuring a fair,
apple-to-apple comparison of regularization paths.

**Multinomial coefficients:** For multinomial models the raw coefficient
matrix is normalized relative to the reference level so that each
comparison reflects log-odds vs. the reference. Custom pairwise
contrasts are not supported in this version; use the `Coefficients` data
frame and compute differences manually.

**Why no p-values?** P-values are not provided because they are
statistically invalid after regularized variable selection. Report
cross-validated performance metrics and coefficient magnitudes instead.

## References

Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths
for Generalized Linear Models via Coordinate Descent. *Journal of
Statistical Software*, 33(1), 1–22.
[doi:10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01) .

Zou, H. and Hastie, T. (2005) Regularization and variable selection via
the elastic net. *Journal of the Royal Statistical Society: Series B*,
67(2), 301–320.

Tibshirani, R. (1996) Regression shrinkage and selection via the lasso.
*Journal of the Royal Statistical Society: Series B*, 58(1), 267–288.

Kuhn, M. (2008) Building predictive models in R using the caret package.
*Journal of Statistical Software*.
[doi:10.18637/jss.v028.i05](https://doi.org/10.18637/jss.v028.i05) .

Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011)
Regularization Paths for Cox's Proportional Hazards Model via Coordinate
Descent. *Journal of Statistical Software*, 39(5), 1–13.
[doi:10.18637/jss.v039.i05](https://doi.org/10.18637/jss.v039.i05) .

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
## ---- Simulate data -------------------------------------------------------
set.seed(42)
n  <- 80
p  <- 30

feat <- as.data.frame(
  matrix(rnorm(n * p), nrow = n,
         dimnames = list(NULL, paste0("feat_", seq_len(p))))
)

meta <- data.frame(
  SampleID = paste0("S", seq_len(n)),
  Group    = sample(c("Control", "Case"), n, replace = TRUE),
  Age      = round(runif(n, 20, 65))
)

## ---- Binary classification (Elastic Net, single alpha) -------------------
res_bin <- run_regreg(
  x             = feat,
  metadata      = meta,
  pred          = "Group",
  ref           = "Control",
  alpha         = 0.5,
  train_percent = 0.8,
  seed          = 123
)
print(res_bin)
res_bin$modelsummary
res_bin$BestModel$Coefficients

## ---- Binary classification (multiple alpha sweep) ------------------------
res_sweep <- run_regreg(
  x             = feat,
  metadata      = meta,
  pred          = "Group",
  ref           = "Control",
  alpha         = c(0, 0.25, 0.5, 0.75, 1),
  lambda        = "1se",
  seed          = 123
)
res_sweep$AlphaComparison

## ---- LASSO with an unpenalized covariate ---------------------------------
res_lasso <- run_regreg(
  x             = feat,
  metadata      = meta,
  pred          = "Group",
  ref           = "Control",
  not_penalized = "Age",
  alpha         = 1,
  seed          = 123
)

## ---- Gaussian regression (numeric pred) ----------------------------------
res_num <- run_regreg(
  x             = feat,
  metadata      = meta,
  pred          = "Age",
  alpha         = c(0, 0.5, 1),
  train_percent = 0.75,
  seed          = 7
)
res_num$BestModel$Performance

## ---- Multinomial classification ------------------------------------------
meta3 <- meta
meta3$Group <- sample(c("A", "B", "C"), n, replace = TRUE)

res_multi <- run_regreg(
  x             = feat,
  metadata      = meta3,
  pred          = "Group",
  ref           = "A",
  alpha         = 0.5,
  seed          = 1
)
res_multi$BestModel$Coefficients
} # }
```
