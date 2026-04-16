# Impute Missing Values in Metabolomics Data

Replaces missing values (NAs) using either a simple deterministic method
(fraction of minimum), a statistical imputation method (quantile
regression imputation of left-censored data, QRILC), or any method
available in [`mice`](https://amices.org/mice/reference/mice.html)
(Multiple Imputation by Chained Equations).

## Usage

``` r
run_mvimpute(
  x,
  method = 0.2,
  tune_sigma = 1,
  m = 5,
  maxit = 5,
  seed = NA,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns. Missing values should be represented as `NA`.

- method:

  Numeric or character. Controls the imputation strategy:

  - **Numeric** — fraction of the smallest positive value per feature
    used for deterministic imputation (e.g., `0.2` = 1/5th of the
    smallest positive value). Recommended range: 0.01–0.5. Values \> 1
    trigger a warning.

  - **`"quantileregression"`** — uses the QRILC algorithm from the
    imputeLCMD package. Assumes Missing Not At Random (MNAR) data.
    Controlled further by `tune_sigma`.

  - **Any `mice` method string** — e.g., `"pmm"`, `"norm"`,
    `"norm.boot"`, `"norm.nob"`, `"mean"`, `"2l.norm"`, `"rf"`,
    `"cart"`, etc. Delegates to
    [`mice`](https://amices.org/mice/reference/mice.html) from the mice
    package. Controlled further by `m`, `maxit`, `seed`, and `...`.

  Default: `0.2` (1/5th of minimum, deterministic).

- tune_sigma:

  Numeric. Only used when `method = "quantileregression"`. Controls the
  standard deviation of the left-censored missing value distribution.
  Smaller values (e.g., `0.5`) assume a tighter distribution; larger
  values (e.g., `2`) allow more spread. Default: `1`.

- m:

  Integer. Only used when `method` routes to mice. Number of multiple
  imputations. The first completed dataset is used as the imputed
  result. Default: `5`.

- maxit:

  Integer. Only used when `method` routes to mice. Number of iterations
  for the chained equations algorithm. Default: `5`.

- seed:

  Integer or `NA`. Only used when `method` routes to mice. Random seed
  passed to [`mice`](https://amices.org/mice/reference/mice.html) for
  reproducibility. Default: `NA` (no seed set).

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

- ...:

  Additional arguments passed to
  [`mice`](https://amices.org/mice/reference/mice.html) when `method`
  routes to the mice engine. Has no effect for numeric or
  `"quantileregression"` methods. A notable argument is
  `predictorMatrix`: by default, `run_mvimpute` builds one automatically
  via [`quickpred`](https://amices.org/mice/reference/quickpred.html) to
  avoid failures caused by constant or collinear features, which are
  common in wide metabolomics/metagenomics data. Supply your own
  `predictorMatrix` via `...` to override this.

## Value

A list with the following elements:

- data:

  Matrix or data frame with imputed values (same class as input `x`).

- n_missing_before:

  Integer. Number of missing values before imputation.

- n_missing_after:

  Integer. Number of missing values after imputation (should be 0).

- imputed_summary:

  Data frame with per-feature imputation statistics.

- parameters:

  List of parameters used, including method, tune_sigma (QRILC only),
  and m/maxit/seed (mice only).

- method_used:

  Character string describing the imputation method applied.

## Details

**Missing Value Imputation Strategies:**

1.  **Deterministic (fraction of minimum):** When `method` is numeric,
    each missing value in a feature is replaced by the smallest strictly
    positive value observed in that feature, multiplied by the specified
    fraction. Zero and negative values are excluded when determining
    this minimum. This assumes missing values represent low abundance
    (MNAR — Missing Not At Random) and is standard practice in
    metabolomics for values below the detection limit.

2.  **Quantile Regression (QRILC):** When
    `method = "quantileregression"`, uses the Quantile Regression
    Imputation of Left-Censored data algorithm from the imputeLCMD
    package. Models missing values from a left-censored distribution,
    controlled by `tune_sigma`. More statistically principled than
    simple fraction-of-minimum for MNAR data. Requires imputeLCMD:
    `BiocManager::install("imputeLCMD")`.

3.  **mice (Multiple Imputation by Chained Equations):** When `method`
    is a recognised mice method string,
    [`mice`](https://amices.org/mice/reference/mice.html) is called
    internally on the feature matrix. The first of the `m` completed
    datasets is extracted via
    [`mice::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
    and returned. Supported `mice` methods include (non-exhaustive):
    `"pmm"` (predictive mean matching), `"norm"`, `"norm.boot"`,
    `"norm.nob"`, `"mean"`, `"rf"` (random forest), `"cart"`,
    `"2l.norm"`, `"midastouch"`, `"sample"`. See
    [`mice`](https://amices.org/mice/reference/mice.html) for the full
    list. Requires mice: `install.packages("mice")`. By default,
    `run_mvimpute` builds a predictor matrix automatically via
    [`quickpred`](https://amices.org/mice/reference/quickpred.html) to
    guard against the common failure where `mice` removes all predictors
    due to constant or collinear columns (frequent in wide, sparse
    metabolomics/metagenomics data). Supply a custom `predictorMatrix`
    via `...` to override. Additional arguments for mice can also be
    passed via `...`.

## References

Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
Imputation by Chained Equations in R. *Journal of Statistical Software*,
45(3), 1–67.
[doi:10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03)

Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
Missing Value Imputation Approach for Mass Spectrometry-based
Metabolomics Data. *Scientific Reports*, 8(1), 663.
[doi:10.1038/s41598-017-19120-0](https://doi.org/10.1038/s41598-017-19120-0)

Lazar, C., Gatto, L., Ferro, M., Bruley, C., & Burger, T. (2016).
Accounting for the Multiple Natures of Missing Values in Label-Free
Quantitative Proteomics Data Sets to Compare Imputation Strategies.
*Journal of Proteome Research*, 15(4), 1116–1125.
[doi:10.1021/acs.jproteome.5b00981](https://doi.org/10.1021/acs.jproteome.5b00981)

## See also

[`mice`](https://amices.org/mice/reference/mice.html),
[`impute.QRILC`](https://rdrr.io/pkg/imputeLCMD/man/impute.QRILC.html)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with missing values
set.seed(123)
x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
x[sample(length(x), 200)] <- NA
colnames(x) <- paste0("Feature", 1:50)

# 1. Deterministic imputation (1/5th of minimum per feature)
result1 <- run_mvimpute(x, method = 0.2)

# 2. Quantile regression imputation (QRILC)
result2 <- run_mvimpute(x, method = "quantileregression", tune_sigma = 1)

# 3. Predictive mean matching via mice
result3 <- run_mvimpute(x, method = "pmm", m = 5, maxit = 5, seed = 42)

# 4. Random forest imputation via mice
result4 <- run_mvimpute(x, method = "rf", m = 5, maxit = 5, seed = 42)

# 5. Passing additional mice arguments via ...
result5 <- run_mvimpute(x, method = "norm.boot", m = 10, printFlag = FALSE)
} # }
```
