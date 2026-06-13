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
  positive_only = TRUE,
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

  - **Numeric** — fraction of the smallest observed value per feature
    used for deterministic imputation (e.g., `0.2` = 1/5th of the
    smallest value). By default (`positive_only = TRUE`), the minimum is
    computed from strictly positive values only, which is the standard
    MNAR convention in metabolomics. Set `positive_only = FALSE` to use
    the global minimum including zeros and negative values. Recommended
    range: 0.01–0.5. Values \> 1 trigger a warning.

  - **`"quantileregression"`** — uses the QRILC algorithm from the
    imputeLCMD package. Assumes Missing Not At Random (MNAR) data.
    Controlled further by `tune_sigma`.

  - **Any `mice` method string** — e.g., `"pmm"`, `"norm"`,
    `"norm.boot"`, `"norm.nob"`, `"mean"`, `"2l.norm"`, `"rf"`,
    `"cart"`, etc. Delegates to
    [`mice`](https://amices.org/mice/reference/mice.html) from the mice
    package. Controlled further by `m`, `maxit`, `seed`, and `...`.

  - **`"none"`** — Skips missing value imputation entirely and returns
    the data exactly as-is. Useful for deferring imputation to a later
    step.

  Default: `0.2` (1/5th of minimum, deterministic).

- positive_only:

  Logical. Only used when `method` is numeric. If `TRUE` (default), the
  per-feature minimum is computed from strictly positive values only
  (zeros and negatives are excluded), which is standard practice for
  MNAR imputation in metabolomics where missing values represent
  abundances below the detection limit. If `FALSE`, the global minimum
  across all observed (non-NA) values is used, which may be appropriate
  when features can legitimately be zero or negative (e.g.,
  log-transformed or mean-centred data). Has no effect when `method` is
  `"quantileregression"` or a mice method string.

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

  List of parameters used, including method, positive_only, tune_sigma
  (QRILC only), and m/maxit/seed (mice only).

- method_used:

  Character string describing the imputation method applied.

## Details

**Missing Value Imputation Strategies:**

1.  **Deterministic (fraction of minimum):** When `method` is numeric,
    each missing value in a feature is replaced by a fraction of the
    per-feature minimum, multiplied by the specified fraction. The
    minimum is determined by `positive_only`:

    - `positive_only = TRUE` (default): uses the smallest strictly
      positive observed value per feature, excluding zeros and
      negatives. This is standard in metabolomics for values below the
      detection limit (MNAR — Missing Not At Random).

    - `positive_only = FALSE`: uses the global minimum across all
      observed (non-NA) values, including zeros and negatives. Suitable
      for data that has been log-transformed or mean-centred.

    If no valid minimum can be found (e.g., a feature is entirely NA or,
    when `positive_only = TRUE`, has no positive values), a fallback of
    `1e-9` is used.

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

[`run_DIpreprocess`](https://jllcalorio.github.io/pondeR/reference/run_DIpreprocess.md)

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

# 1. Deterministic imputation (1/5th of minimum positive value per feature)
result1 <- run_mvimpute(x, method = 0.2)

# 2. Deterministic imputation using global minimum (including zeros/negatives)
result2 <- run_mvimpute(x, method = 0.2, positive_only = FALSE)

# 3. Quantile regression imputation (QRILC)
result3 <- run_mvimpute(x, method = "quantileregression", tune_sigma = 1)

# 4. Predictive mean matching via mice
result4 <- run_mvimpute(x, method = "pmm", m = 5, maxit = 5, seed = 42)

# 5. Random forest imputation via mice
result5 <- run_mvimpute(x, method = "rf", m = 5, maxit = 5, seed = 42)

# 6. Passing additional mice arguments via ...
result6 <- run_mvimpute(x, method = "norm.boot", m = 10, printFlag = FALSE)
} # }
```
