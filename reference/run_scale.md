# Scale Metabolomics Data

Applies various scaling methods to standardize features and prepare data
for multivariate analysis. Includes classical scaling methods as well as
variance-stability (VAST) scaling variants.

## Usage

``` r
run_scale(x, method = "auto", metadata = NULL, group = "Group", verbose = TRUE)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- method:

  Character. Scaling method to apply. Options:

  - `"mean"`: Mean-centering only (centers distribution at zero)

  - `"auto"`: Auto-scaling (mean-centering + unit variance scaling)

  - `"pareto"`: Pareto-scaling (mean-centering + sqrt(SD) scaling)

  - `"ln"`: Logarithm (ln) transformation

  - `"vast"`: VAST scaling (auto-scaling weighted by CV)

  - `"svast"`: Supervised VAST scaling (auto-scaling weighted by mean
    class CV)

  - `"xvast"`: Extended supervised VAST scaling (auto-scaling weighted
    by maximum class CV)

  - `"none"`: No scaling method applied

- metadata:

  Data frame. Sample metadata with `nrow(metadata) == nrow(x)`. Required
  when `method` is `"svast"` or `"xvast"`; ignored otherwise. Default:
  `NULL`.

- group:

  Character. Name of the column in `metadata` that contains group/class
  labels. Required when `method` is `"svast"` or `"xvast"`. Default:
  `"Group"`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

A list of class `"run_scale"` containing:

- data:

  Matrix or data frame of scaled data (same class as input `x`)

- scaling_factors:

  Numeric vector of scaling factors applied to each feature

- center_values:

  Numeric vector of centering values (means) for each feature

- method_used:

  Character describing scaling method applied

- parameters:

  List of parameters used

## Details

**Scaling Methods:**

- **mean**: Mean-centering only. Each feature is centered to have mean =
  0 but retains original variance. Useful when features are already on
  similar scales and you want to preserve relative variance information.

- **auto**: Auto-scaling (also called standardization or unit variance
  scaling). Each feature is mean-centered and divided by its standard
  deviation, resulting in mean = 0 and SD = 1 for all features. This
  gives equal weight to all features regardless of their original
  variance.

  **Recommended for:** PCA, hierarchical clustering, correlation
  analysis, t-tests and ANOVA, and most univariate and multivariate
  methods. Note that auto-scaling can amplify noise in low-variance
  features.

- **pareto**: Pareto-scaling. Each feature is mean-centered and divided
  by the square root of its standard deviation. This provides a
  compromise between no scaling and auto-scaling, reducing the influence
  of large values while preserving some of the original variance
  structure.

  **Recommended for:** PLS-DA, OPLS-DA, sPLS-DA, and other PLS-type
  methods.

- **ln**: Logarithm (ln) transformation. Each feature is transformed
  using the natural logarithm. This is a non-linear transformation
  included for cases where log-transformation is desired as part of the
  scaling step (e.g., to reduce the influence of large values).

- **vast**: VAST (Variable Stability) scaling. Extends auto-scaling by
  additionally weighting each feature by its overall coefficient of
  variation (CV = mean / SD). The scaling factor is \\s_k / ({\bar{x}\_k
  / s_k})\\ = \\s_k^2 / \bar{x}\_k\\. Features with stable measurements
  (low CV) receive higher weights than noisy features. Equivalent to
  dividing auto-scaled data by the overall CV.

  **Recommended for:** NMR-based metabolic profiling; data where
  measurement stability varies across features.

- **svast**: Supervised VAST (s-VAST) scaling. A group-aware extension
  of VAST that weights each feature by the average of its within-class
  CVs (\\(1/n)\sum\_{j=1}^{n} \bar{x}\_{jk}/s\_{jk}\\), rewarding
  features that are stable within classes. Requires `group`.

- **xvast**: Extended supervised VAST (x-VAST) scaling, introduced by
  Yang et al. (2015). Weights each feature by the *maximum* within-class
  ratio \\\max(\bar{x}\_{1k}/s\_{1k}, \ldots, \bar{x}\_{nk}/s\_{nk})\\,
  placing emphasis on the class in which a feature is most stable.
  Requires `group`.

**Choosing a Scaling Method:**

The choice depends on your downstream analysis:

- For exploratory data analysis (PCA, clustering): **auto-scaling**

- For classification/discrimination (PLS-DA): **Pareto-scaling**

- When feature magnitudes carry meaning: **mean-centering**

- For NMR/metabolomics with variable measurement stability: **VAST**

- When class-specific stability is important: **s-VAST** or **x-VAST**

## References

van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A., Smilde, A.K., &
van der Werf, M.J. (2006). Centering, scaling, and transformations:
improving the biological information content of metabolomics data. *BMC
Genomics*, 7, 142.
[doi:10.1186/1471-2164-7-142](https://doi.org/10.1186/1471-2164-7-142)

Keun, H.C., Ebbels, T.M.D., Antti, H., Bollard, M.E., Beckonert, O.,
Holmes, E., et al. (2003). Improved analysis of multivariate data by
variable stability scaling: application to NMR-based metabolic
profiling. *Analytica Chimica Acta*, 490, 265–276.
[doi:10.1016/S0003-2670(03)00094-1](https://doi.org/10.1016/S0003-2670%2803%2900094-1)

Yang, J., Zhao, X., Lu, X., Lin, X., & Xu, G. (2015). A data
preprocessing strategy for metabolomics to reduce the mask effect in
data analysis. *Frontiers in Molecular Biosciences*, 2, 4.
[doi:10.3389/fmolb.2015.00004](https://doi.org/10.3389/fmolb.2015.00004)

Becker, R.A., Chambers, J.M., & Wilks, A.R. (1988). *The New S
Language*. Wadsworth & Brooks/Cole.
[doi:10.1201/9781351074988](https://doi.org/10.1201/9781351074988)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/pondeR/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(rnorm(100 * 50, mean = 100, sd = 50),
            nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)
group <- rep(c("A", "B"), each = 50)

# Auto-scaling for PCA
result1 <- run_scale(x, method = "auto")

# Pareto-scaling for PLS-DA
result2 <- run_scale(x, method = "pareto")

# Log (ln) transformation
result_ln <- run_scale(x, method = "ln")

# Mean-centering only
result3 <- run_scale(x, method = "mean")

# VAST scaling
result4 <- run_scale(x, method = "vast")

# Supervised VAST scaling
result5 <- run_scale(x, method = "svast", group = group)

# Extended supervised VAST scaling (Yang et al., 2015)
result6 <- run_scale(x, method = "xvast", group = group)
} # }
```
