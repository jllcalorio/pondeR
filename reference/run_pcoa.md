# Principal Coordinate Analysis (PCoA) with Optional PERMANOVA and Post-hoc

Performs Principal Coordinate Analysis (PCoA), also known as metric
multidimensional scaling (MDS), as a sequential wrapper around
[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
and [`pcoa`](https://rdrr.io/pkg/ape/man/pcoa.html). Optionally runs a
Permutational Multivariate Analysis of Variance (PERMANOVA) via
[`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html) on
the same dissimilarity matrix. Can also perform pairwise post-hoc
testing using
[`OmicFlow::pairwise_adonis`](https://rdrr.io/pkg/OmicFlow/man/pairwise_adonis.html).
Both analyses share a single dissimilarity matrix so the method is
applied exactly once.

## Usage

``` r
run_pcoa(
  x,
  method = "bray",
  correction = "none",
  metadata = NULL,
  rhs = NULL,
  perms = 9999L,
  groups = NULL,
  p_adjust = "BH",
  ...
)
```

## Arguments

- x:

  A numeric `data.frame`, `tibble`, or `matrix` of observations (rows)
  by variables (columns). Column names may contain special characters.
  Alternatively, a precomputed dissimilarity object of class `"dist"`;
  when `x` is a `"dist"` object the `method` argument is ignored and the
  matrix is used directly for both PCoA and PERMANOVA.

- method:

  A single character string specifying the dissimilarity index passed to
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html).
  Common choices include `"bray"` (Bray–Curtis), `"euclidean"`,
  `"aitchison"`, and `"robust.aitchison"`; see
  [`?vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  for the full list. Ignored when `x` is already a `"dist"` object. The
  same dissimilarity matrix is reused for PERMANOVA — the method is
  never applied twice. Default: `"bray"`.

- correction:

  Correction method for negative eigenvalues passed to
  [`pcoa`](https://rdrr.io/pkg/ape/man/pcoa.html). One of `"lingoes"`,
  `"cailliez"`, or `"none"`. Default: `"none"`.

- metadata:

  A `data.frame` with one row per observation, in the same row order as
  `x`. Required when `rhs` is specified (PERMANOVA); ignored otherwise.
  The grouping variable(s) named in `rhs` must be columns of `metadata`.

- rhs:

  A single character string naming a column of `metadata` to use as the
  right-hand side of the PERMANOVA formula (`dist_matrix ~ rhs`). When
  `NULL` (default), PERMANOVA is skipped entirely.

- perms:

  A positive integer. Number of permutations for PERMANOVA. Passed as
  the `permutations` argument to
  [`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html).
  Default: `9999`.

- groups:

  A single character string naming a column of `metadata` to use for
  pairwise post-hoc testing via
  [`OmicFlow::pairwise_adonis`](https://rdrr.io/pkg/OmicFlow/man/pairwise_adonis.html).
  When `NULL` (default), post-hoc analysis is skipped.

- p_adjust:

  A single character string specifying the p-value adjustment method for
  post-hoc testing. Passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default: `"BH"`
  (Benjamini-Hochberg).

- ...:

  Additional arguments forwarded to
  [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  (e.g., `binary`, `na.rm`) when `x` is not a `"dist"` object, to
  [`pcoa`](https://rdrr.io/pkg/ape/man/pcoa.html) (e.g., `rn`), and to
  [`adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html)
  (e.g., `by`, `parallel`, `strata`) when PERMANOVA is requested, and to
  [`pairwise_adonis`](https://rdrr.io/pkg/OmicFlow/man/pairwise_adonis.html).
  Arguments are routed to each function by matching against their
  respective [`formals()`](https://rdrr.io/r/base/formals.html);
  unrecognised names trigger a warning.

## Value

A named list of class `c("run_pcoa", "list")` with the following
elements:

- `scores`:

  A `data.frame` of principal coordinate scores, columns `PC1`, `PC2`,
  ... Row names match those of `x`.

- `variance_explained`:

  Named numeric vector of percentage variance explained per axis (`PC1`,
  `PC2`, ...).

- `eigenvalues`:

  Raw eigenvalue vector from
  [`ape::pcoa`](https://rdrr.io/pkg/ape/man/pcoa.html).

- `dist_matrix`:

  The `"dist"` object used for both PCoA and PERMANOVA.

- `method`:

  Dissimilarity method used, or `"user-supplied dist"`.

- `correction`:

  Correction method applied.

- `n_obs`:

  Integer. Number of observations.

- `n_axes`:

  Integer. Number of positive principal coordinate axes.

- `pcoa_raw`:

  Complete [`ape::pcoa`](https://rdrr.io/pkg/ape/man/pcoa.html) output
  for downstream compatibility.

- `vegdist_raw`:

  `"dist"` object from
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html),
  or `NULL` when `x` was already a `"dist"` object.

- `permanova`:

  When `rhs` is not `NULL`: a list with elements `table` (the full
  `adonis2` result as a `data.frame`), `R2` (numeric, proportion of
  variance for the grouping term), `F_stat` (numeric, pseudo-F),
  `p_value` (numeric, permutation p-value), `rhs` (character, the
  grouping variable name), and `perms` (integer, number of permutations
  used). `NULL` when PERMANOVA was not requested.

- `pairwise_adonis`:

  When `groups` is specified: the full result table from
  [`OmicFlow::pairwise_adonis`](https://rdrr.io/pkg/OmicFlow/man/pairwise_adonis.html).

- `call`:

  The matched function call.

## Details

### Algorithm

When `x` is a numeric matrix or data frame the function proceeds as:

1.  **Dissimilarity**: `vegan::vegdist(x, method = method, ...)`
    produces a `"dist"` object.

2.  **PCoA**: `ape::pcoa(dist_matrix, correction = correction, ...)`
    decomposes the matrix into principal coordinates.

3.  **PERMANOVA** (when `rhs` is not `NULL`):
    `vegan::adonis2(dist_matrix ~ rhs, data = metadata, ...)` tests
    group differences using the already-computed dissimilarity matrix.
    Because the `"dist"` object is passed directly, `method` is not
    re-applied.

When `x` is already a `"dist"` object, stage 1 is skipped.

### PERMANOVA

PERMANOVA (Anderson, 2001; McArdle & Anderson, 2001) partitions
multivariate variance using a distance matrix and tests significance by
permuting rows. It is robust to non-normality but sensitive to
heterogeneous within-group dispersions (Warton et al., 2012); use
alongside a test of homogeneity of dispersions
([`vegan::betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html))
when groups differ in spread.

Key output fields from `adonis2`:

- **R²**: proportion of total variance explained by the grouping factor.

- **F**: pseudo-F statistic.

- **Pr(\>F)**: permutation p-value.

### Negative eigenvalues

Non-Euclidean dissimilarity coefficients may produce negative
eigenvalues. Apply `correction = "lingoes"` or `"cailliez"` to obtain a
fully Euclidean representation (Lingoes, 1971; Cailliez, 1983).

### Variance explained

Computed as \\\lambda_i / \sum \lambda^+\\, using only positive
eigenvalues, consistent with ape convention.

### Compatibility with `plot_score`

The returned object carries class `"run_pcoa"` and a `scores` element
with columns `PC1`, `PC2`, ... so that
[`pondeR::plot_score`](https://jllcalorio.github.io/pondeR/reference/plot_score.md)
can dispatch correctly. When PERMANOVA results are present,
`plot_score.run_pcoa` formats and injects them automatically as the plot
subtitle (R² and p-value).

## References

Anderson, M.J. (2001) A new method for non-parametric multivariate
analysis of variance. *Austral Ecology*, **26**, 32–46.

Cailliez, F. (1983) The analytical solution of the additive constant
problem. *Psychometrika*, **48**, 305–308.

Excoffier, L., Smouse, P.E., and Quattro, J.M. (1992) Analysis of
molecular variance inferred from metric distances among DNA haplotypes:
Application to human mitochondrial DNA restriction data. *Genetics*,
**131**, 479–491.

Gower, J.C. (1966) Some distance properties of latent root and vector
methods used in multivariate analysis. *Biometrika*, **53**, 325–338.

Gower, J.C. and Legendre, P. (1986) Metric and Euclidean properties of
dissimilarity coefficients. *Journal of Classification*, **3**, 5–48.

Legendre, P. and Anderson, M.J. (1999) Distance-based redundancy
analysis: Testing multispecies responses in multifactorial ecological
experiments. *Ecological Monographs*, **69**, 1–24.

Legendre, P. and Gallagher, E.D. (2001) Ecologically meaningful
transformations for ordination of species data. *Oecologia*, **129**,
271–280.

Legendre, P. and Legendre, L. (1998) *Numerical Ecology*, 2nd English
edition. Amsterdam: Elsevier Science BV.

Lingoes, J.C. (1971) Some boundary conditions for a monotone analysis of
symmetric matrices. *Psychometrika*, **36**, 195–203.

McArdle, B.H. and Anderson, M.J. (2001) Fitting multivariate models to
community data: A comment on distance-based redundancy analysis.
*Ecology*, **82**, 290–297.

Warton, D.I., Wright, T.W., and Wang, Y. (2012) Distance-based
multivariate analyses confound location and dispersion effects. *Methods
in Ecology and Evolution*, **3**, 89–101.

## Author

John Lennon L. Calorio

## Examples

``` r
## ── Example 1: Bray–Curtis PCoA on the dune dataset ──────────────────────
if (requireNamespace("vegan", quietly = TRUE) &&
    requireNamespace("ape",   quietly = TRUE)) {

  data(dune,     package = "vegan")
  data(dune.env, package = "vegan")

  result <- run_pcoa(dune, method = "bray")
  print(result)
  summary(result)
  head(result$scores[, c("PC1", "PC2")])
  result$variance_explained[1:3]
}
#> Warning: 5 negative eigenvalue(s) detected in the PCoA decomposition.
#>   This is common with non-Euclidean dissimilarity indices (e.g., Bray-Curtis).
#>   Consider correction = "lingoes" or "cailliez" for a fully Euclidean representation.
#> Principal Coordinate Analysis (PCoA)
#> ---------------------------------------- 
#> Dissimilarity method : bray 
#> Correction           : none 
#> Observations         : 20 
#> Axes retained        : 14 
#> 
#> Variance explained (first 5 axes):
#>   PC1: 37.36%
#>   PC2: 22.26%
#>   PC3: 10.05%
#>   PC4: 8.32%
#>   PC5: 6.08%
#> 
#> Principal Coordinate Analysis Summary
#> ============================================= 
#> Call:
#>    run_pcoa(x = dune, method = "bray") 
#> 
#> Dissimilarity method : bray 
#> Correction           : none 
#> Observations (n)     : 20 
#> Total axes           : 14 
#> 
#> Eigenvalues and variance explained:
#>  Axis Eigenvalue Variance_pct Cumulative_pct
#>   PC1   1.716266       37.359         37.359
#>   PC2   1.022398       22.255         59.615
#>   PC3   0.461464       10.045         69.660
#>   PC4   0.382249        8.321         77.980
#>   PC5   0.279135        6.076         84.057
#>   PC6   0.236631        5.151         89.207
#>   PC7   0.169120        3.681         92.889
#>   PC8   0.096245        2.095         94.984
#>   PC9   0.074492        1.622         96.605
#>  PC10   0.061712        1.343         97.949
#>  PC11   0.054940        1.196         99.145
#>  PC12   0.019174        0.417         99.562
#>  PC13   0.016119        0.351         99.913
#>  PC14   0.004001        0.087        100.000
#>      PC1      PC2      PC3 
#> 37.35930 22.25533 10.04505 

## ── Example 2: PCoA + PERMANOVA ───────────────────────────────────────────
if (requireNamespace("vegan", quietly = TRUE) &&
    requireNamespace("ape",   quietly = TRUE)) {

  data(dune,     package = "vegan")
  data(dune.env, package = "vegan")

  result_perm <- run_pcoa(
    x        = dune,
    method   = "bray",
    metadata = dune.env,
    rhs      = "Management",
    perms    = 999
  )

  # Extract PERMANOVA results individually
  result_perm$permanova$R2
  result_perm$permanova$p_value
  result_perm$permanova$table
}
#> Warning: 5 negative eigenvalue(s) detected in the PCoA decomposition.
#>   This is common with non-Euclidean dissimilarity indices (e.g., Bray-Curtis).
#>   Consider correction = "lingoes" or "cailliez" for a fully Euclidean representation.
#> Warning: number of items to replace is not a multiple of replacement length
#>            Df SumOfSqs        R2        F Pr(>F)
#> Management  3 1.468592 0.3416107 2.767243  0.002
#> Residual   16 2.830430 0.6583893       NA     NA
#> Total      19 4.299022 1.0000000       NA     NA

## ── Example 3: Lingoes correction ─────────────────────────────────────────
if (requireNamespace("vegan", quietly = TRUE) &&
    requireNamespace("ape",   quietly = TRUE)) {

  data(dune, package = "vegan")
  result_ling <- run_pcoa(dune, method = "bray", correction = "lingoes")
  result_ling$variance_explained[1:3]
}
#>       PC1       PC2       PC3 
#> 29.538396 18.233837  9.095053 

## ── Example 4: Precomputed dist + PERMANOVA ───────────────────────────────
if (requireNamespace("vegan", quietly = TRUE) &&
    requireNamespace("ape",   quietly = TRUE)) {

  data(dune,     package = "vegan")
  data(dune.env, package = "vegan")
  d <- vegan::vegdist(dune, method = "jaccard")

  result_dist <- run_pcoa(
    x        = d,
    metadata = dune.env,
    rhs      = "Management",
    perms    = 999
  )
  result_dist$method          # "user-supplied dist"
  result_dist$permanova$R2
}
#> Warning: 1 negative eigenvalue(s) detected in the PCoA decomposition.
#>   This is common with non-Euclidean dissimilarity indices (e.g., Bray-Curtis).
#>   Consider correction = "lingoes" or "cailliez" for a fully Euclidean representation.
#> Warning: number of items to replace is not a multiple of replacement length
#> [1] 0.2943778
```
