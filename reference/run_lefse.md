# Perform LEfSe (Linear Discriminant Analysis Effect Size) analysis

This function performs LEfSe analysis to identify differentially
abundant features between classes (e.g., biomarkers) with effect size
estimation.

## Usage

``` r
run_lefse(
  x,
  metadata,
  class = "Group",
  subclass = NULL,
  alpha = 0.05,
  wilcox_alpha = 0.05,
  lda_threshold = 2,
  verbose = TRUE
)
```

## Arguments

- x:

  A data frame or matrix of features (samples x features). Should be
  relative abundance data (proportions summing to 1 per sample).

- metadata:

  A data frame containing sample metadata. Must include the column
  specified by `class` and optionally `subclass`.

- class:

  A single character string specifying the column name in `metadata`
  that defines the main classes for comparison (e.g., disease status).

- subclass:

  Optional. A single character string specifying the column name in
  `metadata` that defines subclasses (e.g., age groups, gender) for
  stratified analysis. If provided, the analysis is performed within
  each subclass.

- alpha:

  Numeric value between 0 and 1 specifying the significance level for
  the Kruskal-Wallis test (default: 0.05).

- wilcox_alpha:

  Numeric value between 0 and 1 specifying the significance level for
  the pairwise Wilcoxon test (default: 0.05).

- lda_threshold:

  Numeric value specifying the threshold on the absolute logarithmic LDA
  score for discriminative features (default: 2.0).

- verbose:

  Logical indicating whether to print progress messages (default: TRUE).

## Value

A list of class `run_lefse` containing:

- lda:

  A data frame with LDA results (features, means per class, LDA scores).

- summary_table:

  A data frame summarizing significant features with effect sizes,
  p-values, and classifications.

- parameters:

  A list of the parameters used in the analysis.

- call:

  The matched call.

## Details

The LEfSe algorithm consists of three steps:

1.  Kruskal-Wallis test to detect significantly abundant features
    between classes.

2.  Pairwise Wilcoxon test to check consistency of biological trends
    across subclasses within each class.

3.  Linear Discriminant Analysis to estimate the effect size of each
    significantly abundant feature.

The LDA score of each feature is \\\log\_{10}\\1 + (\|\Delta\bar{x}\| +
\|w \cdot \Delta LD\|)/2\\\\, where \\\Delta\bar{x}\\ is the largest
difference in class means, \\w\\ the unit-length first
linear-discriminant coefficient
([`lda`](https://rdrr.io/pkg/MASS/man/lda.html)) and \\\Delta LD\\ the
class separation along that discriminant (Segata et al., 2011). Unlike
the original LEfSe, the LDA is fitted once (no bootstrap).

Scores depend on the scale of `x`. Original LEfSe scales each sample's
relative abundances to sum to \\10^6\\; do the same (e.g.
`x / rowSums(x) * 1e6`) for the default `lda_threshold = 2` to be
meaningful. Higher values are assumed to indicate higher abundance.

## References

Segata N, Izard J, Waldron L, et al. (2011). Metagenomic biomarker
discovery and explanation. *Genome Biology*, 12(6), R60.
[doi:10.1186/gb-2011-12-6-r60](https://doi.org/10.1186/gb-2011-12-6-r60)

Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance
analysis. *Journal of the American Statistical Association*, 47(260),
583-621.
[doi:10.1080/01621459.1952.10483441](https://doi.org/10.1080/01621459.1952.10483441)

Mann HB, Whitney DR (1947). On a test of whether one of two random
variables is stochastically larger than the other. *The Annals of
Mathematical Statistics*, 18(1), 50-60.
[doi:10.1214/aoms/1177730491](https://doi.org/10.1214/aoms/1177730491)

## Examples

``` r
set.seed(1)
x <- as.data.frame(matrix(rpois(40 * 6, 20), 40,
                          dimnames = list(NULL, paste0("taxon", 1:6))))
meta <- data.frame(Group = rep(c("Control", "Case"), each = 20))
x$taxon1[meta$Group == "Case"] <- x$taxon1[meta$Group == "Case"] + 30

rel <- x / rowSums(x) * 1e6   # LEfSe scale
res <- run_lefse(rel, meta, class = "Group", verbose = FALSE)
#> Warning: variables are collinear
res$summary_table
#>   feature mean_Control mean_Case lda_score       pvalue
#> 1  taxon1     176427.2  338688.4  5.206063 6.292924e-08
#> 3  taxon3     177038.1  125748.7  4.686542 2.750987e-05
#> 5  taxon5     159957.5  125765.5  4.578816 1.187308e-02
#> 6  taxon6     158882.2  131631.8  4.450627 9.402764e-03
#> 2  taxon2     168136.0  144647.1  4.190324 2.844022e-02
#> 4  taxon4     159559.0  133518.4  4.184950 8.024598e-03
```
