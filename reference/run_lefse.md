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

The input data should ideally be relative abundances (e.g., from
`run_relabund` or similar). The function assumes that higher values
indicate higher abundance.

## See also

`plot_lefse` for visualization (if implemented),
[`run_assoc`](https://jllcalorio.github.io/pondeR/reference/run_assoc.md)
for association testing.

## Examples

``` r
# Load example data (if available) or use your own
# data(metagenomicexample, package = "pondeR")
if (FALSE) { # \dontrun{
  # Assuming `otu_table` is a feature table and `meta` is metadata
  lefse_res <- run_lefse(
    x = otu_table,
    metadata = meta,
    class = "Condition",
    subclass = "AgeGroup",
    alpha = 0.05,
    wilcox_alpha = 0.05,
    lda_threshold = 2.0
  )
} # }
```
