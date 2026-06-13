# Run Lawshe's Content Validity Ratio (CVR) Analysis

Computes Lawshe's (1975) Content Validity Ratio (CVR) for one or more
items rated by an expert panel. Each expert classifies an item as
"essential", "useful but not essential", or "not necessary"; CVR
captures the proportion of experts endorsing "essential" relative to
chance. The function returns a tidy data frame of CVR values, critical
values, and significance interpretations for each item.

## Usage

``` r
run_contentvaliditycvr(
  x,
  essential = 1,
  n = dim(x)[1L],
  alpha = 0.05,
  na.rm = FALSE
)
```

## Arguments

- x:

  A numeric data frame or matrix where rows represent experts and
  columns represent items.

- essential:

  A numeric vector indicating the rating value(s) that classify an item
  as "essential". Defaults to `1`. A vector may be supplied (e.g.,
  `c(1, 2)`) when multiple adjacent categories should count as
  essential.

- n:

  An integer specifying the panel size used to compute the minimum CVR
  considered statistically significant via `cvr_critical()`. Defaults to
  `dim(x)[1]` (the number of rows in `x`).

- alpha:

  A numeric value for the one-tailed significance level used in
  `cvr_critical()`. Defaults to `0.05`.

- na.rm:

  Logical. If `TRUE`, missing ratings are excluded when counting experts
  per item. Defaults to `FALSE`.

## Value

A data frame of class `run_contentvaliditycvr` with one row per item and
the following columns:

- `CVR`:

  Lawshe's Content Validity Ratio.

- `Critical Value`:

  Minimum CVR for statistical significance at the specified `alpha` and
  panel size `n`.

- `Interpretation`:

  Character; `"Significant"` if `CVR >= Critical Value`, otherwise
  `"Not significant"`.

## Details

CVR is defined as: \$\$\text{CVR} = \frac{n_e - n/2}{n/2}\$\$ where
\\n_e\\ is the number of experts rating the item as "essential" and
\\n\\ is the total number of experts. CVR ranges from \\-1\\ (no expert
endorses "essential") to \\+1\\ (all experts endorse "essential"). A
value of \\0\\ indicates that exactly half the panel rated the item
"essential".

Statistical significance of each CVR is evaluated against the critical
value returned by `contentValidity::cvr_critical(n, alpha)`, which
implements the recalculated critical values of Wilson, Pan, and Schumsky
(2012).

## References

Lawshe, C. H. (1975). A quantitative approach to content validity.
*Personnel Psychology*, *28*(4), 563–575.
[doi:10.1111/j.1744-6570.1975.tb01393.x](https://doi.org/10.1111/j.1744-6570.1975.tb01393.x)

Wilson, F. R., Pan, W., & Schumsky, D. A. (2012). Recalculation of the
critical values for Lawshe's content validity ratio. *Measurement and
Evaluation in Counseling and Development*, *45*(3), 197–210.
[doi:10.1177/0748175612440286](https://doi.org/10.1177/0748175612440286)

## See also

[`cvr`](https://rdrr.io/pkg/contentValidity/man/cvr.html) for the
underlying CVR computation,
[`cvr_critical`](https://rdrr.io/pkg/contentValidity/man/cvr_critical.html)
for the critical value lookup.

## Examples

``` r
if (FALSE) { # \dontrun{
# Lawshe's (1975) original 3-point scale:
# 1 = essential, 2 = useful but not essential, 3 = not necessary
ratings <- data.frame(
  item1 = c(1, 1, 1, 2, 1),
  item2 = c(3, 3, 1, 3, 3),
  item3 = c(1, 1, 1, 1, 1)
)
result <- run_contentvaliditycvr(ratings, essential = 1)
result

# Treat both "essential" and "useful but not essential" as endorsing the item
result2 <- run_contentvaliditycvr(ratings, essential = c(1, 2))
result2
} # }
```
