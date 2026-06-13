# Run Comprehensive Item-Level and Scale-Level Content Validity Analysis

Computes standard relevance-based content validity indices on a single
ratings matrix and returns a tidy summary. At the item level, the
function calculates the Item-level Content Validity Index (I-CVI),
modified kappa, and Aiken's V. At the scale level, it reports the
Scale-level CVI by averaging (S-CVI/Ave), Scale-level CVI by universal
agreement (S-CVI/UA), and the mean modified kappa across all items.

## Usage

``` r
run_contentvalidity(x, threshold = 3L, low = 1L, high = 4L, na.rm = FALSE)
```

## Arguments

- x:

  A numeric data frame or matrix where rows represent experts and
  columns represent items.

- threshold:

  An integer specifying the minimum rating value considered "relevant".
  Defaults to `3`.

- low:

  An integer specifying the minimum possible rating value on the scale.
  Passed to `aiken_v()`. Defaults to `1`.

- high:

  An integer specifying the maximum possible rating value on the scale.
  Passed to `aiken_v()`. Defaults to `4`.

- na.rm:

  Logical. If `TRUE`, ratings with missing values are removed prior to
  computation. Defaults to `FALSE`.

## Value

A list of class `run_contentvalidity` with two elements:

- `$items`:

  A data frame of item-level indices: I-CVI, modified kappa, Aiken's V,
  and Modified-Kappa Interpretation.

- `$scale`:

  A named numeric vector of scale-level indices: S-CVI/Ave, S-CVI/UA,
  and mean modified kappa.

## Details

The function merges the unrounded numeric output of
[`contentValidity::content_validity()`](https://rdrr.io/pkg/contentValidity/man/content_validity.html)
with the formatted column names and interpretation labels from
[`contentValidity::apa_table()`](https://rdrr.io/pkg/contentValidity/man/apa_table.html),
yielding full-precision estimates alongside human-readable column
headers.

**Modified-Kappa Cutoffs**

Modified kappa values are interpreted using the criteria of Cicchetti
and Sparrow (1981), as adopted by Polit, Beck, and Owen (2007):

|                              |                    |
|------------------------------|--------------------|
| **Range**                    | **Interpretation** |
| \\\kappa \> 0.74\\           | Excellent          |
| \\0.60 \le \kappa \le 0.74\\ | Good               |
| \\0.40 \le \kappa \le 0.59\\ | Fair               |
| \\\kappa \< 0.40\\           | Poor               |

## References

Aiken, L. R. (1985). Three coefficients for analyzing the reliability
and validity of ratings. *Educational and Psychological Measurement*,
*45*(1), 131–142.
[doi:10.1177/0013164485451012](https://doi.org/10.1177/0013164485451012)

Cicchetti, D. V., & Sparrow, S. A. (1981). Developing criteria for
establishing interrater reliability of specific items: Applications to
assessment of adaptive behavior. *American Journal of Mental
Deficiency*, *86*(2), 127–137.

Lynn, M. R. (1986). Determination and quantification of content
validity. *Nursing Research*, *35*(6), 382–385.
[doi:10.1097/00006199-198611000-00017](https://doi.org/10.1097/00006199-198611000-00017)

Polit, D. F., & Beck, C. T. (2006). The content validity index: Are you
sure you know what's being reported? Critique and recommendations.
*Research in Nursing & Health*, *29*(5), 489–497.
[doi:10.1002/nur.20147](https://doi.org/10.1002/nur.20147)

Polit, D. F., Beck, C. T., & Owen, S. V. (2007). Is the CVI an
acceptable indicator of content validity? Appraisal and recommendations.
*Research in Nursing & Health*, *30*(4), 459–467.
[doi:10.1002/nur.20199](https://doi.org/10.1002/nur.20199)

## See also

[`icvi`](https://rdrr.io/pkg/contentValidity/man/icvi.html),
[`scvi_ave`](https://rdrr.io/pkg/contentValidity/man/scvi_ave.html),
[`scvi_ua`](https://rdrr.io/pkg/contentValidity/man/scvi_ua.html),
[`mod_kappa`](https://rdrr.io/pkg/contentValidity/man/mod_kappa.html),
[`aiken_v`](https://rdrr.io/pkg/contentValidity/man/aiken_v.html)

## Examples

``` r
if (FALSE) { # \dontrun{
ratings <- data.frame(
  item1 = c(4, 3, 4, 4, 3),
  item2 = c(2, 3, 2, 1, 2),
  item3 = c(4, 4, 3, 4, 4)
)
result <- run_contentvalidity(ratings)
result$items
result$scale
} # }
```
