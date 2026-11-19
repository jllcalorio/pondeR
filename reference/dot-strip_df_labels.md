# Strip "name = " labels from a .format_df() string and reformat the numbers per n_digits / zero_as_exp, e.g. "df1 = 2, df2 = 27" -\> "2, 27", or "num df = 2, denom df = 24.316" -\> "2, 24.32" (n_digits = 2). Needed because non-integer df (e.g. Welch's ANOVA denom df) otherwise keeps .format_df()'s fixed 4-decimal rounding regardless of n_digits.

Strip "name = " labels from a .format_df() string and reformat the
numbers per n_digits / zero_as_exp, e.g. "df1 = 2, df2 = 27" -\> "2,
27", or "num df = 2, denom df = 24.316" -\> "2, 24.32" (n_digits = 2).
Needed because non-integer df (e.g. Welch's ANOVA denom df) otherwise
keeps .format_df()'s fixed 4-decimal rounding regardless of n_digits.

## Usage

``` r
.strip_df_labels(x, dig = 2, zero_as_exp = TRUE)
```
