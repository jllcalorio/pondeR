# Format a single numeric value, falling back to scientific notation when it would otherwise round to zero (mirrors run_summarytable's zero_as_exp). Whole numbers (e.g. df counts, exact-0/1 CI bounds) are shown bare, with no decimal places, regardless of `dig`.

Format a single numeric value, falling back to scientific notation when
it would otherwise round to zero (mirrors run_summarytable's
zero_as_exp). Whole numbers (e.g. df counts, exact-0/1 CI bounds) are
shown bare, with no decimal places, regardless of `dig`.

## Usage

``` r
.format_num_simplify(v, dig = 2, zero_as_exp = TRUE)
```
