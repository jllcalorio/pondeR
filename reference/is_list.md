# Check if an Object is a Plain List

Tests whether an object is a base R `list` (but not a `data.frame`,
tibble, or `data.table`, which are technically lists under the hood).
Returns a single logical value (`TRUE` or `FALSE`).

## Usage

``` r
is_list(x, warn = TRUE)
```

## Arguments

- x:

  An R object to test.

- warn:

  Logical scalar. If `TRUE` (default), an informative warning is emitted
  when `x` fails the check. Set to `FALSE` to suppress all warnings and
  obtain a silent `FALSE`.

## Value

A single logical value: `TRUE` if `x` is a plain list, `FALSE`
otherwise. Never `NA`.

## Details

`is_list()` returns `TRUE` only when `x` is a plain `list` - i.e.,
`is.list(x)` is `TRUE` *and* `x` does not inherit from any tabular class
(`data.frame`, `tbl_df`, `tbl`, `data.table`, or `matrix`). This
distinction matters because `data.frame` and `data.table` objects are
internally lists but should not be treated as such in the context of
pondeR workflows.

`NULL` always returns `FALSE`.

## Author

John Lennon L. Calorio

## Examples

``` r
# TRUE cases
is_list(list(a = 1, b = "hello"))
#> [1] TRUE
is_list(list(1:3, TRUE, NULL))
#> [1] TRUE

# FALSE cases - tabular objects are excluded
is_list(mtcars)
#> Warning: 'x' has class <data.frame>, which is a tabular object, not a plain list. Use is_tabular() to check for data.frame, tibble, data.table, or matrix objects.
#> [1] FALSE
is_list(matrix(1:4, 2))
#> Warning: 'x' has class <matrix, array>, which is a tabular object, not a plain list. Use is_tabular() to check for data.frame, tibble, data.table, or matrix objects.
#> [1] FALSE

# FALSE - atomic vector
is_list(1:10)
#> Warning: 'x' has class <integer>, which is not a list. A plain list object was expected.
#> [1] FALSE

# FALSE - NULL
is_list(NULL)
#> Warning: 'x' is NULL. A plain list was expected.
#> [1] FALSE

# Silent FALSE for programmatic use
is_list(42, warn = FALSE)
#> [1] FALSE
```
