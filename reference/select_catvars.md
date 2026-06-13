# Select Categorical Variables from a Data Frame

This function selects variables from a data frame that are categorical
(factors, ordered factors, characters, or logicals). It also allows
forcing specific columns to be treated as categorical even if they are
stored as numeric types.

## Usage

``` r
select_catvars(x, force_cat = NULL, include_chr = TRUE, include_lgl = FALSE)
```

## Arguments

- x:

  A data frame or matrix.

- force_cat:

  Character vector of column names in `x` to be forcefully included as
  categorical variables.

- include_chr:

  Logical. If `TRUE` (default), includes character type variables.

- include_lgl:

  Logical. If `FALSE` (default), excludes logical type variables.

## Value

A data frame containing only the selected categorical variables.

## Author

John Lennon L. Calorio

## Examples

``` r
df <- data.frame(a = factor(c("A", "B")), b = 1:2, 
c = c("x", "y"),
d = c(TRUE, FALSE))
select_catvars(df)
#>   a c
#> 1 A x
#> 2 B y
```
