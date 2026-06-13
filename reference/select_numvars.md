# Select Numeric Variables from a Data Frame

This function selects variables from a data frame based on their numeric
types (integer, double, logical, or complex). It also allows forcing
specific columns to be treated as numeric even if they are stored as
other types.

## Usage

``` r
select_numvars(
  x,
  force_num = NULL,
  include_int = TRUE,
  include_dbl = TRUE,
  include_lgl = TRUE,
  include_comp = TRUE
)
```

## Arguments

- x:

  A data frame or matrix.

- force_num:

  Character vector of column names in `x` to be forcefully included as
  numeric variables.

- include_int:

  Logical. If `TRUE` (default), includes integer type variables.

- include_dbl:

  Logical. If `TRUE` (default), includes double type variables.

- include_lgl:

  Logical. If `TRUE` (default), includes logical type variables.

- include_comp:

  Logical. If `TRUE` (default), includes complex type variables.

## Value

A data frame containing only the selected numeric variables.

## Author

John Lennon L. Calorio

## Examples

``` r
df <- data.frame(a = 1:3, b = c(1.1, 2.2, 3.3), 
c = c("x", "y", "z"),
d = c(TRUE, FALSE, TRUE))
select_numvars(df)
#>   a   b     d
#> 1 1 1.1  TRUE
#> 2 2 2.2 FALSE
#> 3 3 3.3  TRUE
```
