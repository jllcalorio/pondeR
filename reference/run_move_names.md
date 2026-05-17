# Move Row or Column Names into the Data Area

This function moves existing row names into a specified column or column
names into a specified row. By default, it prioritizes moving row names
to the first column. If both a column index and a row index are
provided, row name movement takes precedence and a warning is issued.

## Usage

``` r
run_move_names(
  x,
  col_index = 1,
  row_index = NULL,
  col_name = "row_names",
  remove = TRUE
)
```

## Arguments

- x:

  A dataframe, tibble, or matrix.

- col_index:

  Numeric. The column index where to place the row names. Defaults to 1.

- row_index:

  Numeric or NULL. The row index where to place the column names. If not
  NULL and `col_index` is also active, this will be ignored with a
  warning.

- col_name:

  String. The name of the new column created when moving row names.
  Defaults to "row_names".

- remove:

  Logical. If TRUE (default), removes the original row/column names
  after moving them into the data.

## Value

A data frame with the names moved into the data area. When moving row
names to a column, original data types of other columns are preserved.

## Author

John Lennon L. Calorio

## Examples

``` r
# 1. Move row names to the first column (default)
run_move_names(mtcars[1:5, 1:3])
#>           row_names  mpg cyl disp
#> 1         Mazda RX4 21.0   6  160
#> 2     Mazda RX4 Wag 21.0   6  160
#> 3        Datsun 710 22.8   4  108
#> 4    Hornet 4 Drive 21.4   6  258
#> 5 Hornet Sportabout 18.7   8  360

# 2. Move row names to a specific column index
run_move_names(mtcars[1:5, 1:3], col_index = 2)
#>    mpg         row_names cyl disp
#> 1 21.0         Mazda RX4   6  160
#> 2 21.0     Mazda RX4 Wag   6  160
#> 3 22.8        Datsun 710   4  108
#> 4 21.4    Hornet 4 Drive   6  258
#> 5 18.7 Hornet Sportabout   8  360

# 3. Move column names to the first row
# Note: col_index must be NULL to process row_index per requirements
run_move_names(mtcars[1:5, 1:3], col_index = NULL, row_index = 1)
#>                     V1  V2   V3
#> cn                 mpg cyl disp
#> Mazda RX4           21   6  160
#> Mazda RX4 Wag       21   6  160
#> Datsun 710        22.8   4  108
#> Hornet 4 Drive    21.4   6  258
#> Hornet Sportabout 18.7   8  360

# 4. Demonstrate priority warning
run_move_names(mtcars[1:5, 1:3], col_index = 1, row_index = 1)
#> Warning: Ambiguity detected: Both 'col_index' and 'row_index' provided. Prioritizing movement of row names to column. 'row_index' will be ignored.
#>           row_names  mpg cyl disp
#> 1         Mazda RX4 21.0   6  160
#> 2     Mazda RX4 Wag 21.0   6  160
#> 3        Datsun 710 22.8   4  108
#> 4    Hornet 4 Drive 21.4   6  258
#> 5 Hornet Sportabout 18.7   8  360
```
