# Recode and Factorize a Column

This function recodes values in a specified column of a data frame
according to a provided mapping, then converts the column to a factor
with ordered levels. The recoded values appear first in the factor
levels, followed by any remaining unique values in alphabetical order.

## Usage

``` r
run_relabel(x, column, conversion_map)
```

## Arguments

- x:

  A data frame containing the column to be recoded and factorized.

- column:

  A character string specifying the name of the column to recode. The
  column must exist in the data frame.

- conversion_map:

  A named list or vector where names are the original values to be
  replaced and values are the new values. All elements must be atomic
  (character, numeric, or logical).

## Value

A data frame identical to the input except that the specified column has
been recoded according to the conversion map and converted to a factor
with ordered levels.

## Details

The function performs the following steps:

1.  Converts the specified column to character to handle mixed data
    types

2.  Applies the recoding based on the conversion map

3.  Creates factor levels with recoded values first, followed by other
    unique values in alphabetical order

4.  Converts the column to a factor with the specified level ordering

## Examples

``` r
# Basic usage with character recoding
x <- data.frame(
  id = 1:5,
  status = c("A", "B", "C", "A", "D")
)
map <- c("A" = "Active", "B" = "Inactive", "C" = "Pending")
result <- run_relabel(x, "status", map)
levels(result$status)  # "Active" "Inactive" "Pending" "D"
#> [1] "Active"   "Inactive" "Pending"  "D"       

# Usage with numeric to character conversion
x2 <- data.frame(
  score = c(1, 2, 3, 1, 4),
  grade = c(90, 85, 78, 92, 65)
)
grade_map <- c("1" = "Excellent", "2" = "Good", "3" = "Fair")
result2 <- run_relabel(x2, "score", grade_map)
```
