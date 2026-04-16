# Randomize a Data Frame, Matrix, or Vector

Randomly samples rows from a data frame, matrix, or elements from a
vector, with optional group-stratified sampling. Serves as a structured
wrapper around
[`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html),
extending it with flexible per-group control (joint or
sequential/independent stratification), reproducible seeding, and strict
input validation suitable for biological and clinical research
workflows.

## Usage

``` r
run_randomize(
  x,
  n = NULL,
  p = NULL,
  per_group = NULL,
  joint = TRUE,
  joint_behaviour = "cascading",
  replace = FALSE,
  seed = 123L
)
```

## Arguments

- x:

  A data frame, tibble, matrix, or atomic vector. Column names with
  special characters are supported. Matrices are coerced to data frames
  internally and the result is recast to a matrix on return.

- n:

  A single positive integer specifying the number of rows (or elements,
  if `x` is a vector) to return. Defaults to `NULL`. Mutually exclusive
  with `p`; if both are supplied, `n` takes priority and a warning is
  issued. When `per_group` is supplied without an explicit `n`/`p`
  inside the list, this top-level value is used as the fallback.

- p:

  A single numeric value in `[0, 1]` specifying the proportion of rows
  (or elements) to return. Defaults to `NULL` (equivalent to `p = 1`,
  i.e., all rows shuffled). Mutually exclusive with `n`. When
  `per_group` is supplied without an explicit `n`/`p` inside the list,
  this top-level value is used as the fallback.

- per_group:

  An optional list of per-group sampling specs controlling stratified
  sampling. Each element must follow one of two forms:

  - **Single column**: `list("ColumnName", n = <integer>)` or
    `list("ColumnName", p = <number>)`

  - **Multiple columns (joint or independent)**:  
    `list(c("Col1", "Col2"), n = <integer>)` or  
    `list(c("Col1", "Col2"), p = <number>)` or  
    `list(c("Col1", "Col2"), n = c(2, 3))` (per-column values, for
    `joint = FALSE` + `joint_behaviour = "independent"` only)

  The `n`/`p` inside each list element is optional. If omitted, the
  top-level `n` or `p` argument is used as the fallback (top-level `n`
  takes priority over top-level `p` if both are provided). If both `n`
  and `p` are supplied inside a list element, `n` takes priority and a
  warning is issued.

  Requires `x` to be a data frame or matrix. Grouping columns must be
  categorical (character, factor, or logical). See `joint` and
  `joint_behaviour` for how multiple columns are handled.

- joint:

  A single logical value. Defaults to `TRUE`.

  - `TRUE`: all columns listed in a `per_group` entry are treated
    *simultaneously* — rows are grouped by every unique combination of
    levels across those columns (via
    [`interaction`](https://rdrr.io/r/base/interaction.html)), and
    `n`/`p` is enforced for each combination. For example, with `cyl`
    and `gear`, a temporary group key such as `"6:4"` or `"8:3"` is
    created internally and discarded after sampling.

  - `FALSE`: the `joint_behaviour` argument controls whether columns are
    handled in *cascading* or *independent* fashion.

  Ignored when only a single column is listed in a `per_group` entry.

- joint_behaviour:

  A character string, either `"cascading"` (or `"c"`) or `"independent"`
  (or `"i"`). Only applies when `joint = FALSE` and a `per_group` entry
  lists multiple columns. Defaults to `"cascading"`.

  - `"cascading"`: columns are applied sequentially in the order listed.
    The data are first grouped by column 1 and sampled; then the
    *reduced result* is grouped by column 2 and sampled again, and so
    on. Each step's `n`/`p` is applied to whatever rows remain after the
    previous step, which can aggressively reduce the output when `n` is
    small.

  - `"independent"`: each column is sampled independently from the
    *original* `x`. The results are combined by taking the *union* of
    selected row indices, so a row qualifies if it was chosen under
    *any* of the grouping columns. When per-column `n`/`p` vectors are
    supplied (e.g., `n = c(2, 3)`), each value corresponds positionally
    to each column. Vectors must not exceed the number of columns in
    length.

- replace:

  A single logical value. If `FALSE` (default), sampling is performed
  without replacement. If `TRUE`, rows or elements may be repeated in
  the output.

- seed:

  A single integer passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) for reproducibility.
  Defaults to `123L`. Supply `NA` to skip seeding (non-reproducible
  runs).

## Value

An object of the same class as `x`:

- **Data frame / tibble**: a data frame or tibble with sampled rows,
  preserving all columns and their classes. Row names are dropped
  (consistent with
  [`slice_sample`](https://dplyr.tidyverse.org/reference/slice.html)
  behaviour).

- **Matrix**: a matrix with sampled rows, preserving column names and
  storage mode.

- **Vector**: an atomic vector of sampled elements, preserving names if
  present.

## Details

### Sampling modes

`run_randomize()` operates in three mutually exclusive top-level modes,
resolved in the following order of priority:

1.  **Per-group stratified sampling** (when `per_group` is not `NULL`):
    Rows are sampled within groups defined by the column(s) in each
    `per_group` entry. The `joint` and `joint_behaviour` arguments
    control how multiple columns are handled (see below).

2.  **Fixed-count sampling** (`n` supplied, `per_group = NULL`): Exactly
    `n` rows are returned. Setting `n` equal to `nrow(x)` returns all
    rows in a new random order.

3.  **Proportional sampling** (`p` supplied, or both `n` and `p`
    omitted): A fraction `p` of rows is returned. Omitting both defaults
    to `p = 1` (full shuffle). `p = 0` returns a zero-row object of the
    same class.

### Per-group syntax

Each element of `per_group` bundles one or more grouping columns with an
optional sampling quantity:

    # Single column, explicit n
    per_group = list(
      list("Species", n = 10)
    )

    # Two columns, joint (default), explicit n shared across all combinations
    per_group = list(
      list(c("cyl", "gear"), n = 2)
    )

    # Two columns, independent, per-column n vector
    per_group = list(
      list(c("cyl", "gear"), n = c(2, 3))
    )

    # Omit n/p to fall back on the top-level n or p argument
    per_group = list(
      list("Species")
    )

### Behaviour of `joint` and `joint_behaviour`

|           |                     |                                                                                                                                        |
|-----------|---------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| **joint** | **joint_behaviour** | **Effect**                                                                                                                             |
| `TRUE`    | (ignored)           | Groups by all columns simultaneously via [`interaction()`](https://rdrr.io/r/base/interaction.html); one `n`/`p` for every combination |
| `FALSE`   | `"cascading"`       | Applies each column's `n`/`p` sequentially on the progressively reduced data                                                           |
| `FALSE`   | `"independent"`     | Samples each column independently from the original data; returns the union of selected rows                                           |

### Output row order

When `per_group` is used, the final output is randomly shuffled so that
rows from the same group are not clustered together. For example, if
grouping by `Species`, the returned rows will be intermixed across
species rather than appearing as contiguous blocks of *setosa*, then
*versicolor*, and so on. This mimics a fully randomized layout and
avoids order effects in downstream analyses.

## Author

John Lennon L. Calorio

## Examples

``` r
## ------------------------------------------------------------------
## 1. Shuffle all rows (default: p = 1)
## ------------------------------------------------------------------
run_randomize(mtcars)
#>                      mpg cyl  disp  hp drat    wt  qsec vs am gear carb
#> Maserati Bora       15.0   8 301.0 335 3.54 3.570 14.60  0  1    5    8
#> Cadillac Fleetwood  10.4   8 472.0 205 2.93 5.250 17.98  0  0    3    4
#> Honda Civic         30.4   4  75.7  52 4.93 1.615 18.52  1  1    4    2
#> Merc 450SLC         15.2   8 275.8 180 3.07 3.780 18.00  0  0    3    3
#> Datsun 710          22.8   4 108.0  93 3.85 2.320 18.61  1  1    4    1
#> Merc 280            19.2   6 167.6 123 3.92 3.440 18.30  1  0    4    4
#> Fiat 128            32.4   4  78.7  66 4.08 2.200 19.47  1  1    4    1
#> Dodge Challenger    15.5   8 318.0 150 2.76 3.520 16.87  0  0    3    2
#> Merc 280C           17.8   6 167.6 123 3.92 3.440 18.90  1  0    4    4
#> Hornet Sportabout   18.7   8 360.0 175 3.15 3.440 17.02  0  0    3    2
#> Toyota Corolla      33.9   4  71.1  65 4.22 1.835 19.90  1  1    4    1
#> Ford Pantera L      15.8   8 351.0 264 4.22 3.170 14.50  0  1    5    4
#> AMC Javelin         15.2   8 304.0 150 3.15 3.435 17.30  0  0    3    2
#> Ferrari Dino        19.7   6 145.0 175 3.62 2.770 15.50  0  1    5    6
#> Merc 230            22.8   4 140.8  95 3.92 3.150 22.90  1  0    4    2
#> Lotus Europa        30.4   4  95.1 113 3.77 1.513 16.90  1  1    5    2
#> Merc 240D           24.4   4 146.7  62 3.69 3.190 20.00  1  0    4    2
#> Porsche 914-2       26.0   4 120.3  91 4.43 2.140 16.70  0  1    5    2
#> Duster 360          14.3   8 360.0 245 3.21 3.570 15.84  0  0    3    4
#> Volvo 142E          21.4   4 121.0 109 4.11 2.780 18.60  1  1    4    2
#> Fiat X1-9           27.3   4  79.0  66 4.08 1.935 18.90  1  1    4    1
#> Chrysler Imperial   14.7   8 440.0 230 3.23 5.345 17.42  0  0    3    4
#> Hornet 4 Drive      21.4   6 258.0 110 3.08 3.215 19.44  1  0    3    1
#> Mazda RX4           21.0   6 160.0 110 3.90 2.620 16.46  0  1    4    4
#> Camaro Z28          13.3   8 350.0 245 3.73 3.840 15.41  0  0    3    4
#> Toyota Corona       21.5   4 120.1  97 3.70 2.465 20.01  1  0    3    1
#> Pontiac Firebird    19.2   8 400.0 175 3.08 3.845 17.05  0  0    3    2
#> Merc 450SL          17.3   8 275.8 180 3.07 3.730 17.60  0  0    3    3
#> Lincoln Continental 10.4   8 460.0 215 3.00 5.424 17.82  0  0    3    4
#> Mazda RX4 Wag       21.0   6 160.0 110 3.90 2.875 17.02  0  1    4    4
#> Merc 450SE          16.4   8 275.8 180 3.07 4.070 17.40  0  0    3    3
#> Valiant             18.1   6 225.0 105 2.76 3.460 20.22  1  0    3    1

## ------------------------------------------------------------------
## 2. Sample a fixed number of rows
## ------------------------------------------------------------------
run_randomize(mtcars, n = 5)
#>                     mpg cyl  disp  hp drat    wt  qsec vs am gear carb
#> Maserati Bora      15.0   8 301.0 335 3.54 3.570 14.60  0  1    5    8
#> Cadillac Fleetwood 10.4   8 472.0 205 2.93 5.250 17.98  0  0    3    4
#> Honda Civic        30.4   4  75.7  52 4.93 1.615 18.52  1  1    4    2
#> Merc 450SLC        15.2   8 275.8 180 3.07 3.780 18.00  0  0    3    3
#> Datsun 710         22.8   4 108.0  93 3.85 2.320 18.61  1  1    4    1

## ------------------------------------------------------------------
## 3. Sample a proportion of rows
## ------------------------------------------------------------------
run_randomize(mtcars, p = 0.25)
#>                     mpg cyl  disp  hp drat    wt  qsec vs am gear carb
#> Maserati Bora      15.0   8 301.0 335 3.54 3.570 14.60  0  1    5    8
#> Cadillac Fleetwood 10.4   8 472.0 205 2.93 5.250 17.98  0  0    3    4
#> Honda Civic        30.4   4  75.7  52 4.93 1.615 18.52  1  1    4    2
#> Merc 450SLC        15.2   8 275.8 180 3.07 3.780 18.00  0  0    3    3
#> Datsun 710         22.8   4 108.0  93 3.85 2.320 18.61  1  1    4    1
#> Merc 280           19.2   6 167.6 123 3.92 3.440 18.30  1  0    4    4
#> Fiat 128           32.4   4  78.7  66 4.08 2.200 19.47  1  1    4    1
#> Dodge Challenger   15.5   8 318.0 150 2.76 3.520 16.87  0  0    3    2

## ------------------------------------------------------------------
## 4. Sampling with replacement
## ------------------------------------------------------------------
run_randomize(mtcars, n = 10, replace = TRUE)
#>                     mpg cyl  disp  hp drat    wt  qsec vs am gear carb
#> Maserati Bora      15.0   8 301.0 335 3.54 3.570 14.60  0  1    5    8
#> Cadillac Fleetwood 10.4   8 472.0 205 2.93 5.250 17.98  0  0    3    4
#> Honda Civic        30.4   4  75.7  52 4.93 1.615 18.52  1  1    4    2
#> Merc 450SLC        15.2   8 275.8 180 3.07 3.780 18.00  0  0    3    3
#> Datsun 710         22.8   4 108.0  93 3.85 2.320 18.61  1  1    4    1
#> Merc 280           19.2   6 167.6 123 3.92 3.440 18.30  1  0    4    4
#> Fiat 128           32.4   4  78.7  66 4.08 2.200 19.47  1  1    4    1
#> Dodge Challenger   15.5   8 318.0 150 2.76 3.520 16.87  0  0    3    2
#> Merc 280C          17.8   6 167.6 123 3.92 3.440 18.90  1  0    4    4
#> Hornet Sportabout  18.7   8 360.0 175 3.15 3.440 17.02  0  0    3    2

## ------------------------------------------------------------------
## 5. Randomize a vector
## ------------------------------------------------------------------
run_randomize(1:20, n = 8, seed = 42)
#> [1] 17  5  1 10  4  2 20 18

run_randomize(c(a = 1, b = 2, c = 3, d = 4), p = 0.5, seed = 7)
#> b c 
#> 2 3 

## ------------------------------------------------------------------
## 6. Per-group: single column, explicit n
## ------------------------------------------------------------------
run_randomize(
  iris,
  per_group = list(list("Species", n = 10)),
  seed = 2024
)
#> # A tibble: 30 × 5
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species   
#>           <dbl>       <dbl>        <dbl>       <dbl> <fct>     
#>  1          6.7         3.1          4.4         1.4 versicolor
#>  2          6.1         3            4.9         1.8 virginica 
#>  3          7           3.2          4.7         1.4 versicolor
#>  4          6           2.2          5           1.5 virginica 
#>  5          6.7         3.1          5.6         2.4 virginica 
#>  6          6.7         3.3          5.7         2.5 virginica 
#>  7          5.4         3.7          1.5         0.2 setosa    
#>  8          5.4         3.9          1.3         0.4 setosa    
#>  9          5           3.5          1.3         0.3 setosa    
#> 10          6.3         2.7          4.9         1.8 virginica 
#> # ℹ 20 more rows

## ------------------------------------------------------------------
## 7. Per-group: single column, explicit p
## ------------------------------------------------------------------
run_randomize(
  iris,
  per_group = list(list("Species", p = 0.6)),
  seed = 2024
)
#> # A tibble: 90 × 5
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species   
#>           <dbl>       <dbl>        <dbl>       <dbl> <fct>     
#>  1          5           3.2          1.2         0.2 setosa    
#>  2          7           3.2          4.7         1.4 versicolor
#>  3          7.1         3            5.9         2.1 virginica 
#>  4          6.2         3.4          5.4         2.3 virginica 
#>  5          6.5         3            5.2         2   virginica 
#>  6          5.6         2.5          3.9         1.1 versicolor
#>  7          5.1         3.8          1.5         0.3 setosa    
#>  8          6           3            4.8         1.8 virginica 
#>  9          6.5         2.8          4.6         1.5 versicolor
#> 10          7.3         2.9          6.3         1.8 virginica 
#> # ℹ 80 more rows

## ------------------------------------------------------------------
## 8. Per-group: single column, n/p omitted — falls back to top-level n
## ------------------------------------------------------------------
run_randomize(
  iris,
  n = 5,
  per_group = list(list("Species")),
  seed = 2024
)
#> Error: `per_group` element 1 is missing a valid column name.
#> Supply column name(s) as the first (unnamed) element:
#>   list("ColumnName", n = 3)
#>   list(c("Col1", "Col2"), n = 2)

## ------------------------------------------------------------------
## 9. Per-group: two columns, joint = TRUE (default)
##    Samples n = 2 from every unique cyl x gear combination.
##    Expected combinations and row counts in mtcars:
##      4:3 (1), 4:4 (8), 4:5 (2), 6:3 (2), 6:4 (4),
##      6:5 (1), 8:3 (12), 8:5 (2)
##    => returns 1+2+2+2+2+1+2+2 = 14 rows
##       (4:3 and 6:5 are undersized; all their rows are returned)
## ------------------------------------------------------------------
mtcars2 <- mtcars
mtcars2$cyl  <- factor(mtcars2$cyl)
mtcars2$gear <- factor(mtcars2$gear)

run_randomize(
  mtcars2,
  per_group = list(list(c("cyl", "gear"), n = 2)),
  joint = TRUE,
  seed = 99
)
#> Warning: Some joint groups have fewer rows than n = 2:
#>   - '4:3' (1 row(s))
#>   - '6:5' (1 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> # A tibble: 14 × 11
#>      mpg cyl    disp    hp  drat    wt  qsec    vs    am gear   carb
#>    <dbl> <fct> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <fct> <dbl>
#>  1  10.4 8     472     205  2.93  5.25  18.0     0     0 3         4
#>  2  15.8 8     351     264  4.22  3.17  14.5     0     1 5         4
#>  3  15   8     301     335  3.54  3.57  14.6     0     1 5         8
#>  4  17.8 6     168.    123  3.92  3.44  18.9     1     0 4         4
#>  5  21.4 6     258     110  3.08  3.22  19.4     1     0 3         1
#>  6  21   6     160     110  3.9   2.88  17.0     0     1 4         4
#>  7  30.4 4      75.7    52  4.93  1.62  18.5     1     1 4         2
#>  8  26   4     120.     91  4.43  2.14  16.7     0     1 5         2
#>  9  19.7 6     145     175  3.62  2.77  15.5     0     1 5         6
#> 10  19.2 8     400     175  3.08  3.84  17.0     0     0 3         2
#> 11  22.8 4     141.     95  3.92  3.15  22.9     1     0 4         2
#> 12  30.4 4      95.1   113  3.77  1.51  16.9     1     1 5         2
#> 13  18.1 6     225     105  2.76  3.46  20.2     1     0 3         1
#> 14  21.5 4     120.     97  3.7   2.46  20.0     1     0 3         1

## ------------------------------------------------------------------
## 10. Per-group: two columns, joint = FALSE, cascading (default)
##     Step 1 — group by cyl, take n = 2 per level => 6 rows
##     Step 2 — group by gear, take n = 3 per level from those 6 rows
##     Final row count depends on how many gear levels survived step 1
## ------------------------------------------------------------------
run_randomize(
  mtcars2,
  per_group = list(list(c("cyl", "gear"), n = c(2, 3))),
  joint = FALSE,
  joint_behaviour = "cascading",
  seed = 99
)
#> Warning: Some groups in column 'gear' have fewer rows than n = 3:
#>   - '3' (2 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> # A tibble: 5 × 11
#>     mpg cyl    disp    hp  drat    wt  qsec    vs    am gear   carb
#>   <dbl> <fct> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <fct> <dbl>
#> 1  17.8 6     168.    123  3.92  3.44  18.9     1     0 4         4
#> 2  33.9 4      71.1    65  4.22  1.84  19.9     1     1 4         1
#> 3  19.2 6     168.    123  3.92  3.44  18.3     1     0 4         4
#> 4  16.4 8     276.    180  3.07  4.07  17.4     0     0 3         3
#> 5  15.2 8     304     150  3.15  3.44  17.3     0     0 3         2

## ------------------------------------------------------------------
## 11. Per-group: two columns, joint = FALSE, independent
##     cyl sampled with n = 2 per level from original data  => row set A
##     gear sampled with n = 3 per level from original data => row set B
##     Returns union of A and B (rows qualifying under either column)
## ------------------------------------------------------------------
run_randomize(
  mtcars2,
  per_group = list(list(c("cyl", "gear"), n = c(2, 3))),
  joint = FALSE,
  joint_behaviour = "independent",
  seed = 99
)
#>                    mpg cyl  disp  hp drat    wt  qsec vs am gear carb
#> Merc 450SE        16.4   8 275.8 180 3.07 4.070 17.40  0  0    3    3
#> Merc 280          19.2   6 167.6 123 3.92 3.440 18.30  1  0    4    4
#> Lotus Europa      30.4   4  95.1 113 3.77 1.513 16.90  1  1    5    2
#> Ferrari Dino      19.7   6 145.0 175 3.62 2.770 15.50  0  1    5    6
#> Datsun 710        22.8   4 108.0  93 3.85 2.320 18.61  1  1    4    1
#> Merc 450SL        17.3   8 275.8 180 3.07 3.730 17.60  0  0    3    3
#> Hornet Sportabout 18.7   8 360.0 175 3.15 3.440 17.02  0  0    3    2
#> Camaro Z28        13.3   8 350.0 245 3.73 3.840 15.41  0  0    3    4
#> Volvo 142E        21.4   4 121.0 109 4.11 2.780 18.60  1  1    4    2
#> Toyota Corolla    33.9   4  71.1  65 4.22 1.835 19.90  1  1    4    1
#> Merc 280C         17.8   6 167.6 123 3.92 3.440 18.90  1  0    4    4
#> Porsche 914-2     26.0   4 120.3  91 4.43 2.140 16.70  0  1    5    2
#> Merc 240D         24.4   4 146.7  62 3.69 3.190 20.00  1  0    4    2
#> AMC Javelin       15.2   8 304.0 150 3.15 3.435 17.30  0  0    3    2
#> Fiat X1-9         27.3   4  79.0  66 4.08 1.935 18.90  1  1    4    1

## ------------------------------------------------------------------
## 12. Per-group: column name with special characters
## ------------------------------------------------------------------
df <- data.frame(
  `Age Group` = c("Young", "Young", "Middle", "Middle", "Old", "Old"),
  Score       = c(80, 85, 70, 75, 90, 95),
  check.names = FALSE
)

run_randomize(
  df,
  per_group = list(list("Age Group", n = 1)),
  seed = 1L
)
#> # A tibble: 3 × 2
#>   `Age Group` Score
#>   <chr>       <dbl>
#> 1 Middle         70
#> 2 Old            95
#> 3 Young          80

## ------------------------------------------------------------------
## 13. If both n and p are supplied, n takes priority (with warning)
## ------------------------------------------------------------------
run_randomize(iris, n = 10, p = 0.5, seed = 1L)
#> Warning: Both `n` and `p` were supplied. `n` takes priority; `p` is ignored.
#> To use `p`, omit the `n` argument.
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
#> 1           5.8         2.7          4.1         1.0 versicolor
#> 2           6.4         2.8          5.6         2.1  virginica
#> 3           4.4         3.2          1.3         0.2     setosa
#> 4           4.3         3.0          1.1         0.1     setosa
#> 5           7.0         3.2          4.7         1.4 versicolor
#> 6           5.4         3.0          4.5         1.5 versicolor
#> 7           5.4         3.4          1.7         0.2     setosa
#> 8           7.6         3.0          6.6         2.1  virginica
#> 9           6.1         2.8          4.7         1.2 versicolor
#> 10          4.6         3.4          1.4         0.3     setosa
```
