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
  per = c("all", "row", "col"),
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

- per:

  A character string specifying the randomization strategy for data
  frames, matrices, or tibbles.

  - **"all"** (default): Will randomize all of the items in the data
    frame. Any value can be placed anywhere in the randomized data
    frame.

  - **"row"**: Will randomize the items row-wise. The randomization is
    done for each row independently.

  - **"col"**: Will randomize the items column-wise. The randomization
    is done for each column independently.

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
#>                         mpg    cyl   disp      hp    drat      wt    qsec
#> Camaro Z28          440.000   8.00   1.00 275.800   8.000   4.000   3.000
#> Lotus Europa        180.000   6.00   1.00   3.000   3.000   1.000   1.000
#> Duster 360            0.000 146.70   3.62  13.300   5.424   4.000   4.000
#> Fiat 128             27.300  10.40   0.00   8.000   4.000   0.000   2.465
#> Chrysler Imperial     2.000  14.60   2.20   4.000   6.000  16.900  26.000
#> Toyota Corolla        4.000 145.00  16.87   3.000  17.800   2.140   1.000
#> Merc 230              8.000   8.00  15.84  21.500  62.000  19.200   2.760
#> Toyota Corona         0.000   0.00 335.00   0.000   1.000   4.220   0.000
#> Merc 280C            79.000  19.47   3.15 167.600 351.000  22.800   3.900
#> Merc 450SE          120.300   2.00 120.10   3.215   1.000   3.000 123.000
#> Merc 240D             1.000 275.80   4.00   1.000  33.900   2.000   1.000
#> Ford Pantera L        0.000  21.40   4.08  19.440   1.000   0.000   0.000
#> Hornet Sportabout   150.000   4.00 175.00 225.000   0.000 264.000 175.000
#> Pontiac Firebird      6.000   3.54   4.00   1.000   4.000   0.000   0.000
#> Valiant               4.000 360.00 180.00   3.000  18.900   6.000   1.000
#> Merc 280              3.080   3.07 108.00   8.000   6.000   1.513   1.000
#> Maserati Bora         3.000   5.00   0.00   3.190   4.070   1.835   8.000
#> Fiat X1-9             2.930   0.00 230.00   3.840   2.000   2.000   3.150
#> Volvo 142E            3.210   1.00   4.11  19.700   1.000   4.000  20.220
#> Ferrari Dino         15.200  17.30 105.00  21.000   0.000  17.300  20.000
#> Honda Civic          14.300  65.00 110.00  18.900  71.100   4.000   4.000
#> Datsun 710           18.610   3.46  19.20   3.440   0.000   3.070   8.000
#> Dodge Challenger     18.520   1.00  19.90   0.000   3.435  52.000   1.000
#> Merc 450SL            1.000   3.00   3.44   0.000   4.000   1.000   1.000
#> Hornet 4 Drive        4.000  17.05   4.00  18.700 350.000  15.500  95.000
#> Lincoln Continental   3.920   2.00 167.60   5.000   5.345   3.000   2.320
#> Merc 450SLC          17.020   0.00   3.69   0.000 180.000   8.000   2.000
#> Cadillac Fleetwood    1.000   4.00   4.00  75.700  15.500   3.730  15.800
#> AMC Javelin         275.800 318.00   3.78   8.000   1.000   0.000 160.000
#> Mazda RX4 Wag        15.200  18.60 400.00   4.430   0.000   3.730   0.000
#> Porsche 914-2        21.400  17.42 123.00 205.000  17.020  17.600  97.000
#> Mazda RX4             1.615   4.00   0.00   4.000   3.000   6.000   3.000
#>                         vs     am    gear    carb
#> Camaro Z28            3.00   4.93   3.000   4.000
#> Lotus Europa          3.00   0.00  17.400  20.010
#> Duster 360            3.15   4.00  66.000   0.000
#> Fiat 128              8.00 140.80   4.080   3.000
#> Chrysler Imperial    21.00  15.00   4.000   1.000
#> Toyota Corolla        4.00   3.00 175.000 160.000
#> Merc 230             66.00   3.70  93.000   1.000
#> Toyota Corona         1.00   8.00   0.000   3.850
#> Merc 280C           360.00   3.92   8.000   3.080
#> Merc 450SE           18.00   1.00  17.980   4.220
#> Merc 240D            15.41 304.00 121.000 301.000
#> Ford Pantera L        2.00  16.46  30.400 245.000
#> Hornet Sportabout    95.10   5.25  32.400   3.170
#> Pontiac Firebird      2.62 245.00   4.000  18.300
#> Valiant               4.00 113.00   0.000   0.000
#> Merc 280             91.00 109.00   1.935  30.400
#> Maserati Bora         2.00   4.00 258.000  22.900
#> Fiat X1-9             8.00   3.57   4.000   0.000
#> Volvo 142E            1.00   8.00   5.000   0.000
#> Ferrari Dino         17.82   0.00   6.000   4.000
#> Honda Civic         215.00   3.44   1.000  18.100
#> Datsun 710            3.00 460.00   1.000   1.000
#> Dodge Challenger      1.00   2.77  24.400   3.570
#> Merc 450SL           78.70  14.50   3.230   4.000
#> Hornet 4 Drive        3.90   2.78   3.920   3.770
#> Lincoln Continental 110.00   1.00   2.760  10.400
#> Merc 450SLC           3.07   3.52 110.000   2.875
#> Cadillac Fleetwood    0.00   0.00   3.000   5.000
#> AMC Javelin           0.00   4.00 472.000   2.000
#> Mazda RX4 Wag        14.70   0.00  22.800   4.000
#> Porsche 914-2       150.00  16.40   5.000   3.000
#> Mazda RX4             1.00   6.00   3.845  16.700

## ------------------------------------------------------------------
## 2. Sample a fixed number of rows
## ------------------------------------------------------------------
run_randomize(mtcars, n = 5)
#>                     mpg   cyl disp    hp  drat   wt   qsec    vs     am  gear
#> Camaro Z28        440.0   8.0 1.00 275.8 8.000  4.0  3.000  3.00   4.93  3.00
#> Lotus Europa      180.0   6.0 1.00   3.0 3.000  1.0  1.000  3.00   0.00 17.40
#> Duster 360          0.0 146.7 3.62  13.3 5.424  4.0  4.000  3.15   4.00 66.00
#> Fiat 128           27.3  10.4 0.00   8.0 4.000  0.0  2.465  8.00 140.80  4.08
#> Chrysler Imperial   2.0  14.6 2.20   4.0 6.000 16.9 26.000 21.00  15.00  4.00
#>                    carb
#> Camaro Z28         4.00
#> Lotus Europa      20.01
#> Duster 360         0.00
#> Fiat 128           3.00
#> Chrysler Imperial  1.00

## ------------------------------------------------------------------
## 3. Sample a proportion of rows
## ------------------------------------------------------------------
run_randomize(mtcars, p = 0.25)
#>                     mpg   cyl   disp    hp   drat    wt   qsec    vs     am
#> Camaro Z28        440.0   8.0   1.00 275.8  8.000  4.00  3.000  3.00   4.93
#> Lotus Europa      180.0   6.0   1.00   3.0  3.000  1.00  1.000  3.00   0.00
#> Duster 360          0.0 146.7   3.62  13.3  5.424  4.00  4.000  3.15   4.00
#> Fiat 128           27.3  10.4   0.00   8.0  4.000  0.00  2.465  8.00 140.80
#> Chrysler Imperial   2.0  14.6   2.20   4.0  6.000 16.90 26.000 21.00  15.00
#> Toyota Corolla      4.0 145.0  16.87   3.0 17.800  2.14  1.000  4.00   3.00
#> Merc 230            8.0   8.0  15.84  21.5 62.000 19.20  2.760 66.00   3.70
#> Toyota Corona       0.0   0.0 335.00   0.0  1.000  4.22  0.000  1.00   8.00
#>                     gear   carb
#> Camaro Z28          3.00   4.00
#> Lotus Europa       17.40  20.01
#> Duster 360         66.00   0.00
#> Fiat 128            4.08   3.00
#> Chrysler Imperial   4.00   1.00
#> Toyota Corolla    175.00 160.00
#> Merc 230           93.00   1.00
#> Toyota Corona       0.00   3.85

## ------------------------------------------------------------------
## 4. Sampling with replacement
## ------------------------------------------------------------------
run_randomize(mtcars, n = 10, replace = TRUE)
#>                     mpg    cyl   disp    hp    drat    wt   qsec     vs     am
#> Camaro Z28        440.0   8.00   1.00 275.8   8.000  4.00  3.000   3.00   4.93
#> Lotus Europa      180.0   6.00   1.00   3.0   3.000  1.00  1.000   3.00   0.00
#> Duster 360          0.0 146.70   3.62  13.3   5.424  4.00  4.000   3.15   4.00
#> Fiat 128           27.3  10.40   0.00   8.0   4.000  0.00  2.465   8.00 140.80
#> Chrysler Imperial   2.0  14.60   2.20   4.0   6.000 16.90 26.000  21.00  15.00
#> Toyota Corolla      4.0 145.00  16.87   3.0  17.800  2.14  1.000   4.00   3.00
#> Porsche 914-2      21.4  17.42 123.00 205.0  17.020 17.60 97.000 150.00  16.40
#> Merc 230            8.0   8.00  15.84  21.5  62.000 19.20  2.760  66.00   3.70
#> Toyota Corona       0.0   0.00 335.00   0.0   1.000  4.22  0.000   1.00   8.00
#> Merc 280C          79.0  19.47   3.15 167.6 351.000 22.80  3.900 360.00   3.92
#>                     gear   carb
#> Camaro Z28          3.00   4.00
#> Lotus Europa       17.40  20.01
#> Duster 360         66.00   0.00
#> Fiat 128            4.08   3.00
#> Chrysler Imperial   4.00   1.00
#> Toyota Corolla    175.00 160.00
#> Porsche 914-2       5.00   3.00
#> Merc 230           93.00   1.00
#> Toyota Corona       0.00   3.85
#> Merc 280C           8.00   3.08

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
#> Warning: Some groups in column 'Species' have fewer rows than n = 10:
#>   - '0.1' (2 row(s))
#>   - '0.2' (5 row(s))
#>   - '0.3' (1 row(s))
#>   - '0.4' (1 row(s))
#>   - '0.5' (1 row(s))
#>   - '0.6' (1 row(s))
#>   - '1.0' (3 row(s))
#>   - '1.3' (6 row(s))
#>   - '1.4' (3 row(s))
#>   - '1.5' (7 row(s))
#>   - '1.8' (2 row(s))
#>   - '1.9' (2 row(s))
#>   - '2.0' (3 row(s))
#>   - '2.1' (2 row(s))
#>   - '2.2' (1 row(s))
#>   - '2.3' (3 row(s))
#>   - '2.4' (1 row(s))
#>   - '2.5' (2 row(s))
#>   - '2.6' (2 row(s))
#>   - '2.7' (1 row(s))
#>   - '2.8' (1 row(s))
#>   - '2.9' (2 row(s))
#>   - '3.0' (4 row(s))
#>   - '3.1' (1 row(s))
#>   - '3.2' (1 row(s))
#>   - '3.3' (1 row(s))
#>   - '3.4' (3 row(s))
#>   - '3.5' (1 row(s))
#>   - '3.6' (1 row(s))
#>   - '3.7' (2 row(s))
#>   - '3.8' (2 row(s))
#>   - '3.9' (2 row(s))
#>   - '4.0' (1 row(s))
#>   - '4.1' (1 row(s))
#>   - '4.3' (2 row(s))
#>   - '4.5' (4 row(s))
#>   - '4.6' (2 row(s))
#>   - '4.7' (3 row(s))
#>   - '4.8' (1 row(s))
#>   - '4.9' (4 row(s))
#>   - '5.0' (4 row(s))
#>   - '5.1' (2 row(s))
#>   - '5.4' (2 row(s))
#>   - '5.5' (1 row(s))
#>   - '5.6' (3 row(s))
#>   - '5.7' (2 row(s))
#>   - '5.8' (3 row(s))
#>   - '5.9' (1 row(s))
#>   - '6.1' (3 row(s))
#>   - '6.2' (1 row(s))
#>   - '6.3' (1 row(s))
#>   - '6.5' (2 row(s))
#>   - '6.6' (2 row(s))
#>   - '6.7' (3 row(s))
#>   - '6.9' (1 row(s))
#>   - '7.2' (1 row(s))
#>   - '7.9' (1 row(s))
#>   - 'setosa' (8 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> # A tibble: 149 × 5
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species   
#>    <chr>        <chr>       <chr>        <chr>       <chr>     
#>  1 5.1          6.1         1.8          setosa      1.5       
#>  2 3.1          2.3         setosa       5.5         versicolor
#>  3 5.7          setosa      5.7          6.4         1.5       
#>  4 4.4          0.2         6.1          2.5         3.1       
#>  5 5.4          0.2         3.7          5.1         1.5       
#>  6 setosa       2.9         5.2          1.6         versicolor
#>  7 2.3          versicolor  1.6          0.2         4.5       
#>  8 1.8          5.0         6.4          versicolor  4.9       
#>  9 versicolor   6.7         6.3          4.8         4.9       
#> 10 virginica    5.4         7.7          virginica   virginica 
#> # ℹ 139 more rows

## ------------------------------------------------------------------
## 7. Per-group: single column, explicit p
## ------------------------------------------------------------------
run_randomize(
  iris,
  per_group = list(list("Species", p = 0.6)),
  seed = 2024
)
#> # A tibble: 61 × 5
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width Species   
#>    <chr>        <chr>       <chr>        <chr>       <chr>     
#>  1 7.4          3.0         1.5          3.0         2.1       
#>  2 4.9          4.7         virginica    2.8         6.6       
#>  3 setosa       versicolor  virginica    5.0         0.2       
#>  4 1.4          6.8         1.0          setosa      virginica 
#>  5 3.0          0.2         versicolor   6.2         1.5       
#>  6 1.9          2.7         7.2          1.8         5.6       
#>  7 4.0          5.9         2.5          1.4         versicolor
#>  8 7.7          5.6         1.6          virginica   5.8       
#>  9 versicolor   virginica   3.0          3.6         1.9       
#> 10 versicolor   6.1         5.1          setosa      virginica 
#> # ℹ 51 more rows

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
#>   - '4: 66' (1 row(s))
#>   - '1:0' (1 row(s))
#>   - '10.4:0' (1 row(s))
#>   - '360.0:0' (1 row(s))
#>   - ' 62:1' (1 row(s))
#>   - '301.0:1' (1 row(s))
#>   - '3:110' (1 row(s))
#>   - '17.98:13.3' (1 row(s))
#>   - '22.90:14.7' (1 row(s))
#>   - '4:16.4' (1 row(s))
#>   - '4:17.02' (1 row(s))
#>   - '3:18.1' (1 row(s))
#>   - '0:2' (1 row(s))
#>   - '1:2' (1 row(s))
#>   - '123:2' (1 row(s))
#>   - '0:2.140' (1 row(s))
#>   - ' 93:21.4' (1 row(s))
#>   - '8:21.4' (1 row(s))
#>   - '0:21.5' (1 row(s))
#>   - '4.11:215' (1 row(s))
#>   - '120.3:27.3' (1 row(s))
#>   - '245:275.8' (1 row(s))
#>   - '1:3' (1 row(s))
#>   - '1.513:3.23' (1 row(s))
#>   - '0:3.440' (1 row(s))
#>   - '1:3.73' (1 row(s))
#>   - '4:3.730' (1 row(s))
#>   - ' 91:3.90' (1 row(s))
#>   - '2.465:30.4' (1 row(s))
#>   - '110:4' (1 row(s))
#>   - '5:4' (1 row(s))
#>   - '3.845:6' (1 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> # A tibble: 32 × 11
#>    mpg   cyl     disp    hp    drat    wt    qsec  vs    am      gear  carb 
#>    <chr> <chr>   <chr>   <chr> <chr>   <chr> <chr> <chr> <chr>   <chr> <chr>
#>  1 0     "3"     "351.0" 2.93  "5.250" 4     17.60 0     "1"     18.1  3.00 
#>  2 3     "17.98" "0"     26.0  " 71.1" 14.60 2.770 4     "14.50" 13.3  3.92 
#>  3 0     " 91"   "19.7"  2     "3.190" 15.0  16.46 2     " 95.1" 3.90  0    
#>  4 8     "1"     "1"     18.60 "3"     275.8 15.50 0     "0"     2     1    
#>  5 109   "1"     "22.8"  167.6 "4"     110   1     2.200 "16.70" 0     3.440
#>  6 4     " 93"   " 79.0" 0     "0"     1     150   3.780 "19.2"  21.4  175  
#>  7 14.3  "0"     "5"     5.424 "1"     121.0 3.07  3.440 "3.150" 2.140 17.30
#>  8 3.08  "245"   "4"     0     "1"     1     24.4  4     "4"     275.8 8    
#>  9 3.570 "8"     "8"     0     "0"     1     146.7 140.8 "1"     21.4  3.90 
#> 10 3     "2.465" "1"     3.70  "2.620" 3.840 3.170 17.82 "1"     30.4  19.90
#> # ℹ 22 more rows

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
#> Warning: Some groups in column 'cyl' have fewer rows than n = 2:
#>   - ' 62' (1 row(s))
#>   - ' 91' (1 row(s))
#>   - ' 93' (1 row(s))
#>   - '1.513' (1 row(s))
#>   - '10.4' (1 row(s))
#>   - '110' (1 row(s))
#>   - '120.3' (1 row(s))
#>   - '123' (1 row(s))
#>   - '17.98' (1 row(s))
#>   - '2.465' (1 row(s))
#>   - '22.90' (1 row(s))
#>   - '245' (1 row(s))
#>   - '3.845' (1 row(s))
#>   - '301.0' (1 row(s))
#>   - '360.0' (1 row(s))
#>   - '4.11' (1 row(s))
#>   - '5' (1 row(s))
#>   - '8' (1 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> Warning: Some groups in column 'gear' have fewer rows than n = 3:
#>   - '0' (2 row(s))
#>   - '1' (2 row(s))
#>   - '110' (1 row(s))
#>   - '13.3' (1 row(s))
#>   - '14.7' (1 row(s))
#>   - '16.4' (1 row(s))
#>   - '18.1' (1 row(s))
#>   - '2' (2 row(s))
#>   - '21.4' (2 row(s))
#>   - '215' (1 row(s))
#>   - '27.3' (1 row(s))
#>   - '275.8' (1 row(s))
#>   - '3' (1 row(s))
#>   - '3.23' (1 row(s))
#>   - '3.440' (1 row(s))
#>   - '3.73' (1 row(s))
#>   - '3.730' (1 row(s))
#>   - '3.90' (1 row(s))
#>   - '30.4' (1 row(s))
#>   - '4' (2 row(s))
#>   - '6' (1 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> # A tibble: 26 × 11
#>    mpg   cyl     disp    hp    drat  wt    qsec  vs      am      gear  carb   
#>    <chr> <chr>   <chr>   <chr> <chr> <chr> <chr> <chr>   <chr>   <chr> <chr>  
#>  1 0     "3"     "351.0" 2.93  5.250 4     17.60 "0"     "1"     18.1  "3.00" 
#>  2 4     "0"     "1.615" 1     0     0     1     "4"     "160.0" 3.440 " 65"  
#>  3 2     "5"     "18.30" 1.835 400.0 17.40 1     "3.69"  "230"   4     "4"    
#>  4 0     "120.3" "17.42" 108.0 180   4     4     " 75.7" "180"   27.3  "3.07" 
#>  5 318.0 "4"     "4.93"  4     18.90 17.3  4     "3.15"  "180"   3.730 "20.01"
#>  6 0     " 91"   "19.7"  2     3.190 15.0  16.46 "2"     " 95.1" 3.90  "0"    
#>  7 4     "110"   "4"     3.92  16.90 6     0     "4"     "2.875" 4     "1"    
#>  8 350.0 "1"     "20.00" 6     8     3     1     "18.90" "4"     3.73  "15.2" 
#>  9 4     " 93"   " 79.0" 0     0     1     150   "3.780" "19.2"  21.4  "175"  
#> 10 145.0 "1.513" "440.0" 1     0     4     3     "17.02" "4.08"  3.23  "3.85" 
#> # ℹ 16 more rows

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
#> Warning: Some groups in column 'cyl' have fewer rows than n = 2:
#>   - ' 62' (1 row(s))
#>   - ' 91' (1 row(s))
#>   - ' 93' (1 row(s))
#>   - '1.513' (1 row(s))
#>   - '10.4' (1 row(s))
#>   - '110' (1 row(s))
#>   - '120.3' (1 row(s))
#>   - '123' (1 row(s))
#>   - '17.98' (1 row(s))
#>   - '2.465' (1 row(s))
#>   - '22.90' (1 row(s))
#>   - '245' (1 row(s))
#>   - '3.845' (1 row(s))
#>   - '301.0' (1 row(s))
#>   - '360.0' (1 row(s))
#>   - '4.11' (1 row(s))
#>   - '5' (1 row(s))
#>   - '8' (1 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#> Warning: Some groups in column 'gear' have fewer rows than n = 3:
#>   - ' 66' (1 row(s))
#>   - '1' (2 row(s))
#>   - '110' (1 row(s))
#>   - '13.3' (1 row(s))
#>   - '14.7' (1 row(s))
#>   - '16.4' (1 row(s))
#>   - '17.02' (1 row(s))
#>   - '18.1' (1 row(s))
#>   - '2.140' (1 row(s))
#>   - '21.4' (2 row(s))
#>   - '21.5' (1 row(s))
#>   - '215' (1 row(s))
#>   - '27.3' (1 row(s))
#>   - '275.8' (1 row(s))
#>   - '3' (1 row(s))
#>   - '3.23' (1 row(s))
#>   - '3.440' (1 row(s))
#>   - '3.73' (1 row(s))
#>   - '3.730' (1 row(s))
#>   - '3.90' (1 row(s))
#>   - '30.4' (1 row(s))
#>   - '4' (2 row(s))
#>   - '6' (1 row(s))
#> All rows from these groups will be returned.
#> Set `replace = TRUE` to allow oversampling.
#>                       mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear
#> Honda Civic         15.84 301.0   264     5 18.00  15.2 120.1  30.4 18.61     1
#> Merc 450SE          350.0     1 20.00     6     8     3     1 18.90     4  3.73
#> Toyota Corolla      145.0 1.513 440.0     1     0     4     3 17.02  4.08  3.23
#> Maserati Bora       3.570     8     8     0     0     1 146.7 140.8     1  21.4
#> Toyota Corona           0 120.3 17.42 108.0   180     4     4  75.7   180  27.3
#> Pontiac Firebird        8     1     6     6     3     0     3     3   123     3
#> AMC Javelin             4   123 2.780   335  15.8     3     3  18.7     0     2
#> Datsun 710           15.5  4.11   175 4.070     4     8 3.520     1   150   215
#> Ferrari Dino            4    93  79.0     0     0     1   150 3.780  19.2  21.4
#> Fiat X1-9               8     1     1 18.60     3 275.8 15.50     0     0     2
#> Chrysler Imperial    14.3     0     5 5.424     1 121.0  3.07 3.440 3.150 2.140
#> Lotus Europa            4   110     4  3.92 16.90     6     0     4 2.875     4
#> Valiant               109     1  22.8 167.6     4   110     1 2.200 16.70     0
#> Mazda RX4 Wag           3 17.98     0  26.0  71.1 14.60 2.770     4 14.50  13.3
#> Mazda RX4               4     0 1.615     1     0     0     1     4 160.0 3.440
#> Merc 240D            33.9     0  17.8 360.0    95     4 16.87 15.41     4     2
#> Merc 230            258.0 22.90     0     4     2     1    52     4 460.0  14.7
#> Porsche 914-2       1.935  10.4    66 167.6     0     1     0     8     2     0
#> Hornet Sportabout       0    91  19.7     2 3.190  15.0 16.46     2  95.1  3.90
#> Merc 280                0    62     1 17.05  32.4 19.44    97     2 225.0     1
#> Merc 280C           19.47     4     6  22.8  4.22   245  3.77  21.0     6    66
#> Hornet 4 Drive          0 360.0     3     2  3.15  78.7     5 3.435  3.21     0
#> Ford Pantera L       3.07     0  4.43     3 472.0 160.0  21.0     1     8  21.5
#> Merc 450SLC         318.0     4  4.93     4 18.90  17.3     4  3.15   180 3.730
#> Cadillac Fleetwood   3.08   245     4     0     1     1  24.4     4     4 275.8
#> Lincoln Continental 275.8     4  2.76 3.215  4.22 2.320  3.54     8  19.2 17.02
#> Dodge Challenger        1     4   105     0     0 304.0     3     0     8  16.4
#> Fiat 128                2     5 18.30 1.835 400.0 17.40     1  3.69   230     4
#> Duster 360              0     3 351.0  2.93 5.250     4 17.60     0     1  18.1
#> Camaro Z28              3 2.465     1  3.70 2.620 3.840 3.170 17.82     1  30.4
#> Merc 450SL           3.62     3     3     1     8     1     0 20.22     8   110
#> Volvo 142E          3.460 3.845   205 3.570     8  3.92     4 18.52     3     6
#>                      carb
#> Honda Civic             1
#> Merc 450SE           15.2
#> Toyota Corolla       3.85
#> Maserati Bora        3.90
#> Toyota Corona        3.07
#> Pontiac Firebird    5.345
#> AMC Javelin             8
#> Datsun 710           4.08
#> Ferrari Dino          175
#> Fiat X1-9               1
#> Chrysler Imperial   17.30
#> Lotus Europa            1
#> Valiant             3.440
#> Mazda RX4 Wag        3.92
#> Mazda RX4              65
#> Merc 240D            3.08
#> Merc 230              113
#> Porsche 914-2        2.76
#> Hornet Sportabout       0
#> Merc 280                0
#> Merc 280C             175
#> Hornet 4 Drive       10.4
#> Ford Pantera L          5
#> Merc 450SLC         20.01
#> Cadillac Fleetwood      8
#> Lincoln Continental     6
#> Dodge Challenger        1
#> Fiat 128                4
#> Duster 360           3.00
#> Camaro Z28          19.90
#> Merc 450SL              0
#> Volvo 142E              4

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
#> # A tibble: 5 × 2
#>   `Age Group` Score 
#>   <chr>       <chr> 
#> 1 70          Middle
#> 2 Middle      85    
#> 3 Young       95    
#> 4 Old         75    
#> 5 80          Old   

## ------------------------------------------------------------------
## 13. If both n and p are supplied, n takes priority (with warning)
## ------------------------------------------------------------------
run_randomize(iris, n = 10, p = 0.5, seed = 1L)
#> Warning: Both `n` and `p` were supplied. `n` takes priority; `p` is ignored.
#> To use `p`, omit the `n` argument.
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
#> 1           3.4      setosa          2.0         3.4        5.0
#> 2           3.8  versicolor          6.1         2.4        5.2
#> 3           5.7         4.8          2.7         0.2        1.4
#> 4           3.5         5.7          1.8         4.7     setosa
#> 5    versicolor         3.5   versicolor  versicolor versicolor
#> 6           1.4         1.0          3.6  versicolor versicolor
#> 7           0.4  versicolor       setosa         1.5        3.7
#> 8        setosa         1.9          4.9         0.1        1.3
#> 9    versicolor         3.0          2.9  versicolor        2.5
#> 10          6.7         5.5          1.5         7.7        0.2
```
