# Perform Correlation Analysis

Computes pairwise or targeted correlation analysis between numeric,
ordinal, or binary variables in a data frame, tibble, or matrix.
Supports automatic method selection (Pearson, Spearman, Kendall,
point-biserial, phi, Cramer's V, polychoric, and polyserial) based on
data characteristics. Returns tidy output tables of correlation
estimates, masked correlations, confidence intervals, and p-values.

## Usage

``` r
run_correl(
  x,
  y = NULL,
  z = NULL,
  metadata = NULL,
  remove = NULL,
  method = "auto",
  force_test = NULL,
  no_kendall = FALSE,
  remove_y = TRUE,
  sig_threshold = 0.05,
  normality_method = "auto",
  alpha_normality = 0.05,
  interpretation = NULL,
  summary_table_only = FALSE,
  method_in_row = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A data frame, tibble, or matrix whose columns are the variables to
  correlate. Column names may contain special characters.

- y:

  A character vector naming one or more columns in `x` to use as primary
  variables. When both `y` and `z` are `NULL`, all pairwise correlations
  among numeric columns of `x` are computed. When only `y` is supplied,
  correlations among variables in `y` and between `y` and all other
  columns are returned.

- z:

  A single character string naming a second column in `x`. When supplied
  together with `y`, correlations between each element of `y` and `z`
  are computed.

- metadata:

  A data frame with the same number of rows as `x` containing
  sample-level metadata. Appended to the returned tidy tables when
  provided but not used in correlation computation.

- remove:

  A named list of one-sided formulas of the form
  `list("colname" ~ c("level1", "level2"))` specifying levels of a
  column in `x` to exclude before analysis. The left-hand side must be a
  single column name; the right-hand side may be a character vector of
  one or more levels.

- method:

  A single character string specifying the correlation method. One of
  `"auto"` (default), `"pearson"`, `"spearman"`, `"kendall"`,
  `"pointbiserial"`, `"phi"`, `"cramerv"`, or `"poly"`. When `"auto"`,
  the method is selected per variable pair based on data type,
  normality, and sample size.

- force_test:

  A list of two-element lists, each specifying a variable pair and the
  correlation method to use for that pair, overriding the global
  `method` argument. Each element must be of the form
  `list(c("var1", "var2"), "method")`. For example:

      force_test = list(
        list(c("var1", "var2"), "pearson"),
        list(c("var1", "var3"), "spearman")
      )

  Valid methods are all values accepted by `method` except `"auto"`.
  Column names are order-insensitive within each pair. If a pair appears
  more than once, the last entry takes precedence with a warning. Pairs
  not listed here are analyzed using the global `method`. Default `NULL`
  (no overrides).

- no_kendall:

  Logical. If `TRUE`, the automatic method selection avoids Kendall's
  tau and uses Spearman's rho instead when assumptions for Pearson's are
  violated and sample size is small or ties are frequent. Note that
  specific tests provided in `force_test` take precedence. Default
  `FALSE`.

- remove_y:

  Logical. If `TRUE`, correlations between variables within the `y`
  vector are excluded, focusing the analysis only on `y` vs. other
  variables. Default `TRUE`.

- sig_threshold:

  A single numeric value in `(0, 1)` specifying the significance
  threshold for p-values. Correlations with `p >= sig_threshold` are
  masked in the masked table. Default `0.05`.

- normality_method:

  A single character string, one of `"auto"` (default), `"shapiro"`, or
  `"lilliefors"`, controlling the normality test used in automatic
  method selection. `"auto"` applies Shapiro-Wilk for \\n \< 50\\ and
  Lilliefors otherwise.

- alpha_normality:

  A single numeric value in `(0, 1)` specifying the significance
  threshold used internally for normality and variance tests when
  `method = "auto"`. Default `0.05`.

- interpretation:

  Character string or `NULL`. If not `NULL` and
  `summary_table_only = TRUE`, adds a qualitative interpretation column
  based on `"psychology"`, `"politics"`, or `"medicine"` thresholds.
  Default `NULL`.

- summary_table_only:

  Logical. If `TRUE`, returns only a single data frame with all result
  columns. Default is `FALSE`.

- method_in_row:

  Logical. When summary_table_only = TRUE and method_in_row = FALSE, the
  method column is removed. Instead, correlation estimates are converted
  to strings with alphabetical superscripts. Default `TRUE`.

- verbose:

  Logical. When `TRUE`, diagnostic messages about assumption checks and
  method selection are printed. Default `FALSE`.

- ...:

  Additional arguments passed to
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html).

## Value

A named list of class `"run_correl"` with the following elements:

- `correlations`:

  A tidy data frame of pairwise correlation estimates with columns
  `var1`, `var2`, `estimate`, `method_used`, and `n`.

- `correlations_masked`:

  Same as `correlations` but `estimate` is `NA` for pairs where
  `p_value >= sig_threshold`.

- `p_values`:

  A tidy data frame with columns `var1`, `var2`, and `p_value`.

- `confidence_intervals`:

  A tidy data frame with columns `var1`, `var2`, `ci_lower`, and
  `ci_upper`.

- `sig_threshold`:

  The significance threshold used.

- `method`:

  The method argument as supplied.

- `call`:

  The matched call.

## Details

**Automatic method selection** (`method = "auto"`):

|                |                                                                              |
|----------------|------------------------------------------------------------------------------|
| **Method**     | **Conditions**                                                               |
| Pearson        | Both continuous, normal, homoscedastic                                       |
| Spearman       | Ordinal or continuous, non-normal or heteroscedastic, \\n \ge 30\\, few ties |
| Kendall        | Ordinal or continuous, non-normal, \\n \< 30\\ or many ties                  |
| Point-biserial | One continuous, one binary variable                                          |
| Phi            | Both binary                                                                  |
| Cramer's V     | Both nominal (multi-category) factors                                        |
| Polychoric     | Both ordered factors (requires polycor)                                      |
| Polyserial     | One ordered factor, one continuous (requires polycor)                        |

Normality is assessed via Shapiro-Wilk (\\n \< 50\\) or Lilliefors (\\n
\ge 50\\; requires nortest). Homoscedasticity is assessed via Levene's
test (requires car) with Bartlett's test as fallback.

Confidence intervals for Cramer's V and polychoric/polyserial
correlations are not available and are returned as `NA`. When
`summary_table_only = TRUE`, the function returns a flat data frame
instead of the structured list.

**Interpretation thresholds** (`interpretation`):

*Psychology*:

|                          |                    |
|--------------------------|--------------------|
| **Absolute Correlation** | **Interpretation** |
| 1.0                      | Perfect            |
| 0.7 to \\\<\\ 1.0        | Strong             |
| 0.4 to \\\<\\ 0.7        | Moderate           |
| 0.1 to \\\<\\ 0.4        | Weak               |
| \\\<\\ 0.1               | Zero               |

*Politics*:

|                          |                    |
|--------------------------|--------------------|
| **Absolute Correlation** | **Interpretation** |
| 1.0                      | Perfect            |
| 0.7 to \\\<\\ 1.0        | Very Strong        |
| 0.4 to \\\<\\ 0.7        | Strong             |
| 0.3 to \\\<\\ 0.4        | Moderate           |
| 0.2 to \\\<\\ 0.3        | Weak               |
| 0.1 to \\\<\\ 0.2        | Negligible         |
| \\\<\\ 0.1               | None               |

*Medicine*:

|                          |                    |
|--------------------------|--------------------|
| **Absolute Correlation** | **Interpretation** |
| 1.0                      | Perfect            |
| 0.8 to \\\<\\ 1.0        | Very Strong        |
| 0.6 to \\\<\\ 0.8        | Moderate           |
| 0.3 to \\\<\\ 0.6        | Fair               |
| 0.1 to \\\<\\ 0.3        | Poor               |
| \\\<\\ 0.1               | None               |

*Categorical (Phi and Cramer's V)*:

|             |                    |
|-------------|--------------------|
| **Value**   | **Interpretation** |
| \\\>\\ 0.25 | Very strong        |
| \\\>\\ 0.15 | Strong             |
| \\\>\\ 0.10 | Moderate           |
| \\\>\\ 0.05 | Weak               |
| 0 to 0.05   | No or very weak    |

## References

Akoglu H. User's guide to correlation coefficients. Turk J Emerg Med.
2018 Aug 7;18(3):91-93. doi: 10.1016/j.tjem.2018.08.001. PMID: 30191186;
PMCID: PMC6107969.

## Author

John Lennon L. Calorio

## Examples

``` r
## Basic pairwise correlations (all numeric columns)
run_correl(mtcars)
#> [1] correlations         correlations_masked  p_values            
#> [4] confidence_intervals sig_threshold        method              
#> [7] call                
#> <0 rows> (or 0-length row.names)

## Correlate one variable against all others
run_correl(mtcars, y = "mpg")
#> [1] correlations         correlations_masked  p_values            
#> [4] confidence_intervals sig_threshold        method              
#> [7] call                
#> <0 rows> (or 0-length row.names)

## Single pair
run_correl(mtcars, y = "mpg", z = "wt")
#> [1] correlations         correlations_masked  p_values            
#> [4] confidence_intervals sig_threshold        method              
#> [7] call                
#> <0 rows> (or 0-length row.names)

## Specify method explicitly
run_correl(mtcars, method = "spearman")
#> [1] correlations         correlations_masked  p_values            
#> [4] confidence_intervals sig_threshold        method              
#> [7] call                
#> <0 rows> (or 0-length row.names)

## Exclude specific levels before analysis
df <- mtcars
df$cyl <- as.numeric(factor(df$cyl))
run_correl(df, remove = list("cyl" ~ c("4", "6")))
#> [1] correlations         correlations_masked  p_values            
#> [4] confidence_intervals sig_threshold        method              
#> [7] call                
#> <0 rows> (or 0-length row.names)

## Access individual output tables
res <- run_correl(mtcars)
res$correlations
#>    var1 var2    estimate   method_used  n
#> 1   mpg  cyl -0.79531341       kendall 32
#> 2   mpg disp -0.76813115       kendall 32
#> 3   mpg   hp -0.74281251       kendall 32
#> 4   mpg drat  0.65145546      spearman 32
#> 5   mpg   wt -0.88642203      spearman 32
#> 6   mpg qsec  0.46693575      spearman 32
#> 7   mpg   vs  0.66403892 pointbiserial 32
#> 8   mpg   am  0.59983243 pointbiserial 32
#> 9   mpg gear  0.43315089       kendall 32
#> 10  mpg carb -0.50439455       kendall 32
#> 11  cyl disp  0.81442625       kendall 32
#> 12  cyl   hp  0.78518650       kendall 32
#> 13  cyl drat -0.55131785       kendall 32
#> 14  cyl   wt  0.72826111       kendall 32
#> 15  cyl qsec -0.44896982       kendall 32
#> 16  cyl   vs -0.81081180 pointbiserial 32
#> 17  cyl   am -0.52260705 pointbiserial 32
#> 18  cyl gear -0.51254349       kendall 32
#> 19  cyl carb  0.46542994       kendall 32
#> 20 disp   hp  0.66599874       kendall 32
#> 21 disp drat -0.49898277       kendall 32
#> 22 disp   wt  0.74338240       kendall 32
#> 23 disp qsec -0.30081549       kendall 32
#> 24 disp   vs -0.71041589 pointbiserial 32
#> 25 disp   am -0.59122704 pointbiserial 32
#> 26 disp gear -0.47597955       kendall 32
#> 27 disp carb  0.41373600       kendall 32
#> 28   hp drat -0.38262689       kendall 32
#> 29   hp   wt  0.61130810       kendall 32
#> 30   hp qsec -0.47290613       kendall 32
#> 31   hp   vs -0.72309674 pointbiserial 32
#> 32   hp   am -0.24320426 pointbiserial 32
#> 33   hp gear -0.27944584       kendall 32
#> 34   hp carb  0.59598416       kendall 32
#> 35 drat   wt -0.75039041      spearman 32
#> 36 drat qsec  0.09186863      spearman 32
#> 37 drat   vs  0.44027846 pointbiserial 32
#> 38 drat   am  0.71271113 pointbiserial 32
#> 39 drat gear  0.58392476       kendall 32
#> 40 drat carb -0.09535193       kendall 32
#> 41   wt qsec -0.22540120      spearman 32
#> 42   wt   vs -0.55491568 pointbiserial 32
#> 43   wt   am -0.69249526 pointbiserial 32
#> 44   wt gear -0.54359562       kendall 32
#> 45   wt carb  0.37137413       kendall 32
#> 46 qsec   vs  0.74453544 pointbiserial 32
#> 47 qsec   am -0.22986086 pointbiserial 32
#> 48 qsec gear -0.09126069       kendall 32
#> 49 qsec carb -0.50643945       kendall 32
#> 50   vs   am  0.16834512           phi 32
#> 51   vs gear  0.20602335 pointbiserial 32
#> 52   vs carb -0.56960714 pointbiserial 32
#> 53   am gear  0.79405876 pointbiserial 32
#> 54   am carb  0.05753435 pointbiserial 32
#> 55 gear carb  0.09801487       kendall 32
res$p_values
#>    var1 var2      p_value
#> 1   mpg  cyl 2.253620e-08
#> 2   mpg disp 1.006950e-09
#> 3   mpg   hp 4.331605e-09
#> 4   mpg drat 5.381347e-05
#> 5   mpg   wt 1.487595e-11
#> 6   mpg qsec 7.055765e-03
#> 7   mpg   vs 3.415937e-05
#> 8   mpg   am 2.850207e-04
#> 9   mpg gear 2.495397e-03
#> 10  mpg carb 2.284747e-04
#> 11  cyl disp 9.956527e-09
#> 12  cyl   hp 3.996666e-08
#> 13  cyl drat 1.139715e-04
#> 14  cyl   wt 2.829114e-07
#> 15  cyl qsec 1.519758e-03
#> 16  cyl   vs 1.843018e-08
#> 17  cyl   am 2.151207e-03
#> 18  cyl gear 1.604329e-03
#> 19  cyl carb 2.682129e-03
#> 20 disp   hp 1.371289e-07
#> 21 disp drat 7.784941e-05
#> 22 disp   wt 3.070059e-09
#> 23 disp qsec 1.627155e-02
#> 24 disp   vs 5.235012e-06
#> 25 disp   am 3.662114e-04
#> 26 disp gear 8.802070e-04
#> 27 disp carb 2.478905e-03
#> 28   hp drat 2.602539e-03
#> 29   hp   wt 1.266240e-06
#> 30   hp qsec 1.737737e-04
#> 31   hp   vs 2.940896e-06
#> 32   hp   am 1.798309e-01
#> 33   hp gear 5.231256e-02
#> 34   hp carb 1.479835e-05
#> 35 drat   wt 7.593194e-07
#> 36 drat qsec 6.170251e-01
#> 37 drat   vs 1.167553e-02
#> 38 drat   am 4.726790e-06
#> 39 drat gear 4.941861e-05
#> 40 drat carb 4.879301e-01
#> 41   wt qsec 2.148388e-01
#> 42   wt   vs 9.798492e-04
#> 43   wt   am 1.125440e-05
#> 44   wt gear 1.413727e-04
#> 45   wt carb 6.507985e-03
#> 46 qsec   vs 1.029669e-06
#> 47 qsec   am 2.056621e-01
#> 48 qsec gear 5.221507e-01
#> 49 qsec carb 2.015873e-04
#> 50   vs   am 3.570439e-01
#> 51   vs gear 2.579439e-01
#> 52   vs carb 6.670496e-04
#> 53   am gear 5.834043e-08
#> 54   am carb 7.544526e-01
#> 55 gear carb 5.302044e-01
res$correlations_masked
#>    var1 var2   estimate   method_used  n
#> 1   mpg  cyl -0.7953134       kendall 32
#> 2   mpg disp -0.7681311       kendall 32
#> 3   mpg   hp -0.7428125       kendall 32
#> 4   mpg drat  0.6514555      spearman 32
#> 5   mpg   wt -0.8864220      spearman 32
#> 6   mpg qsec  0.4669358      spearman 32
#> 7   mpg   vs  0.6640389 pointbiserial 32
#> 8   mpg   am  0.5998324 pointbiserial 32
#> 9   mpg gear  0.4331509       kendall 32
#> 10  mpg carb -0.5043945       kendall 32
#> 11  cyl disp  0.8144263       kendall 32
#> 12  cyl   hp  0.7851865       kendall 32
#> 13  cyl drat -0.5513178       kendall 32
#> 14  cyl   wt  0.7282611       kendall 32
#> 15  cyl qsec -0.4489698       kendall 32
#> 16  cyl   vs -0.8108118 pointbiserial 32
#> 17  cyl   am -0.5226070 pointbiserial 32
#> 18  cyl gear -0.5125435       kendall 32
#> 19  cyl carb  0.4654299       kendall 32
#> 20 disp   hp  0.6659987       kendall 32
#> 21 disp drat -0.4989828       kendall 32
#> 22 disp   wt  0.7433824       kendall 32
#> 23 disp qsec -0.3008155       kendall 32
#> 24 disp   vs -0.7104159 pointbiserial 32
#> 25 disp   am -0.5912270 pointbiserial 32
#> 26 disp gear -0.4759795       kendall 32
#> 27 disp carb  0.4137360       kendall 32
#> 28   hp drat -0.3826269       kendall 32
#> 29   hp   wt  0.6113081       kendall 32
#> 30   hp qsec -0.4729061       kendall 32
#> 31   hp   vs -0.7230967 pointbiserial 32
#> 32   hp   am         NA pointbiserial 32
#> 33   hp gear         NA       kendall 32
#> 34   hp carb  0.5959842       kendall 32
#> 35 drat   wt -0.7503904      spearman 32
#> 36 drat qsec         NA      spearman 32
#> 37 drat   vs  0.4402785 pointbiserial 32
#> 38 drat   am  0.7127111 pointbiserial 32
#> 39 drat gear  0.5839248       kendall 32
#> 40 drat carb         NA       kendall 32
#> 41   wt qsec         NA      spearman 32
#> 42   wt   vs -0.5549157 pointbiserial 32
#> 43   wt   am -0.6924953 pointbiserial 32
#> 44   wt gear -0.5435956       kendall 32
#> 45   wt carb  0.3713741       kendall 32
#> 46 qsec   vs  0.7445354 pointbiserial 32
#> 47 qsec   am         NA pointbiserial 32
#> 48 qsec gear         NA       kendall 32
#> 49 qsec carb -0.5064394       kendall 32
#> 50   vs   am         NA           phi 32
#> 51   vs gear         NA pointbiserial 32
#> 52   vs carb -0.5696071 pointbiserial 32
#> 53   am gear  0.7940588 pointbiserial 32
#> 54   am carb         NA pointbiserial 32
#> 55 gear carb         NA       kendall 32
```
