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
  exact = NULL,
  interpretation = NULL,
  summary_table_only = FALSE,
  p_adjust_method = "BH",
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

- exact:

  Logical or `NULL`. When `method` is `"spearman"` or `"kendall"`,
  controls whether an exact p-value should be computed. Passed to
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html). Default `NULL`
  (let `cor.test` decide).

- interpretation:

  Character string or `NULL`. If not `NULL` and
  `summary_table_only = TRUE`, adds a qualitative interpretation column
  based on `"psychology"`, `"politics"`, or `"medicine"` thresholds.
  Default `NULL`.

- summary_table_only:

  Logical. If `TRUE`, returns only a single data frame with all result
  columns. Default is `FALSE`.

- p_adjust_method:

  P-value adjustment method for multiple testing correction. One of
  `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"` (default),
  `"BY"`, `"fdr"`, or `"none"`. Applied only when
  `summary_table_only = TRUE` and the number of variable pairs is
  greater than 1. An `Adjusted p-value` column is added to the returned
  data frame when applicable. Default `"BH"`.

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

When `summary_table_only = FALSE`, returns a named list of class
`"run_correl"` with the following elements:

- `correlations`:

  A report-ready data frame of pairwise correlations with columns
  `Variable 1`, `Variable 2`, `Correlation`, `Method`, and `n`.

- `correlations_masked`:

  Same as `correlations` but `Correlation` is `NA` for pairs that are
  not significant at `sig_threshold`.

- `p_values`:

  A data frame with columns `Variable 1`, `Variable 2`, and `p-value`.

- `confidence_intervals`:

  A data frame with columns `Variable 1`, `Variable 2`, `95% CI lower`,
  and `95% CI upper`.

- `sig_threshold`:

  The significance threshold used.

- `method`:

  The method argument as supplied.

- `call`:

  The matched call.

When `summary_table_only = TRUE`, returns a data frame of class
`"run_correl"` with columns `Variable 1`, `Variable 2`, `Correlation`,
`95% CI lower`, `95% CI upper`, `p-value`, `Method` and `n`. When the
number of variable pairs is greater than 1, an `Adjusted p-value` column
is added, and an `Interpretation` column when requested.

## Details

**Automatic method selection** (`method = "auto"`):

|  |  |
|----|----|
| **Method** | **Conditions** |
| Pearson | Both continuous and normal |
| Spearman | Ordinal or continuous, non-normal, \\n \ge 30\\, few ties |
| Kendall | Ordinal or continuous, non-normal, \\n \< 30\\ or many ties |
| Point-biserial | One continuous, one binary variable |
| Phi | Both binary |
| Cramer's V | Both nominal (multi-category) factors |
| Polychoric | Both ordered factors (requires polycor) |
| Polyserial | One ordered factor, one continuous (requires polycor) |

Normality is assessed via Shapiro-Wilk (\\n \< 50\\) or Lilliefors (\\n
\ge 50\\; requires nortest). Confidence intervals for Spearman's rho and
Kendall's tau use the Fisher z transformation with the standard errors
of Fieller, Hartley and Pearson (1957).

Confidence intervals for Cramer's V and polychoric/polyserial
correlations are not available and are returned as `NA`. When
`summary_table_only = TRUE`, the function returns a flat data frame
instead of the structured list.

**P-value adjustment:** When `summary_table_only = TRUE` and the number
of variable pairs is greater than 1, p-values are adjusted using the
method specified in `p_adjust_method`, with results stored in an
`Adjusted p-value` column. This accounts for multiple comparisons across
pairs.

**Exact p-values:** For Kendall's tau and Spearman's rho, the `exact`
parameter determines whether an exact p-value is computed. When `NULL`
(default), `cor.test` automatically decides based on sample size and
ties; when `TRUE` or `FALSE`, the choice is forced. See
[`cor.test`](https://rdrr.io/r/stats/cor.test.html) details for the
meaning of `NULL`.

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

Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank correlation
coefficients. I. *Biometrika*, 44(3/4), 470-481.
[doi:10.1093/biomet/44.3-4.470](https://doi.org/10.1093/biomet/44.3-4.470)

## Author

John Lennon L. Calorio

## Examples

``` r
## Basic pairwise correlations (all numeric columns)
run_correl(mtcars)
#>    Variable 1 Variable 2 Correlation         Method  n
#> 1         mpg        cyl -0.79531341        Kendall 32
#> 2         mpg       disp -0.76813115        Kendall 32
#> 3         mpg         hp -0.74281251        Kendall 32
#> 4         mpg       drat  0.68117191        Pearson 32
#> 5         mpg         wt -0.86765938        Pearson 32
#> 6         mpg       qsec  0.41868403        Pearson 32
#> 7         mpg         vs  0.66403892 Point-biserial 32
#> 8         mpg         am  0.59983243 Point-biserial 32
#> 9         mpg       gear  0.43315089        Kendall 32
#> 10        mpg       carb -0.50439455        Kendall 32
#> 11        cyl       disp  0.81442625        Kendall 32
#> 12        cyl         hp  0.78518650        Kendall 32
#> 13        cyl       drat -0.55131785        Kendall 32
#> 14        cyl         wt  0.72826111        Kendall 32
#> 15        cyl       qsec -0.44896982        Kendall 32
#> 16        cyl         vs -0.81081180 Point-biserial 32
#> 17        cyl         am -0.52260705 Point-biserial 32
#> 18        cyl       gear -0.51254349        Kendall 32
#> 19        cyl       carb  0.46542994        Kendall 32
#> 20       disp         hp  0.66599874        Kendall 32
#> 21       disp       drat -0.49898277        Kendall 32
#> 22       disp         wt  0.74338240        Kendall 32
#> 23       disp       qsec -0.30081549        Kendall 32
#> 24       disp         vs -0.71041589 Point-biserial 32
#> 25       disp         am -0.59122704 Point-biserial 32
#> 26       disp       gear -0.47597955        Kendall 32
#> 27       disp       carb  0.41373600        Kendall 32
#> 28         hp       drat -0.38262689        Kendall 32
#> 29         hp         wt  0.61130810        Kendall 32
#> 30         hp       qsec -0.47290613        Kendall 32
#> 31         hp         vs -0.72309674 Point-biserial 32
#> 32         hp         am -0.24320426 Point-biserial 32
#> 33         hp       gear -0.27944584        Kendall 32
#> 34         hp       carb  0.59598416        Kendall 32
#> 35       drat         wt -0.71244065        Pearson 32
#> 36       drat       qsec  0.09120476        Pearson 32
#> 37       drat         vs  0.44027846 Point-biserial 32
#> 38       drat         am  0.71271113 Point-biserial 32
#> 39       drat       gear  0.58392476        Kendall 32
#> 40       drat       carb -0.09535193        Kendall 32
#> 41         wt       qsec -0.17471588        Pearson 32
#> 42         wt         vs -0.55491568 Point-biserial 32
#> 43         wt         am -0.69249526 Point-biserial 32
#> 44         wt       gear -0.54359562        Kendall 32
#> 45         wt       carb  0.37137413        Kendall 32
#> 46       qsec         vs  0.74453544 Point-biserial 32
#> 47       qsec         am -0.22986086 Point-biserial 32
#> 48       qsec       gear -0.09126069        Kendall 32
#> 49       qsec       carb -0.50643945        Kendall 32
#> 50         vs         am  0.16834512            Phi 32
#> 51         vs       gear  0.20602335 Point-biserial 32
#> 52         vs       carb -0.56960714 Point-biserial 32
#> 53         am       gear  0.79405876 Point-biserial 32
#> 54         am       carb  0.05753435 Point-biserial 32
#> 55       gear       carb  0.09801487        Kendall 32

## Correlate one variable against all others
run_correl(mtcars, y = "mpg")
#>    Variable 1 Variable 2 Correlation         Method  n
#> 1         mpg        cyl  -0.7953134        Kendall 32
#> 2         mpg       disp  -0.7681311        Kendall 32
#> 3         mpg         hp  -0.7428125        Kendall 32
#> 4         mpg       drat   0.6811719        Pearson 32
#> 5         mpg         wt  -0.8676594        Pearson 32
#> 6         mpg       qsec   0.4186840        Pearson 32
#> 7         mpg         vs   0.6640389 Point-biserial 32
#> 8         mpg         am   0.5998324 Point-biserial 32
#> 9         mpg       gear   0.4331509        Kendall 32
#> 10        mpg       carb  -0.5043945        Kendall 32

## Single pair
run_correl(mtcars, y = "mpg", z = "wt")
#>   Variable 1 Variable 2 Correlation  Method  n
#> 1        mpg         wt  -0.8676594 Pearson 32

## Specify method explicitly
run_correl(mtcars, method = "spearman")
#>    Variable 1 Variable 2 Correlation   Method  n
#> 1         mpg        cyl -0.91080131 Spearman 32
#> 2         mpg       disp -0.90888236 Spearman 32
#> 3         mpg         hp -0.89466465 Spearman 32
#> 4         mpg       drat  0.65145546 Spearman 32
#> 5         mpg         wt -0.88642203 Spearman 32
#> 6         mpg       qsec  0.46693575 Spearman 32
#> 7         mpg         vs  0.70659679 Spearman 32
#> 8         mpg         am  0.56200569 Spearman 32
#> 9         mpg       gear  0.54278158 Spearman 32
#> 10        mpg       carb -0.65749764 Spearman 32
#> 11        cyl       disp  0.92765158 Spearman 32
#> 12        cyl         hp  0.90179094 Spearman 32
#> 13        cyl       drat -0.67888119 Spearman 32
#> 14        cyl         wt  0.85772816 Spearman 32
#> 15        cyl       qsec -0.57235095 Spearman 32
#> 16        cyl         vs -0.81378895 Spearman 32
#> 17        cyl         am -0.52207118 Spearman 32
#> 18        cyl       gear -0.56431047 Spearman 32
#> 19        cyl       carb  0.58006798 Spearman 32
#> 20       disp         hp  0.85104263 Spearman 32
#> 21       disp       drat -0.68359210 Spearman 32
#> 22       disp         wt  0.89770644 Spearman 32
#> 23       disp       qsec -0.45978176 Spearman 32
#> 24       disp         vs -0.72366435 Spearman 32
#> 25       disp         am -0.62406767 Spearman 32
#> 26       disp       gear -0.59447030 Spearman 32
#> 27       disp       carb  0.53977806 Spearman 32
#> 28         hp       drat -0.52012499 Spearman 32
#> 29         hp         wt  0.77467673 Spearman 32
#> 30         hp       qsec -0.66660602 Spearman 32
#> 31         hp         vs -0.75159339 Spearman 32
#> 32         hp         am -0.36232756 Spearman 32
#> 33         hp       gear -0.33140155 Spearman 32
#> 34         hp       carb  0.73337937 Spearman 32
#> 35       drat         wt -0.75039041 Spearman 32
#> 36       drat       qsec  0.09186863 Spearman 32
#> 37       drat         vs  0.44745745 Spearman 32
#> 38       drat         am  0.68657079 Spearman 32
#> 39       drat       gear  0.74481617 Spearman 32
#> 40       drat       carb -0.12522294 Spearman 32
#> 41         wt       qsec -0.22540120 Spearman 32
#> 42         wt         vs -0.58701619 Spearman 32
#> 43         wt         am -0.73771259 Spearman 32
#> 44         wt       gear -0.67612839 Spearman 32
#> 45         wt       carb  0.49981205 Spearman 32
#> 46       qsec         vs  0.79157148 Spearman 32
#> 47       qsec         am -0.20333211 Spearman 32
#> 48       qsec       gear -0.14819967 Spearman 32
#> 49       qsec       carb -0.65871814 Spearman 32
#> 50         vs         am  0.16834512 Spearman 32
#> 51         vs       gear  0.28266170 Spearman 32
#> 52         vs       carb -0.63369482 Spearman 32
#> 53         am       gear  0.80768800 Spearman 32
#> 54         am       carb -0.06436525 Spearman 32
#> 55       gear       carb  0.11488698 Spearman 32

## Exclude specific levels before analysis
df <- mtcars
df$cyl <- as.numeric(factor(df$cyl))
run_correl(df, remove = list("cyl" ~ c("4", "6")))
#>    Variable 1 Variable 2 Correlation         Method  n
#> 1         mpg        cyl -0.79531341        Kendall 32
#> 2         mpg       disp -0.76813115        Kendall 32
#> 3         mpg         hp -0.74281251        Kendall 32
#> 4         mpg       drat  0.68117191        Pearson 32
#> 5         mpg         wt -0.86765938        Pearson 32
#> 6         mpg       qsec  0.41868403        Pearson 32
#> 7         mpg         vs  0.66403892 Point-biserial 32
#> 8         mpg         am  0.59983243 Point-biserial 32
#> 9         mpg       gear  0.43315089        Kendall 32
#> 10        mpg       carb -0.50439455        Kendall 32
#> 11        cyl       disp  0.81442625        Kendall 32
#> 12        cyl         hp  0.78518650        Kendall 32
#> 13        cyl       drat -0.55131785        Kendall 32
#> 14        cyl         wt  0.72826111        Kendall 32
#> 15        cyl       qsec -0.44896982        Kendall 32
#> 16        cyl         vs -0.81081180 Point-biserial 32
#> 17        cyl         am -0.52260705 Point-biserial 32
#> 18        cyl       gear -0.51254349        Kendall 32
#> 19        cyl       carb  0.46542994        Kendall 32
#> 20       disp         hp  0.66599874        Kendall 32
#> 21       disp       drat -0.49898277        Kendall 32
#> 22       disp         wt  0.74338240        Kendall 32
#> 23       disp       qsec -0.30081549        Kendall 32
#> 24       disp         vs -0.71041589 Point-biserial 32
#> 25       disp         am -0.59122704 Point-biserial 32
#> 26       disp       gear -0.47597955        Kendall 32
#> 27       disp       carb  0.41373600        Kendall 32
#> 28         hp       drat -0.38262689        Kendall 32
#> 29         hp         wt  0.61130810        Kendall 32
#> 30         hp       qsec -0.47290613        Kendall 32
#> 31         hp         vs -0.72309674 Point-biserial 32
#> 32         hp         am -0.24320426 Point-biserial 32
#> 33         hp       gear -0.27944584        Kendall 32
#> 34         hp       carb  0.59598416        Kendall 32
#> 35       drat         wt -0.71244065        Pearson 32
#> 36       drat       qsec  0.09120476        Pearson 32
#> 37       drat         vs  0.44027846 Point-biserial 32
#> 38       drat         am  0.71271113 Point-biserial 32
#> 39       drat       gear  0.58392476        Kendall 32
#> 40       drat       carb -0.09535193        Kendall 32
#> 41         wt       qsec -0.17471588        Pearson 32
#> 42         wt         vs -0.55491568 Point-biserial 32
#> 43         wt         am -0.69249526 Point-biserial 32
#> 44         wt       gear -0.54359562        Kendall 32
#> 45         wt       carb  0.37137413        Kendall 32
#> 46       qsec         vs  0.74453544 Point-biserial 32
#> 47       qsec         am -0.22986086 Point-biserial 32
#> 48       qsec       gear -0.09126069        Kendall 32
#> 49       qsec       carb -0.50643945        Kendall 32
#> 50         vs         am  0.16834512            Phi 32
#> 51         vs       gear  0.20602335 Point-biserial 32
#> 52         vs       carb -0.56960714 Point-biserial 32
#> 53         am       gear  0.79405876 Point-biserial 32
#> 54         am       carb  0.05753435 Point-biserial 32
#> 55       gear       carb  0.09801487        Kendall 32

## Access individual output tables
res <- run_correl(mtcars)
res$correlations
#>    Variable 1 Variable 2 Correlation         Method  n
#> 1         mpg        cyl -0.79531341        Kendall 32
#> 2         mpg       disp -0.76813115        Kendall 32
#> 3         mpg         hp -0.74281251        Kendall 32
#> 4         mpg       drat  0.68117191        Pearson 32
#> 5         mpg         wt -0.86765938        Pearson 32
#> 6         mpg       qsec  0.41868403        Pearson 32
#> 7         mpg         vs  0.66403892 Point-biserial 32
#> 8         mpg         am  0.59983243 Point-biserial 32
#> 9         mpg       gear  0.43315089        Kendall 32
#> 10        mpg       carb -0.50439455        Kendall 32
#> 11        cyl       disp  0.81442625        Kendall 32
#> 12        cyl         hp  0.78518650        Kendall 32
#> 13        cyl       drat -0.55131785        Kendall 32
#> 14        cyl         wt  0.72826111        Kendall 32
#> 15        cyl       qsec -0.44896982        Kendall 32
#> 16        cyl         vs -0.81081180 Point-biserial 32
#> 17        cyl         am -0.52260705 Point-biserial 32
#> 18        cyl       gear -0.51254349        Kendall 32
#> 19        cyl       carb  0.46542994        Kendall 32
#> 20       disp         hp  0.66599874        Kendall 32
#> 21       disp       drat -0.49898277        Kendall 32
#> 22       disp         wt  0.74338240        Kendall 32
#> 23       disp       qsec -0.30081549        Kendall 32
#> 24       disp         vs -0.71041589 Point-biserial 32
#> 25       disp         am -0.59122704 Point-biserial 32
#> 26       disp       gear -0.47597955        Kendall 32
#> 27       disp       carb  0.41373600        Kendall 32
#> 28         hp       drat -0.38262689        Kendall 32
#> 29         hp         wt  0.61130810        Kendall 32
#> 30         hp       qsec -0.47290613        Kendall 32
#> 31         hp         vs -0.72309674 Point-biserial 32
#> 32         hp         am -0.24320426 Point-biserial 32
#> 33         hp       gear -0.27944584        Kendall 32
#> 34         hp       carb  0.59598416        Kendall 32
#> 35       drat         wt -0.71244065        Pearson 32
#> 36       drat       qsec  0.09120476        Pearson 32
#> 37       drat         vs  0.44027846 Point-biserial 32
#> 38       drat         am  0.71271113 Point-biserial 32
#> 39       drat       gear  0.58392476        Kendall 32
#> 40       drat       carb -0.09535193        Kendall 32
#> 41         wt       qsec -0.17471588        Pearson 32
#> 42         wt         vs -0.55491568 Point-biserial 32
#> 43         wt         am -0.69249526 Point-biserial 32
#> 44         wt       gear -0.54359562        Kendall 32
#> 45         wt       carb  0.37137413        Kendall 32
#> 46       qsec         vs  0.74453544 Point-biserial 32
#> 47       qsec         am -0.22986086 Point-biserial 32
#> 48       qsec       gear -0.09126069        Kendall 32
#> 49       qsec       carb -0.50643945        Kendall 32
#> 50         vs         am  0.16834512            Phi 32
#> 51         vs       gear  0.20602335 Point-biserial 32
#> 52         vs       carb -0.56960714 Point-biserial 32
#> 53         am       gear  0.79405876 Point-biserial 32
#> 54         am       carb  0.05753435 Point-biserial 32
#> 55       gear       carb  0.09801487        Kendall 32
res$p_values
#>    Variable 1 Variable 2      p-value
#> 1         mpg        cyl 2.253620e-08
#> 2         mpg       disp 1.006950e-09
#> 3         mpg         hp 4.331605e-09
#> 4         mpg       drat 1.776240e-05
#> 5         mpg         wt 1.293959e-10
#> 6         mpg       qsec 1.708199e-02
#> 7         mpg         vs 3.415937e-05
#> 8         mpg         am 2.850207e-04
#> 9         mpg       gear 2.495397e-03
#> 10        mpg       carb 2.284747e-04
#> 11        cyl       disp 9.956527e-09
#> 12        cyl         hp 3.996666e-08
#> 13        cyl       drat 1.139715e-04
#> 14        cyl         wt 2.829114e-07
#> 15        cyl       qsec 1.519758e-03
#> 16        cyl         vs 1.843018e-08
#> 17        cyl         am 2.151207e-03
#> 18        cyl       gear 1.604329e-03
#> 19        cyl       carb 2.682129e-03
#> 20       disp         hp 1.371289e-07
#> 21       disp       drat 7.784941e-05
#> 22       disp         wt 3.070059e-09
#> 23       disp       qsec 1.627155e-02
#> 24       disp         vs 5.235012e-06
#> 25       disp         am 3.662114e-04
#> 26       disp       gear 8.802070e-04
#> 27       disp       carb 2.478905e-03
#> 28         hp       drat 2.602539e-03
#> 29         hp         wt 1.266240e-06
#> 30         hp       qsec 1.737737e-04
#> 31         hp         vs 2.940896e-06
#> 32         hp         am 1.798309e-01
#> 33         hp       gear 5.231256e-02
#> 34         hp       carb 1.479835e-05
#> 35       drat         wt 4.784260e-06
#> 36       drat       qsec 6.195826e-01
#> 37       drat         vs 1.167553e-02
#> 38       drat         am 4.726790e-06
#> 39       drat       gear 4.941861e-05
#> 40       drat       carb 4.879301e-01
#> 41         wt       qsec 3.388683e-01
#> 42         wt         vs 9.798492e-04
#> 43         wt         am 1.125440e-05
#> 44         wt       gear 1.413727e-04
#> 45         wt       carb 6.507985e-03
#> 46       qsec         vs 1.029669e-06
#> 47       qsec         am 2.056621e-01
#> 48       qsec       gear 5.221507e-01
#> 49       qsec       carb 2.015873e-04
#> 50         vs         am 3.570439e-01
#> 51         vs       gear 2.579439e-01
#> 52         vs       carb 6.670496e-04
#> 53         am       gear 5.834043e-08
#> 54         am       carb 7.544526e-01
#> 55       gear       carb 5.302044e-01
res$correlations_masked
#>    Variable 1 Variable 2 Correlation         Method  n
#> 1         mpg        cyl  -0.7953134        Kendall 32
#> 2         mpg       disp  -0.7681311        Kendall 32
#> 3         mpg         hp  -0.7428125        Kendall 32
#> 4         mpg       drat   0.6811719        Pearson 32
#> 5         mpg         wt  -0.8676594        Pearson 32
#> 6         mpg       qsec   0.4186840        Pearson 32
#> 7         mpg         vs   0.6640389 Point-biserial 32
#> 8         mpg         am   0.5998324 Point-biserial 32
#> 9         mpg       gear   0.4331509        Kendall 32
#> 10        mpg       carb  -0.5043945        Kendall 32
#> 11        cyl       disp   0.8144263        Kendall 32
#> 12        cyl         hp   0.7851865        Kendall 32
#> 13        cyl       drat  -0.5513178        Kendall 32
#> 14        cyl         wt   0.7282611        Kendall 32
#> 15        cyl       qsec  -0.4489698        Kendall 32
#> 16        cyl         vs  -0.8108118 Point-biserial 32
#> 17        cyl         am  -0.5226070 Point-biserial 32
#> 18        cyl       gear  -0.5125435        Kendall 32
#> 19        cyl       carb   0.4654299        Kendall 32
#> 20       disp         hp   0.6659987        Kendall 32
#> 21       disp       drat  -0.4989828        Kendall 32
#> 22       disp         wt   0.7433824        Kendall 32
#> 23       disp       qsec  -0.3008155        Kendall 32
#> 24       disp         vs  -0.7104159 Point-biserial 32
#> 25       disp         am  -0.5912270 Point-biserial 32
#> 26       disp       gear  -0.4759795        Kendall 32
#> 27       disp       carb   0.4137360        Kendall 32
#> 28         hp       drat  -0.3826269        Kendall 32
#> 29         hp         wt   0.6113081        Kendall 32
#> 30         hp       qsec  -0.4729061        Kendall 32
#> 31         hp         vs  -0.7230967 Point-biserial 32
#> 32         hp         am          NA Point-biserial 32
#> 33         hp       gear          NA        Kendall 32
#> 34         hp       carb   0.5959842        Kendall 32
#> 35       drat         wt  -0.7124406        Pearson 32
#> 36       drat       qsec          NA        Pearson 32
#> 37       drat         vs   0.4402785 Point-biserial 32
#> 38       drat         am   0.7127111 Point-biserial 32
#> 39       drat       gear   0.5839248        Kendall 32
#> 40       drat       carb          NA        Kendall 32
#> 41         wt       qsec          NA        Pearson 32
#> 42         wt         vs  -0.5549157 Point-biserial 32
#> 43         wt         am  -0.6924953 Point-biserial 32
#> 44         wt       gear  -0.5435956        Kendall 32
#> 45         wt       carb   0.3713741        Kendall 32
#> 46       qsec         vs   0.7445354 Point-biserial 32
#> 47       qsec         am          NA Point-biserial 32
#> 48       qsec       gear          NA        Kendall 32
#> 49       qsec       carb  -0.5064394        Kendall 32
#> 50         vs         am          NA            Phi 32
#> 51         vs       gear          NA Point-biserial 32
#> 52         vs       carb  -0.5696071 Point-biserial 32
#> 53         am       gear   0.7940588 Point-biserial 32
#> 54         am       carb          NA Point-biserial 32
#> 55       gear       carb          NA        Kendall 32
```
