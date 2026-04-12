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
  sig_threshold = 0.05,
  normality_method = "auto",
  alpha_normality = 0.05,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A data frame, tibble, or matrix whose columns are the variables to
  correlate. Column names may contain special characters.

- y:

  A single character string naming a column in `x` to use as the primary
  variable. When both `y` and `z` are `NULL`, all pairwise correlations
  among numeric columns of `x` are computed. When only `y` is supplied,
  correlations of `y` against all other columns are returned.

- z:

  A single character string naming a second column in `x`. When supplied
  together with `y`, only the correlation between `y` and `z` is
  computed.

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

- verbose:

  Logical. When `TRUE`, diagnostic messages about assumption checks and
  method selection are printed. Default `FALSE`.

- ...:

  Additional arguments passed to
  [`cor.test`](https://rdrr.io/r/stats/cor.test.html).

## Value

A named list of class `"pondeR_correl"` with the following elements:

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
correlations are not available and are returned as `NA`.

## Author

John Lennon L. Calorio

## Examples

``` r
## Basic pairwise correlations (all numeric columns)
run_correl(mtcars)
#> $correlations
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
#> 
#> $correlations_masked
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
#> 
#> $p_values
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
#> 
#> $confidence_intervals
#>    var1 var2    ci_lower    ci_upper
#> 1   mpg  cyl -0.89563030 -0.61800527
#> 2   mpg disp -0.88088618 -0.57287850
#> 3   mpg   hp -0.86695324 -0.53188409
#> 4   mpg drat  0.39175113  0.81501237
#> 5   mpg   wt -0.94349203 -0.77828961
#> 6   mpg qsec  0.14123780  0.70142529
#> 7   mpg   vs  0.41036301  0.82232624
#> 8   mpg   am  0.31755830  0.78445202
#> 9   mpg gear  0.09948247  0.67925253
#> 10  mpg carb -0.72548984 -0.18892934
#> 11  cyl disp  0.65045093  0.90586708
#> 12  cyl   hp  0.60105535  0.89016294
#> 13  cyl drat -0.75489071 -0.25084701
#> 14  cyl   wt  0.50876240  0.85885677
#> 15  cyl qsec -0.68969061 -0.11888816
#> 16  cyl   vs -0.90393935 -0.64426886
#> 17  cyl   am -0.73699794 -0.21266745
#> 18  cyl gear -0.73065432 -0.19950526
#> 19  cyl carb  0.13935172  0.70044660
#> 20 disp   hp  0.41328088  0.82346065
#> 21 disp drat -0.72204626 -0.18194600
#> 22 disp   wt  0.53279607  0.86726900
#> 23 disp qsec -0.58784853  0.05348908
#> 24 disp   vs -0.84883771 -0.48083271
#> 25 disp   am -0.77926901 -0.30551779
#> 26 disp gear -0.70728464 -0.15261532
#> 27 disp carb  0.07600741  0.66630411
#> 28   hp drat -0.64523394 -0.03915727
#> 29   hp   wt  0.33375774  0.79132390
#> 30   hp qsec -0.70529697 -0.14873917
#> 31   hp   vs -0.85596751 -0.50063181
#> 32   hp   am -0.54562696  0.11526460
#> 33   hp gear -0.57236769  0.07672422
#> 34   hp carb  0.31216260  0.78213740
#> 35 drat   wt -0.87114384 -0.54405098
#> 36 drat qsec -0.26532476  0.42688767
#> 37 drat   vs  0.10819483  0.68396800
#> 38 drat   am  0.48439908  0.85013194
#> 39 drat gear  0.29537173  0.77485059
#> 40 drat carb -0.42975706  0.26205501
#> 41   wt qsec -0.53226152  0.13380968
#> 42   wt   vs -0.75711174 -0.25569823
#> 43   wt   am -0.83867523 -0.45324615
#> 44   wt gear -0.75010777 -0.24048527
#> 45   wt carb  0.02605405  0.63751279
#> 46 qsec   vs  0.53464277  0.86790757
#> 47 qsec   am -0.53562398  0.12918764
#> 48 qsec gear -0.42638618  0.26589459
#> 49 qsec carb -0.72678816 -0.19157641
#> 50   vs   am -0.19159569  0.48837125
#> 51   vs gear -0.15371324  0.51753786
#> 52   vs carb -0.76613289 -0.27566545
#> 53   am gear  0.61589632  0.89495457
#> 54   am carb -0.29712041  0.39823890
#> 55 gear carb -0.25954979  0.43194613
#> 
#> $sig_threshold
#> [1] 0.05
#> 
#> $method
#> [1] "auto"
#> 
#> $call
#> run_correl(x = mtcars)
#> 
#> attr(,"class")
#> [1] "pondeR_correl" "list"         

## Correlate one variable against all others
run_correl(mtcars, y = "mpg")
#> $correlations
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
#> 
#> $correlations_masked
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
#> 
#> $p_values
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
#> 
#> $confidence_intervals
#>    var1 var2    ci_lower   ci_upper
#> 1   mpg  cyl -0.89563030 -0.6180053
#> 2   mpg disp -0.88088618 -0.5728785
#> 3   mpg   hp -0.86695324 -0.5318841
#> 4   mpg drat  0.39175113  0.8150124
#> 5   mpg   wt -0.94349203 -0.7782896
#> 6   mpg qsec  0.14123780  0.7014253
#> 7   mpg   vs  0.41036301  0.8223262
#> 8   mpg   am  0.31755830  0.7844520
#> 9   mpg gear  0.09948247  0.6792525
#> 10  mpg carb -0.72548984 -0.1889293
#> 
#> $sig_threshold
#> [1] 0.05
#> 
#> $method
#> [1] "auto"
#> 
#> $call
#> run_correl(x = mtcars, y = "mpg")
#> 
#> attr(,"class")
#> [1] "pondeR_correl" "list"         

## Single pair
run_correl(mtcars, y = "mpg", z = "wt")
#> $correlations
#>   var1 var2  estimate method_used  n
#> 1  mpg   wt -0.886422    spearman 32
#> 
#> $correlations_masked
#>   var1 var2  estimate method_used  n
#> 1  mpg   wt -0.886422    spearman 32
#> 
#> $p_values
#>   var1 var2      p_value
#> 1  mpg   wt 1.487595e-11
#> 
#> $confidence_intervals
#>   var1 var2  ci_lower   ci_upper
#> 1  mpg   wt -0.943492 -0.7782896
#> 
#> $sig_threshold
#> [1] 0.05
#> 
#> $method
#> [1] "auto"
#> 
#> $call
#> run_correl(x = mtcars, y = "mpg", z = "wt")
#> 
#> attr(,"class")
#> [1] "pondeR_correl" "list"         

## Specify method explicitly
run_correl(mtcars, method = "spearman")
#> $correlations
#>    var1 var2    estimate method_used  n
#> 1   mpg  cyl -0.91080131    spearman 32
#> 2   mpg disp -0.90888236    spearman 32
#> 3   mpg   hp -0.89466465    spearman 32
#> 4   mpg drat  0.65145546    spearman 32
#> 5   mpg   wt -0.88642203    spearman 32
#> 6   mpg qsec  0.46693575    spearman 32
#> 7   mpg   vs  0.70659679    spearman 32
#> 8   mpg   am  0.56200569    spearman 32
#> 9   mpg gear  0.54278158    spearman 32
#> 10  mpg carb -0.65749764    spearman 32
#> 11  cyl disp  0.92765158    spearman 32
#> 12  cyl   hp  0.90179094    spearman 32
#> 13  cyl drat -0.67888119    spearman 32
#> 14  cyl   wt  0.85772816    spearman 32
#> 15  cyl qsec -0.57235095    spearman 32
#> 16  cyl   vs -0.81378895    spearman 32
#> 17  cyl   am -0.52207118    spearman 32
#> 18  cyl gear -0.56431047    spearman 32
#> 19  cyl carb  0.58006798    spearman 32
#> 20 disp   hp  0.85104263    spearman 32
#> 21 disp drat -0.68359210    spearman 32
#> 22 disp   wt  0.89770644    spearman 32
#> 23 disp qsec -0.45978176    spearman 32
#> 24 disp   vs -0.72366435    spearman 32
#> 25 disp   am -0.62406767    spearman 32
#> 26 disp gear -0.59447030    spearman 32
#> 27 disp carb  0.53977806    spearman 32
#> 28   hp drat -0.52012499    spearman 32
#> 29   hp   wt  0.77467673    spearman 32
#> 30   hp qsec -0.66660602    spearman 32
#> 31   hp   vs -0.75159339    spearman 32
#> 32   hp   am -0.36232756    spearman 32
#> 33   hp gear -0.33140155    spearman 32
#> 34   hp carb  0.73337937    spearman 32
#> 35 drat   wt -0.75039041    spearman 32
#> 36 drat qsec  0.09186863    spearman 32
#> 37 drat   vs  0.44745745    spearman 32
#> 38 drat   am  0.68657079    spearman 32
#> 39 drat gear  0.74481617    spearman 32
#> 40 drat carb -0.12522294    spearman 32
#> 41   wt qsec -0.22540120    spearman 32
#> 42   wt   vs -0.58701619    spearman 32
#> 43   wt   am -0.73771259    spearman 32
#> 44   wt gear -0.67612839    spearman 32
#> 45   wt carb  0.49981205    spearman 32
#> 46 qsec   vs  0.79157148    spearman 32
#> 47 qsec   am -0.20333211    spearman 32
#> 48 qsec gear -0.14819967    spearman 32
#> 49 qsec carb -0.65871814    spearman 32
#> 50   vs   am  0.16834512    spearman 32
#> 51   vs gear  0.28266170    spearman 32
#> 52   vs carb -0.63369482    spearman 32
#> 53   am gear  0.80768800    spearman 32
#> 54   am carb -0.06436525    spearman 32
#> 55 gear carb  0.11488698    spearman 32
#> 
#> $correlations_masked
#>    var1 var2   estimate method_used  n
#> 1   mpg  cyl -0.9108013    spearman 32
#> 2   mpg disp -0.9088824    spearman 32
#> 3   mpg   hp -0.8946646    spearman 32
#> 4   mpg drat  0.6514555    spearman 32
#> 5   mpg   wt -0.8864220    spearman 32
#> 6   mpg qsec  0.4669358    spearman 32
#> 7   mpg   vs  0.7065968    spearman 32
#> 8   mpg   am  0.5620057    spearman 32
#> 9   mpg gear  0.5427816    spearman 32
#> 10  mpg carb -0.6574976    spearman 32
#> 11  cyl disp  0.9276516    spearman 32
#> 12  cyl   hp  0.9017909    spearman 32
#> 13  cyl drat -0.6788812    spearman 32
#> 14  cyl   wt  0.8577282    spearman 32
#> 15  cyl qsec -0.5723509    spearman 32
#> 16  cyl   vs -0.8137890    spearman 32
#> 17  cyl   am -0.5220712    spearman 32
#> 18  cyl gear -0.5643105    spearman 32
#> 19  cyl carb  0.5800680    spearman 32
#> 20 disp   hp  0.8510426    spearman 32
#> 21 disp drat -0.6835921    spearman 32
#> 22 disp   wt  0.8977064    spearman 32
#> 23 disp qsec -0.4597818    spearman 32
#> 24 disp   vs -0.7236643    spearman 32
#> 25 disp   am -0.6240677    spearman 32
#> 26 disp gear -0.5944703    spearman 32
#> 27 disp carb  0.5397781    spearman 32
#> 28   hp drat -0.5201250    spearman 32
#> 29   hp   wt  0.7746767    spearman 32
#> 30   hp qsec -0.6666060    spearman 32
#> 31   hp   vs -0.7515934    spearman 32
#> 32   hp   am -0.3623276    spearman 32
#> 33   hp gear         NA    spearman 32
#> 34   hp carb  0.7333794    spearman 32
#> 35 drat   wt -0.7503904    spearman 32
#> 36 drat qsec         NA    spearman 32
#> 37 drat   vs  0.4474575    spearman 32
#> 38 drat   am  0.6865708    spearman 32
#> 39 drat gear  0.7448162    spearman 32
#> 40 drat carb         NA    spearman 32
#> 41   wt qsec         NA    spearman 32
#> 42   wt   vs -0.5870162    spearman 32
#> 43   wt   am -0.7377126    spearman 32
#> 44   wt gear -0.6761284    spearman 32
#> 45   wt carb  0.4998120    spearman 32
#> 46 qsec   vs  0.7915715    spearman 32
#> 47 qsec   am         NA    spearman 32
#> 48 qsec gear         NA    spearman 32
#> 49 qsec carb -0.6587181    spearman 32
#> 50   vs   am         NA    spearman 32
#> 51   vs gear         NA    spearman 32
#> 52   vs carb -0.6336948    spearman 32
#> 53   am gear  0.8076880    spearman 32
#> 54   am carb         NA    spearman 32
#> 55 gear carb         NA    spearman 32
#> 
#> $p_values
#>    var1 var2      p_value
#> 1   mpg  cyl 4.690287e-13
#> 2   mpg disp 6.370336e-13
#> 3   mpg   hp 5.085969e-12
#> 4   mpg drat 5.381347e-05
#> 5   mpg   wt 1.487595e-11
#> 6   mpg qsec 7.055765e-03
#> 7   mpg   vs 6.191450e-06
#> 8   mpg   am 8.156989e-04
#> 9   mpg gear 1.328681e-03
#> 10  mpg carb 4.337570e-05
#> 11  cyl disp 2.275443e-14
#> 12  cyl   hp 1.867686e-12
#> 13  cyl drat 1.943342e-05
#> 14  cyl   wt 3.574157e-10
#> 15  cyl qsec 6.195832e-04
#> 16  cyl   vs 1.484058e-08
#> 17  cyl   am 2.178046e-03
#> 18  cyl gear 7.678209e-04
#> 19  cyl carb 5.016643e-04
#> 20 disp   hp 6.791338e-10
#> 21 disp drat 1.613884e-05
#> 22 disp   wt 3.346362e-12
#> 23 disp qsec 8.108019e-03
#> 24 disp   vs 2.863870e-06
#> 25 disp   am 1.352011e-04
#> 26 disp gear 3.334775e-04
#> 27 disp carb 1.430209e-03
#> 28   hp drat 2.277988e-03
#> 29   hp   wt 1.953795e-07
#> 30   hp qsec 3.105344e-05
#> 31   hp   vs 7.125286e-07
#> 32   hp   am 4.155768e-02
#> 33   hp gear 6.390322e-02
#> 34   hp carb 1.799847e-06
#> 35 drat   wt 7.593194e-07
#> 36 drat qsec 6.170251e-01
#> 37 drat   vs 1.023343e-02
#> 38 drat   am 1.432515e-05
#> 39 drat gear 1.014930e-06
#> 40 drat carb 4.946824e-01
#> 41   wt qsec 2.148388e-01
#> 42   wt   vs 4.129434e-04
#> 43   wt   am 1.453656e-06
#> 44   wt gear 2.162837e-05
#> 45   wt carb 3.583063e-03
#> 46 qsec   vs 6.860828e-08
#> 47 qsec   am 2.643506e-01
#> 48 qsec gear 4.182425e-01
#> 49 qsec carb 4.150300e-05
#> 50   vs   am 3.570439e-01
#> 51   vs gear 1.169934e-01
#> 52   vs carb 9.878823e-05
#> 53   am gear 2.304063e-08
#> 54   am carb 7.263524e-01
#> 55 gear carb 5.312358e-01
#> 
#> $confidence_intervals
#>    var1 var2    ci_lower    ci_upper
#> 1   mpg  cyl -0.95590768 -0.82371024
#> 2   mpg disp -0.95493623 -0.82009409
#> 3   mpg   hp -0.94770776 -0.79352070
#> 4   mpg drat  0.39175113  0.81501237
#> 5   mpg   wt -0.94349203 -0.77828961
#> 6   mpg qsec  0.14123780  0.70142529
#> 7   mpg   vs  0.47491533  0.84668051
#> 8   mpg   am  0.26530253  0.76147494
#> 9   mpg gear  0.23939700  0.74960232
#> 10  mpg carb -0.81853078 -0.40066152
#> 11  cyl disp  0.85577081  0.96439580
#> 12  cyl   hp  0.80679193  0.95133767
#> 13  cyl drat -0.83088628 -0.43259090
#> 14  cyl   wt  0.72624203  0.92867092
#> 15  cyl qsec -0.76780916 -0.27942289
#> 16  cyl   vs -0.90552746 -0.64935933
#> 17  cyl   am -0.73666110 -0.21196374
#> 18  cyl gear -0.76288943 -0.26843737
#> 19  cyl carb  0.29003912  0.77250940
#> 20 disp   hp  0.71432792  0.92518482
#> 21 disp drat -0.83358827 -0.43970916
#> 22 disp   wt  0.79917366  0.94925883
#> 23 disp qsec -0.69676775 -0.13229826
#> 24 disp   vs -0.85628547 -0.50152352
#> 25 disp   am -0.79891158 -0.35196414
#> 26 disp gear -0.78122545 -0.31004500
#> 27 disp carb  0.23538825  0.74773528
#> 28   hp drat -0.73543686 -0.20941067
#> 29   hp   wt  0.58363805  0.88445673
#> 30   hp qsec -0.82381191 -0.41418609
#> 31   hp   vs -0.87180747 -0.54599046
#> 32   hp   am -0.63126642 -0.01560516
#> 33   hp gear -0.60964620  0.01955178
#> 34   hp carb  0.51685919  0.86171205
#> 35 drat   wt -0.87114384 -0.54405098
#> 36 drat qsec -0.26532476  0.42688767
#> 37 drat   vs  0.11702191  0.68869698
#> 38 drat   am  0.44422587  0.83529302
#> 39 drat gear  0.53509267  0.86806298
#> 40 drat carb -0.45408933  0.23367430
#> 41   wt qsec -0.53226152  0.13380968
#> 42   wt   vs -0.77672341 -0.29965915
#> 43   wt   am -0.86412306 -0.52374454
#> 44   wt gear -0.82930406 -0.42844555
#> 45   wt carb  0.18301403  0.72257466
#> 46 qsec   vs  0.61172296  0.89361361
#> 47 qsec   am -0.51547831  0.15645460
#> 48 qsec gear -0.47247773  0.21141969
#> 49 qsec carb -0.81924002 -0.40246730
#> 50   vs   am -0.19159569  0.48837125
#> 51   vs gear -0.07325218  0.57471078
#> 52   vs carb -0.80459992 -0.36583858
#> 53   am gear  0.63894348  0.90227025
#> 54   am carb -0.40399211  0.29085665
#> 55 gear carb -0.24356427  0.44572503
#> 
#> $sig_threshold
#> [1] 0.05
#> 
#> $method
#> [1] "spearman"
#> 
#> $call
#> run_correl(x = mtcars, method = "spearman")
#> 
#> attr(,"class")
#> [1] "pondeR_correl" "list"         

## Exclude specific levels before analysis
df <- mtcars
df$cyl <- as.numeric(factor(df$cyl))
run_correl(df, remove = list("cyl" ~ c("4", "6")))
#> $correlations
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
#> 
#> $correlations_masked
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
#> 
#> $p_values
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
#> 
#> $confidence_intervals
#>    var1 var2    ci_lower    ci_upper
#> 1   mpg  cyl -0.89563030 -0.61800527
#> 2   mpg disp -0.88088618 -0.57287850
#> 3   mpg   hp -0.86695324 -0.53188409
#> 4   mpg drat  0.39175113  0.81501237
#> 5   mpg   wt -0.94349203 -0.77828961
#> 6   mpg qsec  0.14123780  0.70142529
#> 7   mpg   vs  0.41036301  0.82232624
#> 8   mpg   am  0.31755830  0.78445202
#> 9   mpg gear  0.09948247  0.67925253
#> 10  mpg carb -0.72548984 -0.18892934
#> 11  cyl disp  0.65045093  0.90586708
#> 12  cyl   hp  0.60105535  0.89016294
#> 13  cyl drat -0.75489071 -0.25084701
#> 14  cyl   wt  0.50876240  0.85885677
#> 15  cyl qsec -0.68969061 -0.11888816
#> 16  cyl   vs -0.90393935 -0.64426886
#> 17  cyl   am -0.73699794 -0.21266745
#> 18  cyl gear -0.73065432 -0.19950526
#> 19  cyl carb  0.13935172  0.70044660
#> 20 disp   hp  0.41328088  0.82346065
#> 21 disp drat -0.72204626 -0.18194600
#> 22 disp   wt  0.53279607  0.86726900
#> 23 disp qsec -0.58784853  0.05348908
#> 24 disp   vs -0.84883771 -0.48083271
#> 25 disp   am -0.77926901 -0.30551779
#> 26 disp gear -0.70728464 -0.15261532
#> 27 disp carb  0.07600741  0.66630411
#> 28   hp drat -0.64523394 -0.03915727
#> 29   hp   wt  0.33375774  0.79132390
#> 30   hp qsec -0.70529697 -0.14873917
#> 31   hp   vs -0.85596751 -0.50063181
#> 32   hp   am -0.54562696  0.11526460
#> 33   hp gear -0.57236769  0.07672422
#> 34   hp carb  0.31216260  0.78213740
#> 35 drat   wt -0.87114384 -0.54405098
#> 36 drat qsec -0.26532476  0.42688767
#> 37 drat   vs  0.10819483  0.68396800
#> 38 drat   am  0.48439908  0.85013194
#> 39 drat gear  0.29537173  0.77485059
#> 40 drat carb -0.42975706  0.26205501
#> 41   wt qsec -0.53226152  0.13380968
#> 42   wt   vs -0.75711174 -0.25569823
#> 43   wt   am -0.83867523 -0.45324615
#> 44   wt gear -0.75010777 -0.24048527
#> 45   wt carb  0.02605405  0.63751279
#> 46 qsec   vs  0.53464277  0.86790757
#> 47 qsec   am -0.53562398  0.12918764
#> 48 qsec gear -0.42638618  0.26589459
#> 49 qsec carb -0.72678816 -0.19157641
#> 50   vs   am -0.19159569  0.48837125
#> 51   vs gear -0.15371324  0.51753786
#> 52   vs carb -0.76613289 -0.27566545
#> 53   am gear  0.61589632  0.89495457
#> 54   am carb -0.29712041  0.39823890
#> 55 gear carb -0.25954979  0.43194613
#> 
#> $sig_threshold
#> [1] 0.05
#> 
#> $method
#> [1] "auto"
#> 
#> $call
#> run_correl(x = df, remove = list("cyl" ~ c("4", "6")))
#> 
#> attr(,"class")
#> [1] "pondeR_correl" "list"         

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
