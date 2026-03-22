# Compute WHO Anthropometric Z-Scores and Nutritional Status

Computes WHO Child Growth Standards z-scores for up to eight
anthropometric indices—weight-for-age (WFA), length/height-for-age
(HFA), weight-for-length/height (WFH), BMI-for-age (BMIFA), head
circumference-for-age (HCFA), arm circumference-for-age (ACFA), triceps
skinfold-for-age (TSFA), and subscapular skinfold-for-age (SSFA)—using
[`anthro_zscores`](https://rdrr.io/pkg/anthro/man/anthro_zscores.html)
(WHO Child Growth Standards 2006) for ages 0 to \< 5 years, and
[`anthroplus_zscores`](https://rdrr.io/pkg/anthroplus/man/anthroplus_zscores.html)
(WHO Growth Reference 2007) for ages 5 to \< 20 years where applicable.

For BMIFA, the function covers the full age range:

- **0 to \< 5 years**: z-scores from
  [`anthro_zscores`](https://rdrr.io/pkg/anthro/man/anthro_zscores.html),
  WHO Child Growth Standards 2006.

- **5 to \< 20 years**: z-scores from
  [`anthroplus_zscores`](https://rdrr.io/pkg/anthroplus/man/anthroplus_zscores.html),
  WHO Growth Reference 2007. The anthroplus package must be installed
  separately when participants in this age range are present.

- **20 years and above**: nutritional status classified directly from
  the computed BMI value using WHO adult cut-offs (no z-score).

A `Nutritional Status` column is appended to each result according to
WHO classification thresholds.

## Usage

``` r
run_anthroindex(
  x,
  sex_col,
  male_code = NULL,
  female_code = NULL,
  age_col,
  age_unit = c("month", "day", "week", "year"),
  weight_col = NULL,
  weight_unit = c("kg", "g", "oz", "lb"),
  height_col = NULL,
  height_unit = c("cm", "m", "ft", "in"),
  height_measured = "L",
  headcircumference_col = NULL,
  headcircumference_unit = c("cm", "m", "ft", "in"),
  armcircumference_col = NULL,
  armcircumference_unit = c("cm", "m", "ft", "in"),
  tricepskinfold_col = NULL,
  tricepskinfold_unit = c("mm", "cm", "m", "ft", "in"),
  subscapularskinfold_col = NULL,
  subscapularskinfold_unit = c("mm", "cm", "m", "ft", "in"),
  oedema_col = NULL,
  return = c("wfa", "hfa", "wfh", "bmifa", "hcfa", "acfa", "tsfa", "ssfa")
)
```

## Arguments

- x:

  A `data.frame`, `tibble`, or named `matrix` containing the
  anthropometric measurements. Column names may include special
  characters.

- sex_col:

  A single character string giving the column name of the sex variable
  in `x`. The column may be character, factor, or numeric. See
  **Details** for automatic and manual coding rules.

- male_code:

  Optional single character string indicating how male participants are
  coded in `sex_col`. Required when automatic detection fails.

- female_code:

  Optional single character string indicating how female participants
  are coded in `sex_col`. Required when automatic detection fails.

- age_col:

  A single character string giving the column name of the age variable
  in `x`.

- age_unit:

  A single character string specifying the unit of `age_col`. One of
  `"day"`, `"week"`, `"month"`, or `"year"`.

- weight_col:

  A single character string giving the column name of the weight
  variable in `x`.

- weight_unit:

  A single character string specifying the weight unit. One of `"kg"`
  (default), `"g"`, `"oz"`, or `"lb"`.

- height_col:

  A single character string giving the column name of the length/height
  variable in `x`.

- height_unit:

  A single character string specifying the height unit. One of `"cm"`
  (default), `"m"`, `"ft"`, or `"in"`.

- height_measured:

  A single valid measurement-mode code (`"L"`, `"l"`, `"Length"`,
  `"length"` for recumbent; `"H"`, `"h"`, `"Height"`, `"height"` for
  standing), *or* a character vector of length `nrow(x)` giving the mode
  for each row, *or* the name of a column in `x` containing those codes.

- headcircumference_col:

  Optional single character string giving the column name of the head
  circumference variable in `x`.

- headcircumference_unit:

  A single character string specifying the head circumference unit. One
  of `"cm"` (default), `"m"`, `"ft"`, or `"in"`.

- armcircumference_col:

  Optional single character string giving the column name of the
  mid-upper arm circumference (MUAC) variable in `x`.

- armcircumference_unit:

  A single character string specifying the MUAC unit. One of `"cm"`
  (default), `"m"`, `"ft"`, or `"in"`.

- tricepskinfold_col:

  Optional single character string giving the column name of the triceps
  skinfold variable in `x`.

- tricepskinfold_unit:

  A single character string specifying the triceps skinfold unit. One of
  `"mm"` (default), `"cm"`, `"m"`, `"ft"`, or `"in"`.

- subscapularskinfold_col:

  Optional single character string giving the column name of the
  subscapular skinfold variable in `x`.

- subscapularskinfold_unit:

  A single character string specifying the subscapular skinfold unit.
  One of `"mm"` (default), `"cm"`, `"m"`, `"ft"`, or `"in"`.

- oedema_col:

  Optional single character string giving the column name of the
  bilateral pitting oedema variable in `x`. Valid cell values are `"n"`,
  `"N"`, or `"2"` (no oedema) and `"y"`, `"Y"`, or `"1"` (oedema
  present).

- return:

  A character string or vector specifying which anthropometric indices
  to compute. One or more of: `"wfa"`, `"hfa"`, `"wfh"`, `"bmifa"`,
  `"hcfa"`, `"acfa"`, `"tsfa"`, `"ssfa"`. Defaults to all eight.

## Value

Each result is an object of class `"run_anthroindex"` and
`"data.frame"`.

The columns are ordered as follows:

1.  All original columns from `x` (returned as-is, in their original
    order).

2.  Supporting columns from `anthro_zscores` (e.g., `age_in_days`,
    `age_in_months`, `clenhei`, `cbmi`).

3.  `z_score` — the primary z-score for the requested index. For BMIFA:
    sourced from anthro (`zbmi`) for ages 0–\< 5 years and from
    anthroplus (`zbfa`) for ages 5–\< 20 years. For ages \\\geq\\ 20
    years, `z_score` is `NA` (classification uses raw BMI instead).

4.  Additional flag and auxiliary columns for the index (e.g., `fbmi`,
    `zbmia`, `sbmi`).

5.  `Nutritional Status` — a character column with the WHO
    classification label based on z-score (ages \< 20 years) or raw BMI
    (ages \\\geq\\ 20 years, BMIFA only).

If a single index is requested via `return`, the object is returned
directly. If multiple indices are requested, a named `list` of such
objects is returned, with element names matching the index codes (e.g.,
`"wfa"`, `"bmifa"`).

## Details

### Sex coding

Internally, sex is encoded as `1` = male and `2` = female, as required
by both
[`anthro_zscores`](https://rdrr.io/pkg/anthro/man/anthro_zscores.html)
and
[`anthroplus_zscores`](https://rdrr.io/pkg/anthroplus/man/anthroplus_zscores.html).
If `male_code` and `female_code` are not supplied, the following values
are detected automatically: `m`, `M`, `male`, `Male` \\\to\\ male; `f`,
`F`, `female`, `Female` \\\to\\ female. All values are coerced to
character before matching. If the column uses a different coding scheme,
supply `male_code` and `female_code` explicitly.

### Age conversion

All age units are converted to months before use. When
`age_unit = "day"`, the raw values are passed to `anthro_zscores()` with
`is_age_in_month = FALSE`; all other units are converted to months and
`is_age_in_month = TRUE` is used.
[`anthroplus::anthroplus_zscores()`](https://rdrr.io/pkg/anthroplus/man/anthroplus_zscores.html)
always receives age in months, regardless of the original input unit.

### Measurement units

All measurements are converted to the units expected internally before
analysis: weight \\\to\\ kg, lengths/circumferences \\\to\\ cm,
skinfolds \\\to\\ mm.

### Age coverage and z-score sources

The two WHO packages cover different age windows. The table below shows
which package supplies the z-score for each index and age group:

|                             |                 |                     |                   |
|-----------------------------|-----------------|---------------------|-------------------|
| Index                       | 0 to \< 5 years | 5 to \< 20 years    | \\\geq\\ 20 years |
| BMIFA                       | anthro (`zbmi`) | anthroplus (`zbfa`) | raw BMI           |
| WFA                         | anthro (`zwei`) | anthroplus (`zwfa`) | —                 |
| HFA                         | anthro (`zlen`) | anthroplus (`zhfa`) | —                 |
| WFH, HCFA, ACFA, TSFA, SSFA | anthro          | not available       | —                 |

All z-scores are stored in a unified `z_score` column in the output.

### Partial results

When a requested index cannot be computed because its required
measurement column was not supplied, a warning is issued and that index
is silently skipped. Required columns per index:

- `wfa` — `weight_col`

- `hfa` — `height_col`

- `wfh` — `weight_col` + `height_col`

- `bmifa` — `weight_col` + `height_col`

- `hcfa` — `headcircumference_col`

- `acfa` — `armcircumference_col`

- `tsfa` — `tricepskinfold_col`

- `ssfa` — `subscapularskinfold_col`

### WHO nutritional status classifications

#### Birth to 5 years

|                             |               |                                           |
|-----------------------------|---------------|-------------------------------------------|
| Nutritional status          | Index         | Threshold                                 |
| Obese                       | WFH or BMIFA  | \> +3 SD                                  |
| Overweight                  | WFH or BMIFA  | \> +2 SD and \\\leq\\ +3 SD               |
| Normal                      | WFH or BMIFA  | \\\geq\\ \\-\\2 SD and \\\leq\\ +2 SD     |
| Moderately wasted           | WFH           | \\\leq\\ \\-\\2 SD and \\\geq\\ \\-\\3 SD |
| Severely wasted             | WFH           | \< \\-\\3 SD                              |
| Moderate acute malnutrition | BMIFA or ACFA | \\\leq\\ \\-\\2 SD and \\\geq\\ \\-\\3 SD |
| Severe acute malnutrition   | BMIFA or ACFA | \< \\-\\3 SD                              |
| Moderately underweight      | WFA           | \< \\-\\2 SD and \\\geq\\ \\-\\3 SD       |
| Severely underweight        | WFA           | \< \\-\\3 SD                              |
| Normal                      | WFA           | \\\geq\\ \\-\\2 SD                        |
| Moderately stunted          | HFA           | \\\leq\\ \\-\\2 SD and \\\geq\\ \\-\\3 SD |
| Severely stunted            | HFA           | \< \\-\\3 SD                              |
| Normal                      | HFA           | \> \\-\\2 SD                              |
| Normal                      | ACFA          | \> \\-\\2 SD                              |

#### 5 to 19 years — BMI-for-age only

|                    |                                       |
|--------------------|---------------------------------------|
| Nutritional status | Z-score threshold                     |
| Obesity            | \> +2 SD                              |
| Overweight         | \> +1 SD and \\\leq\\ +2 SD           |
| Normal             | \\\geq\\ \\-\\2 SD and \\\leq\\ +1 SD |
| Thinness           | \< \\-\\2 SD and \\\geq\\ \\-\\3 SD   |
| Severe thinness    | \< \\-\\3 SD                          |

WHO nutritional status classifications for ages 5–19 years are currently
only available for BMIFA. For all other indices in participants aged 5
years and above, `Nutritional Status` will be `NA` with an informational
message.

#### 20 years and above — BMI (absolute value) only

For participants aged 20 years and above, nutritional status is derived
from the computed BMI value directly (not from a z-score), and is only
populated when `return` includes `"bmifa"`.

|                              |                                      |
|------------------------------|--------------------------------------|
| Nutritional status           | BMI threshold                        |
| Moderate and severe thinness | \< 17.0 kg/m\\^2\\                   |
| Underweight                  | \\\geq\\ 17.0 and \< 18.5 kg/m\\^2\\ |
| Normal weight                | \\\geq\\ 18.5 and \< 25.0 kg/m\\^2\\ |
| Overweight                   | \\\geq\\ 25.0 and \< 30.0 kg/m\\^2\\ |
| Stage 1 Obesity              | \\\geq\\ 30.0 and \< 35.0 kg/m\\^2\\ |
| Stage 2 Obesity              | \\\geq\\ 35.0 and \< 40.0 kg/m\\^2\\ |
| Stage 3 Obesity              | \\\geq\\ 40.0 kg/m\\^2\\             |

## References

WHO Multicentre Growth Reference Study Group (2006). *WHO Child Growth
Standards: Length/height-for-age, weight-for-age, weight-for-length,
weight-for-height and body mass index-for-age: Methods and development*.
Geneva: World Health Organization; pp 312.
<http://www.who.int/childgrowth/publications/en/>

WHO Multicentre Growth Reference Study Group (2007). *WHO Child Growth
Standards: Head circumference-for-age, arm circumference-for-age,
triceps skinfold-for-age and subscapular skinfold-for-age: Methods and
development*. Geneva: World Health Organization; pp 217.
<http://www.who.int/childgrowth/publications/en/>

de Onis M, Onyango AW, Borghi E, Siyam A, Nishida C, Siekmann J (2007).
Development of a WHO growth reference for school-aged children and
adolescents. *Bulletin of the World Health Organization*, 85, 660-667.
[doi:10.2471/BLT.07.043497](https://doi.org/10.2471/BLT.07.043497)

World Health Organization (2017). *Guideline: Assessing and Managing
Children at Primary Health-Care Facilities to Prevent Overweight and
Obesity in the Context of the Double Burden of Malnutrition: Updates for
the Integrated Management of Childhood Illness (IMCI)*. Table 1: WHO
classification of nutritional status of infants and children. Geneva:
WHO. <https://www.ncbi.nlm.nih.gov/books/NBK487900/table/fm.s1.t1/>

World Health Organization (2006). *WHO child growth standards: methods
and development. Length/height-for-age, weight-for-age,
weight-for-length, weight-for-height and body mass index-for-age*.
Geneva: WHO.
<http://www.who.int/nutrition/publications/childgrowthstandards_technical_report_1/en/>

Note: Weight-for-length is used for infants and young children aged 0-23
months; weight-for-height is used for children aged 24 months and older.

World Health Organization. Nutrition Landscape Information System
(NLiS): Malnutrition in Women — Moderate and severe thinness,
underweight, overweight, obesity. Geneva: WHO.
<https://apps.who.int/nutrition/landscape/help.aspx?menu=0&helpid=420>

Centers for Disease Control and Prevention (CDC). About Adult BMI: BMI
categories. Atlanta: CDC.
<https://www.cdc.gov/bmi/adult-calculator/bmi-categories.html>

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
library(pondeR)

if (requireNamespace("anthro", quietly = TRUE)) {

  # --- Synthetic dataset (no external data file required) ----------------
  set.seed(42)
  n <- 50
  children <- data.frame(
    child_id = seq_len(n),
    sex      = sample(c("Male", "Female"), n, replace = TRUE),
    age_mo   = runif(n, 1, 60),
    wt_kg    = runif(n, 3, 20),
    ht_cm    = runif(n, 50, 115),
    oedema   = sample(c("n", "y"), n, replace = TRUE, prob = c(0.95, 0.05))
  )

  # --- Multiple indices — returns a named list ---------------------------
  res <- run_anthroindex(
    x               = children,
    sex_col         = "sex",
    age_col         = "age_mo",
    age_unit        = "month",
    weight_col      = "wt_kg",
    weight_unit     = "kg",
    height_col      = "ht_cm",
    height_unit     = "cm",
    height_measured = "L",
    oedema_col      = "oedema",
    return          = c("wfa", "hfa", "wfh", "bmifa")
  )

  head(res$wfa[,  c("z_score", "Nutritional Status")])
  head(res$bmifa[, c("cbmi", "z_score", "Nutritional Status")])

  # --- Single index — returns a plain data frame -------------------------
  wfa_only <- run_anthroindex(
    x               = children,
    sex_col         = "sex",
    age_col         = "age_mo",
    age_unit        = "month",
    weight_col      = "wt_kg",
    weight_unit     = "kg",
    height_col      = "ht_cm",
    height_unit     = "cm",
    height_measured = "L",
    return          = "wfa"
  )
  head(wfa_only)

  # --- Mixed-age dataset spanning all three age brackets ----------------
  if (requireNamespace("anthroplus", quietly = TRUE)) {
    mixed <- data.frame(
      sex    = c("Male", "Female", "Male"),
      age_yr = c(3, 14, 25),
      wt_kg  = c(14, 38.8, 80),
      ht_cm  = c(95, 151.13, 175)
    )
    bmi_res <- run_anthroindex(
      x               = mixed,
      sex_col         = "sex",
      age_col         = "age_yr",
      age_unit        = "year",
      weight_col      = "wt_kg",
      height_col      = "ht_cm",
      height_measured = "H",
      return          = "bmifa"
    )
    bmi_res[, c("age_yr", "cbmi", "z_score", "Nutritional Status")]
  }
}
} # }
```
