# ============================================================
#  run_anthroindex — WHO Anthropometric Z-score Calculator
#  Part of the pondeR package
# ============================================================

# ---- Internal helpers (not exported) -----------------------

#' @keywords internal
.anth_convert_weight <- function(w, unit) {
  switch(unit,
    kg  = w,
    g   = w / 1000,
    oz  = w * 0.0283495,
    lb  = w * 0.453592,
    stop("'weight_unit' must be one of: 'kg', 'g', 'oz', 'lb'.", call. = FALSE)
  )
}

#' @keywords internal
.anth_convert_length <- function(l, unit) {
  switch(unit,
    cm   = l,
    m    = l * 100,
    ft   = l * 30.48,
    `in` = l * 2.54,
    stop("Height/circumference unit must be one of: 'cm', 'm', 'ft', 'in'.",
         call. = FALSE)
  )
}

#' @keywords internal
.anth_convert_skinfold <- function(s, unit) {
  switch(unit,
    mm   = s,
    cm   = s * 10,
    m    = s * 1000,
    ft   = s * 304.8,
    `in` = s * 25.4,
    stop("Skinfold unit must be one of: 'mm', 'cm', 'm', 'ft', 'in'.",
         call. = FALSE)
  )
}

#' @keywords internal
.anth_convert_age_to_months <- function(age, unit) {
  switch(unit,
    day   = age,
    week  = age * 7 / 30.4375,
    month = age,
    year  = age * 12,
    stop("'age_unit' must be one of: 'day', 'week', 'month', 'year'.",
         call. = FALSE)
  )
}

#' @keywords internal
.anth_resolve_sex <- function(col_vals, male_code, female_code) {
  # Both anthro::anthro_zscores() and anthroplus::anthroplus_zscores() require
  # numeric sex coded as 1 = male, 2 = female. This helper always returns
  # an integer vector using that convention.
  u <- unique(as.character(col_vals[!is.na(col_vals)]))

  # Default recognisable male/female labels (character and legacy numeric 0/1)
  male_defaults   <- c("m", "M", "male", "Male", "MALE", "0")
  female_defaults <- c("f", "F", "female", "Female", "FEMALE", "1")

  if (!is.null(male_code) && !is.null(female_code)) {
    mc <- as.character(male_code)
    fc <- as.character(female_code)
    if (!(mc %in% as.character(u)))
      warning("'male_code' (\"", mc, "\") was not found in 'sex_col'.",
              call. = FALSE)
    if (!(fc %in% as.character(u)))
      warning("'female_code' (\"", fc, "\") was not found in 'sex_col'.",
              call. = FALSE)
    # Encode as 1 = male, 2 = female (required by anthro and anthroplus)
    out <- ifelse(as.character(col_vals) == mc, 1L,
           ifelse(as.character(col_vals) == fc, 2L, NA_integer_))
    return(out)
  }

  str_vals <- as.character(col_vals)
  if (all(str_vals[!is.na(str_vals)] %in% c(male_defaults, female_defaults))) {
    # Encode as 1 = male, 2 = female (required by anthro and anthroplus)
    out <- ifelse(str_vals %in% male_defaults, 1L,
           ifelse(str_vals %in% female_defaults, 2L, NA_integer_))
    return(out)
  }

  stop(
    "Automatic sex detection failed. The values in 'sex_col' (",
    paste(u, collapse = ", "),
    ") are not among the recognised defaults ",
    "(m/M/male/Male for male; f/F/female/Female for female). ",
    "Please supply 'male_code' and 'female_code' explicitly.",
    call. = FALSE
  )
}

#' @keywords internal
.anth_resolve_measure <- function(x, height_measured, n) {
  valid_L <- c("L", "l", "Length", "length")
  valid_H <- c("H", "h", "Height", "height")
  valid   <- c(valid_L, valid_H)

  if (length(height_measured) == 1L && height_measured %in% names(x)) {
    col <- as.character(x[[height_measured]])
    bad <- !col %in% valid
    if (any(bad))
      stop("Column '", height_measured, "' contains invalid values: ",
           paste(unique(col[bad]), collapse = ", "), ".", call. = FALSE)
    return(col)
  }

  if (length(height_measured) == 1L) {
    hm <- as.character(height_measured)
    if (!hm %in% valid)
      stop("'height_measured' must be one of: ",
           paste(valid, collapse = ", "), ".", call. = FALSE)
    return(rep(hm, n))
  }

  if (length(height_measured) == n) {
    bad <- !height_measured %in% valid
    if (any(bad))
      stop("'height_measured' vector contains invalid values: ",
           paste(unique(height_measured[bad]), collapse = ", "),
           ". Accepted: ", paste(valid, collapse = ", "), ".", call. = FALSE)
    return(as.character(height_measured))
  }

  stop(
    "'height_measured' must be a single valid code, a column name in 'x', ",
    "or a character vector of the same length as nrow(x).",
    call. = FALSE
  )
}

#' @keywords internal
.anth_classify_0_5 <- function(zscore, index) {
  # Label assignment rationale:
  #   wfh   — wasting terminology ("wasted") is specific to weight-for-length/height
  #   bmifa — acute malnutrition terminology used when wfh is unavailable or
  #            BMI-for-age is the primary index; mirrors wfh cut-offs
  #   wfa   — underweight terminology
  #   hfa   — stunting terminology
  #   acfa  — acute malnutrition via MUAC
  if (all(is.na(zscore))) return(rep(NA_character_, length(zscore)))

  if (index == "wfa") {
    return(ifelse(zscore < -3,  "Severely underweight",
           ifelse(zscore < -2,  "Moderately underweight",
                                "Normal")))
  }
  if (index == "wfh") {
    return(ifelse(zscore < -3,  "Severely wasted",
           ifelse(zscore <= -2, "Moderately wasted",
           ifelse(zscore > 3,   "Obese",
           ifelse(zscore > 2,   "Overweight",
                                "Normal")))))
  }
  if (index == "bmifa") {
    return(ifelse(zscore < -3,  "Severe acute malnutrition",
           ifelse(zscore <= -2, "Moderate acute malnutrition",
           ifelse(zscore > 3,   "Obese",
           ifelse(zscore > 2,   "Overweight",
                                "Normal")))))
  }
  if (index == "hfa") {
    return(ifelse(zscore < -3,  "Severely stunted",
           ifelse(zscore <= -2, "Moderately stunted",
                                "Normal")))
  }
  if (index == "acfa") {
    return(ifelse(zscore < -3,  "Severe acute malnutrition",
           ifelse(zscore <= -2, "Moderate acute malnutrition",
                                "Normal")))
  }
  rep(NA_character_, length(zscore))
}

#' @keywords internal
.anth_classify_5_19_bmi <- function(zscore) {
  ifelse(is.na(zscore),  NA_character_,
  ifelse(zscore < -3,    "Severe thinness",
  ifelse(zscore < -2,    "Thinness",
  ifelse(zscore <= 1,    "Normal",
  ifelse(zscore <= 2,    "Overweight",
                         "Obesity")))))
}

#' @keywords internal
.anth_classify_adult_bmi <- function(bmi) {
  ifelse(is.na(bmi),  NA_character_,
  ifelse(bmi < 17.0,  "Moderate and severe thinness",
  ifelse(bmi < 18.5,  "Underweight",
  ifelse(bmi < 25.0,  "Normal weight",
  ifelse(bmi < 30.0,  "Overweight",
  ifelse(bmi < 35.0,  "Stage 1 Obesity",
  ifelse(bmi < 40.0,  "Stage 2 Obesity",
                      "Stage 3 Obesity")))))))
}

#' @keywords internal
.anth_zscore_5_19 <- function(sex_int, age_months, weight_kg, height_cm) {
  # Computes z-scores for ages 5-19 years (60-228 months) using the WHO
  # Growth Reference 2007 via anthroplus::anthroplus_zscores().
  #
  # anthroplus::anthroplus_zscores() argument conventions:
  #   sex           : 1 = male, 2 = female  (integer, same as anthro)
  #   age_in_months : age in months          (must be in months)
  #   weight_in_kg  : weight in kilograms
  #   height_in_cm  : standing height in centimetres
  #
  # Returned z-score columns (parallel to anthro equivalents):
  #   zhfa  <->  zlen  (length/height-for-age)
  #   zwfa  <->  zwei  (weight-for-age)
  #   zbfa  <->  zbmi  (BMI-for-age)
  #
  # Returns a named list of numeric vectors (one per z-score), each of length
  # equal to the full input length; rows outside 60-228 months or with missing
  # inputs carry NA.

  if (!requireNamespace("anthroplus", quietly = TRUE))
    stop(
      "Package 'anthroplus' is required for participants aged 5-19 years but ",
      "is not installed.\nInstall it with: install.packages(\"anthroplus\")",
      call. = FALSE
    )

  n     <- length(age_months)
  # Initialise all three z-score vectors to NA
  z_out <- list(
    zbfa = rep(NA_real_, n),
    zhfa = rep(NA_real_, n),
    zwfa = rep(NA_real_, n)
  )

  valid <- !is.na(sex_int) & !is.na(age_months) &
           !is.na(weight_kg) & !is.na(height_cm) &
           age_months >= 60 & age_months <= 228

  if (!any(valid)) return(z_out)

  res <- anthroplus::anthroplus_zscores(
    sex           = sex_int[valid],
    age_in_months = age_months[valid],   # anthroplus requires months
    weight_in_kg  = weight_kg[valid],
    height_in_cm  = height_cm[valid]
  )

  # Map anthroplus column names to our output slots
  col_slots <- c(zbfa = "zbfa", zhfa = "zhfa", zwfa = "zwfa")
  for (slot in names(col_slots)) {
    src <- col_slots[[slot]]
    if (src %in% names(res)) {
      z_out[[slot]][valid] <- res[[src]]
    } else {
      warning(
        "anthroplus::anthroplus_zscores() did not return column '", src, "'. ",
        "Columns returned: ", paste(names(res), collapse = ", "),
        ".\nz_score will be NA for ages 5-19 for this index. ",
        "Check your anthroplus package version.",
        call. = FALSE
      )
    }
  }

  z_out
}

#' @keywords internal
.anth_add_status <- function(df, index, age_months, sex_int,
                              weight_kg, height_cm) {
  # Map each index to:
  #   zcol_anthro  — column from anthro_zscores()    (valid for 0-<5 yr)
  #   zcol_plus    — column from anthroplus_zscores() (valid for 5-<20 yr)
  #                  NULL means anthroplus does not cover this index
  index_zcols <- list(
    wfa   = list(anthro = "zwei", plus = "zwfa"),
    hfa   = list(anthro = "zlen", plus = "zhfa"),
    wfh   = list(anthro = "zwfl", plus = NULL),
    bmifa = list(anthro = "zbmi", plus = "zbfa"),
    hcfa  = list(anthro = "zhc",  plus = NULL),
    acfa  = list(anthro = "zac",  plus = NULL),
    tsfa  = list(anthro = "zts",  plus = NULL),
    ssfa  = list(anthro = "zss",  plus = NULL)
  )

  zcols <- index_zcols[[index]]
  zcol  <- zcols$anthro

  if (is.null(zcol) || !zcol %in% names(df)) {
    df[["z_score"]]           <- NA_real_
    df[["Nutritional Status"]] <- NA_character_
    return(df)
  }

  # Rename the anthro z-score column to "z_score"; it holds valid values for
  # 0-<5 yr and NA beyond — the 5-19 yr rows are filled in below.
  names(df)[names(df) == zcol] <- "z_score"
  z      <- df[["z_score"]]
  status <- rep(NA_character_, length(z))

  # Age strata (in months)
  i_u5  <- !is.na(age_months) & age_months <  60
  i_519 <- !is.na(age_months) & age_months >= 60  & age_months < 240
  i_a20 <- !is.na(age_months) & age_months >= 240

  # ---- 0 to <5 years: classify from anthro z-scores -----------------------
  if (index %in% c("wfa", "hfa", "wfh", "bmifa", "acfa")) {
    if (any(i_u5))
      status[i_u5] <- .anth_classify_0_5(z[i_u5], index)
  }

  # ---- 5 to <20 years: fill z-scores from anthroplus, then classify --------
  if (any(i_519) && !is.null(zcols$plus)) {
    z_plus <- .anth_zscore_5_19(
      sex_int    = sex_int[i_519],
      age_months = age_months[i_519],   # already in months
      weight_kg  = weight_kg[i_519],
      height_cm  = height_cm[i_519]
    )

    slot <- switch(index,
      wfa   = "zwfa",
      hfa   = "zhfa",
      bmifa = "zbfa",
      NULL
    )

    if (!is.null(slot) && slot %in% names(z_plus)) {
      z[i_519]       <- z_plus[[slot]]
      df[["z_score"]] <- z

      if (index == "bmifa") {
        status[i_519] <- .anth_classify_5_19_bmi(z_plus[["zbfa"]])
      } else if (index %in% c("wfa", "hfa")) {
        status[i_519] <- .anth_classify_0_5(z_plus[[slot]], index)
      }
    }
  } else if (any(i_519) && is.null(zcols$plus)) {
    # Indices not covered by anthroplus above 5 yr — leave status NA
    if (!index %in% c("hcfa", "tsfa", "ssfa"))
      message(
        "No WHO z-score reference is available for index '", index,
        "' for ages 5-19 years. 'Nutritional Status' will be NA for those rows."
      )
  }

  # ---- 20+ years: BMI-for-age only; classify from raw cbmi -----------------
  if (index == "bmifa" && any(i_a20)) {
    if ("cbmi" %in% names(df)) {
      status[i_a20] <- .anth_classify_adult_bmi(df[["cbmi"]][i_a20])
    } else {
      warning(
        "Adult BMI classification (age \u2265 20 years) requires the computed ",
        "BMI column 'cbmi' from anthro_zscores(). Column not found; ",
        "'Nutritional Status' will be NA for these rows.",
        call. = FALSE
      )
    }
  }

  # hcfa, tsfa, ssfa: no WHO classification at any age
  if (index %in% c("hcfa", "tsfa", "ssfa") && !all(is.na(z)))
    message(
      "No WHO nutritional status classification is currently available for '",
      index, "'. 'Nutritional Status' will be NA."
    )

  df[["Nutritional Status"]] <- status
  df
}


# ---- Main exported function ------------------------------------------------

#' @title Compute WHO Anthropometric Z-Scores and Nutritional Status
#'
#' @description
#' Computes WHO Child Growth Standards z-scores for up to eight anthropometric
#' indices—weight-for-age (WFA), length/height-for-age (HFA),
#' weight-for-length/height (WFH), BMI-for-age (BMIFA), head
#' circumference-for-age (HCFA), arm circumference-for-age (ACFA), triceps
#' skinfold-for-age (TSFA), and subscapular skinfold-for-age (SSFA)—using
#' \code{\link[anthro]{anthro_zscores}} (WHO Child Growth Standards 2006) for
#' ages 0 to < 5 years, and
#' \code{\link[anthroplus]{anthroplus_zscores}} (WHO Growth Reference 2007)
#' for ages 5 to < 20 years where applicable.
#'
#' For BMIFA, the function covers the full age range:
#' \itemize{
#'   \item \strong{0 to < 5 years}: z-scores from
#'     \code{\link[anthro]{anthro_zscores}}, WHO Child Growth Standards 2006.
#'   \item \strong{5 to < 20 years}: z-scores from
#'     \code{\link[anthroplus]{anthroplus_zscores}}, WHO Growth Reference 2007.
#'     The \pkg{anthroplus} package must be installed separately when
#'     participants in this age range are present.
#'   \item \strong{20 years and above}: nutritional status classified directly
#'     from the computed BMI value using WHO adult cut-offs (no z-score).
#' }
#' A \code{Nutritional Status} column is appended to each result according to
#' WHO classification thresholds.
#'
#' @param x A \code{data.frame}, \code{tibble}, or named \code{matrix}
#'   containing the anthropometric measurements. Column names may include
#'   special characters.
#' @param sex_col A single character string giving the column name of the sex
#'   variable in \code{x}. The column may be character, factor, or numeric.
#'   See \strong{Details} for automatic and manual coding rules.
#' @param male_code Optional single character string indicating how male
#'   participants are coded in \code{sex_col}. Required when automatic
#'   detection fails.
#' @param female_code Optional single character string indicating how female
#'   participants are coded in \code{sex_col}. Required when automatic
#'   detection fails.
#' @param age_col A single character string giving the column name of the age
#'   variable in \code{x}.
#' @param age_unit A single character string specifying the unit of
#'   \code{age_col}. One of \code{"day"}, \code{"week"}, \code{"month"}, or
#'   \code{"year"}.
#' @param weight_col A single character string giving the column name of the
#'   weight variable in \code{x}.
#' @param weight_unit A single character string specifying the weight unit.
#'   One of \code{"kg"} (default), \code{"g"}, \code{"oz"}, or \code{"lb"}.
#' @param height_col A single character string giving the column name of the
#'   length/height variable in \code{x}.
#' @param height_unit A single character string specifying the height unit.
#'   One of \code{"cm"} (default), \code{"m"}, \code{"ft"}, or \code{"in"}.
#' @param height_measured A single valid measurement-mode code (\code{"L"},
#'   \code{"l"}, \code{"Length"}, \code{"length"} for recumbent;
#'   \code{"H"}, \code{"h"}, \code{"Height"}, \code{"height"} for standing),
#'   \emph{or} a character vector of length \code{nrow(x)} giving the mode
#'   for each row, \emph{or} the name of a column in \code{x} containing
#'   those codes.
#' @param headcircumference_col Optional single character string giving the
#'   column name of the head circumference variable in \code{x}.
#' @param headcircumference_unit A single character string specifying the head
#'   circumference unit. One of \code{"cm"} (default), \code{"m"},
#'   \code{"ft"}, or \code{"in"}.
#' @param armcircumference_col Optional single character string giving the
#'   column name of the mid-upper arm circumference (MUAC) variable in
#'   \code{x}.
#' @param armcircumference_unit A single character string specifying the MUAC
#'   unit. One of \code{"cm"} (default), \code{"m"}, \code{"ft"}, or
#'   \code{"in"}.
#' @param tricepskinfold_col Optional single character string giving the column
#'   name of the triceps skinfold variable in \code{x}.
#' @param tricepskinfold_unit A single character string specifying the triceps
#'   skinfold unit. One of \code{"mm"} (default), \code{"cm"}, \code{"m"},
#'   \code{"ft"}, or \code{"in"}.
#' @param subscapularskinfold_col Optional single character string giving the
#'   column name of the subscapular skinfold variable in \code{x}.
#' @param subscapularskinfold_unit A single character string specifying the
#'   subscapular skinfold unit. One of \code{"mm"} (default), \code{"cm"},
#'   \code{"m"}, \code{"ft"}, or \code{"in"}.
#' @param oedema_col Optional single character string giving the column name of
#'   the bilateral pitting oedema variable in \code{x}. Valid cell values are
#'   \code{"n"}, \code{"N"}, or \code{"2"} (no oedema) and \code{"y"},
#'   \code{"Y"}, or \code{"1"} (oedema present).
#' @param return A character string or vector specifying which anthropometric
#'   indices to compute. One or more of: \code{"wfa"}, \code{"hfa"},
#'   \code{"wfh"}, \code{"bmifa"}, \code{"hcfa"}, \code{"acfa"},
#'   \code{"tsfa"}, \code{"ssfa"}. Defaults to all eight.
#'
#' @details
#' ## Sex coding
#' Internally, sex is encoded as \code{1} = male and \code{2} = female, as
#' required by both \code{\link[anthro]{anthro_zscores}} and
#' \code{\link[anthroplus]{anthroplus_zscores}}. If \code{male_code} and
#' \code{female_code} are not supplied, the following values are detected
#' automatically: \code{m}, \code{M}, \code{male}, \code{Male} \eqn{\to}
#' male; \code{f}, \code{F}, \code{female}, \code{Female} \eqn{\to} female.
#' All values are coerced to character before matching. If the column uses a
#' different coding scheme, supply \code{male_code} and \code{female_code}
#' explicitly.
#'
#' ## Age conversion
#' All age units are converted to months before use.
#' When \code{age_unit = "day"}, the raw values are passed to
#' \code{anthro_zscores()} with \code{is_age_in_month = FALSE}; all other
#' units are converted to months and \code{is_age_in_month = TRUE} is used.
#' \code{anthroplus::anthroplus_zscores()} always receives age in months,
#' regardless of the original input unit.
#'
#' ## Measurement units
#' All measurements are converted to the units expected internally before
#' analysis: weight \eqn{\to} kg, lengths/circumferences \eqn{\to} cm,
#' skinfolds \eqn{\to} mm.
#'
#' ## Age coverage and z-score sources
#' The two WHO packages cover different age windows. The table below shows
#' which package supplies the z-score for each index and age group:
#'
#' | Index | 0 to < 5 years | 5 to < 20 years | \eqn{\geq} 20 years |
#' |---|---|---|---|
#' | BMIFA | \pkg{anthro} (\code{zbmi}) | \pkg{anthroplus} (\code{zbfa}) | raw BMI |
#' | WFA   | \pkg{anthro} (\code{zwei}) | \pkg{anthroplus} (\code{zwfa}) | — |
#' | HFA   | \pkg{anthro} (\code{zlen}) | \pkg{anthroplus} (\code{zhfa}) | — |
#' | WFH, HCFA, ACFA, TSFA, SSFA | \pkg{anthro} | not available | — |
#'
#' All z-scores are stored in a unified \code{z_score} column in the output.
#'
#' ## Partial results
#' When a requested index cannot be computed because its required measurement
#' column was not supplied, a warning is issued and that index is silently
#' skipped. Required columns per index:
#' \itemize{
#'   \item \code{wfa}   — \code{weight_col}
#'   \item \code{hfa}   — \code{height_col}
#'   \item \code{wfh}   — \code{weight_col} + \code{height_col}
#'   \item \code{bmifa} — \code{weight_col} + \code{height_col}
#'   \item \code{hcfa}  — \code{headcircumference_col}
#'   \item \code{acfa}  — \code{armcircumference_col}
#'   \item \code{tsfa}  — \code{tricepskinfold_col}
#'   \item \code{ssfa}  — \code{subscapularskinfold_col}
#' }
#'
#' ## WHO nutritional status classifications
#'
#' ### Birth to 5 years
#'
#' | Nutritional status | Index | Threshold |
#' |---|---|---|
#' | Obese | WFH or BMIFA | > +3 SD |
#' | Overweight | WFH or BMIFA | > +2 SD and \eqn{\leq} +3 SD |
#' | Normal | WFH or BMIFA | \eqn{\geq} \eqn{-}2 SD and \eqn{\leq} +2 SD |
#' | Moderately wasted | WFH | \eqn{\leq} \eqn{-}2 SD and \eqn{\geq} \eqn{-}3 SD |
#' | Severely wasted | WFH | < \eqn{-}3 SD |
#' | Moderate acute malnutrition | BMIFA or ACFA | \eqn{\leq} \eqn{-}2 SD and \eqn{\geq} \eqn{-}3 SD |
#' | Severe acute malnutrition | BMIFA or ACFA | < \eqn{-}3 SD |
#' | Moderately underweight | WFA | < \eqn{-}2 SD and \eqn{\geq} \eqn{-}3 SD |
#' | Severely underweight | WFA | < \eqn{-}3 SD |
#' | Normal | WFA | \eqn{\geq} \eqn{-}2 SD |
#' | Moderately stunted | HFA | \eqn{\leq} \eqn{-}2 SD and \eqn{\geq} \eqn{-}3 SD |
#' | Severely stunted | HFA | < \eqn{-}3 SD |
#' | Normal | HFA | > \eqn{-}2 SD |
#' | Normal | ACFA | > \eqn{-}2 SD |
#'
#' ### 5 to 19 years — BMI-for-age only
#'
#' | Nutritional status | Z-score threshold |
#' |---|---|
#' | Obesity | > +2 SD |
#' | Overweight | > +1 SD and \eqn{\leq} +2 SD |
#' | Normal | \eqn{\geq} \eqn{-}2 SD and \eqn{\leq} +1 SD |
#' | Thinness | < \eqn{-}2 SD and \eqn{\geq} \eqn{-}3 SD |
#' | Severe thinness | < \eqn{-}3 SD |
#'
#' WHO nutritional status classifications for ages 5–19 years are currently
#' only available for BMIFA. For all other indices in participants aged 5 years
#' and above, \code{Nutritional Status} will be \code{NA} with an informational
#' message.
#'
#' ### 20 years and above — BMI (absolute value) only
#'
#' For participants aged 20 years and above, nutritional status is derived from
#' the computed BMI value directly (not from a z-score), and is only populated
#' when \code{return} includes \code{"bmifa"}.
#'
#' | Nutritional status | BMI threshold |
#' |---|---|
#' | Moderate and severe thinness | < 17.0 kg/m\eqn{^2} |
#' | Underweight | \eqn{\geq} 17.0 and < 18.5 kg/m\eqn{^2} |
#' | Normal weight | \eqn{\geq} 18.5 and < 25.0 kg/m\eqn{^2} |
#' | Overweight | \eqn{\geq} 25.0 and < 30.0 kg/m\eqn{^2} |
#' | Stage 1 Obesity | \eqn{\geq} 30.0 and < 35.0 kg/m\eqn{^2} |
#' | Stage 2 Obesity | \eqn{\geq} 35.0 and < 40.0 kg/m\eqn{^2} |
#' | Stage 3 Obesity | \eqn{\geq} 40.0 kg/m\eqn{^2} |
#'
#' @return
#' Each result is an object of class \code{"run_anthroindex"} and
#' \code{"data.frame"}.
#'
#' The columns are ordered as follows:
#' \enumerate{
#'   \item All original columns from \code{x} (returned as-is, in their
#'     original order).
#'   \item Supporting columns from \code{anthro_zscores} (e.g.,
#'     \code{age_in_days}, \code{age_in_months}, \code{clenhei},
#'     \code{cbmi}).
#'   \item \code{z_score} — the primary z-score for the requested index.
#'     For BMIFA: sourced from \pkg{anthro} (\code{zbmi}) for ages 0–< 5 years
#'     and from \pkg{anthroplus} (\code{zbfa}) for ages 5–< 20 years. For ages
#'     \eqn{\geq} 20 years, \code{z_score} is \code{NA} (classification uses
#'     raw BMI instead).
#'   \item Additional flag and auxiliary columns for the index (e.g.,
#'     \code{fbmi}, \code{zbmia}, \code{sbmi}).
#'   \item \code{Nutritional Status} — a character column with the WHO
#'     classification label based on z-score (ages < 20 years) or raw BMI
#'     (ages \eqn{\geq} 20 years, BMIFA only).
#' }
#'
#' If a single index is requested via \code{return}, the object is returned
#' directly. If multiple indices are requested, a named \code{list} of such
#' objects is returned, with element names matching the index codes (e.g.,
#' \code{"wfa"}, \code{"bmifa"}).
#'
#' @examples
#' \dontrun{
#' library(pondeR)
#'
#' if (requireNamespace("anthro", quietly = TRUE)) {
#'
#'   # --- Synthetic dataset (no external data file required) ----------------
#'   set.seed(42)
#'   n <- 50
#'   children <- data.frame(
#'     child_id = seq_len(n),
#'     sex      = sample(c("Male", "Female"), n, replace = TRUE),
#'     age_mo   = runif(n, 1, 60),
#'     wt_kg    = runif(n, 3, 20),
#'     ht_cm    = runif(n, 50, 115),
#'     oedema   = sample(c("n", "y"), n, replace = TRUE, prob = c(0.95, 0.05))
#'   )
#'
#'   # --- Multiple indices — returns a named list ---------------------------
#'   res <- run_anthroindex(
#'     x               = children,
#'     sex_col         = "sex",
#'     age_col         = "age_mo",
#'     age_unit        = "month",
#'     weight_col      = "wt_kg",
#'     weight_unit     = "kg",
#'     height_col      = "ht_cm",
#'     height_unit     = "cm",
#'     height_measured = "L",
#'     oedema_col      = "oedema",
#'     return          = c("wfa", "hfa", "wfh", "bmifa")
#'   )
#'
#'   head(res$wfa[,  c("z_score", "Nutritional Status")])
#'   head(res$bmifa[, c("cbmi", "z_score", "Nutritional Status")])
#'
#'   # --- Single index — returns a plain data frame -------------------------
#'   wfa_only <- run_anthroindex(
#'     x               = children,
#'     sex_col         = "sex",
#'     age_col         = "age_mo",
#'     age_unit        = "month",
#'     weight_col      = "wt_kg",
#'     weight_unit     = "kg",
#'     height_col      = "ht_cm",
#'     height_unit     = "cm",
#'     height_measured = "L",
#'     return          = "wfa"
#'   )
#'   head(wfa_only)
#'
#'   # --- Mixed-age dataset spanning all three age brackets ----------------
#'   if (requireNamespace("anthroplus", quietly = TRUE)) {
#'     mixed <- data.frame(
#'       sex    = c("Male", "Female", "Male"),
#'       age_yr = c(3, 14, 25),
#'       wt_kg  = c(14, 38.8, 80),
#'       ht_cm  = c(95, 151.13, 175)
#'     )
#'     bmi_res <- run_anthroindex(
#'       x               = mixed,
#'       sex_col         = "sex",
#'       age_col         = "age_yr",
#'       age_unit        = "year",
#'       weight_col      = "wt_kg",
#'       height_col      = "ht_cm",
#'       height_measured = "H",
#'       return          = "bmifa"
#'     )
#'     bmi_res[, c("age_yr", "cbmi", "z_score", "Nutritional Status")]
#'   }
#' }
#' }
#'
#' @references
#' WHO Multicentre Growth Reference Study Group (2006). \emph{WHO Child Growth
#' Standards: Length/height-for-age, weight-for-age, weight-for-length,
#' weight-for-height and body mass index-for-age: Methods and development}.
#' Geneva: World Health Organization; pp 312.
#' \url{http://www.who.int/childgrowth/publications/en/}
#'
#' WHO Multicentre Growth Reference Study Group (2007). \emph{WHO Child Growth
#' Standards: Head circumference-for-age, arm circumference-for-age, triceps
#' skinfold-for-age and subscapular skinfold-for-age: Methods and development}.
#' Geneva: World Health Organization; pp 217.
#' \url{http://www.who.int/childgrowth/publications/en/}
#'
#' de Onis M, Onyango AW, Borghi E, Siyam A, Nishida C, Siekmann J (2007).
#' Development of a WHO growth reference for school-aged children and
#' adolescents. \emph{Bulletin of the World Health Organization}, 85, 660-667.
#' \doi{10.2471/BLT.07.043497}
#'
#' World Health Organization (2017). \emph{Guideline: Assessing and Managing
#' Children at Primary Health-Care Facilities to Prevent Overweight and Obesity
#' in the Context of the Double Burden of Malnutrition: Updates for the
#' Integrated Management of Childhood Illness (IMCI)}. Table 1: WHO
#' classification of nutritional status of infants and children. Geneva: WHO.
#' \url{https://www.ncbi.nlm.nih.gov/books/NBK487900/table/fm.s1.t1/}
#'
#' World Health Organization (2006). \emph{WHO child growth standards: methods
#' and development. Length/height-for-age, weight-for-age, weight-for-length,
#' weight-for-height and body mass index-for-age}. Geneva: WHO.
#' \url{http://www.who.int/nutrition/publications/childgrowthstandards_technical_report_1/en/}
#'
#' Note: Weight-for-length is used for infants and young children aged 0-23
#' months; weight-for-height is used for children aged 24 months and older.
#'
#' World Health Organization. Nutrition Landscape Information System (NLiS):
#' Malnutrition in Women — Moderate and severe thinness, underweight,
#' overweight, obesity. Geneva: WHO.
#' \url{https://apps.who.int/nutrition/landscape/help.aspx?menu=0&helpid=420}
#'
#' Centers for Disease Control and Prevention (CDC). About Adult BMI:
#' BMI categories. Atlanta: CDC.
#' \url{https://www.cdc.gov/bmi/adult-calculator/bmi-categories.html}
#'
#' @author John Lennon L. Calorio
#' @export
run_anthroindex <- function(
    x,
    sex_col,
    male_code                = NULL,
    female_code              = NULL,
    age_col,
    age_unit                 = c("month", "day", "week", "year"),
    weight_col               = NULL,
    weight_unit              = c("kg", "g", "oz", "lb"),
    height_col               = NULL,
    height_unit              = c("cm", "m", "ft", "in"),
    height_measured          = "L",
    headcircumference_col    = NULL,
    headcircumference_unit   = c("cm", "m", "ft", "in"),
    armcircumference_col     = NULL,
    armcircumference_unit    = c("cm", "m", "ft", "in"),
    tricepskinfold_col       = NULL,
    tricepskinfold_unit      = c("mm", "cm", "m", "ft", "in"),
    subscapularskinfold_col  = NULL,
    subscapularskinfold_unit = c("mm", "cm", "m", "ft", "in"),
    oedema_col               = NULL,
    return = c("wfa", "hfa", "wfh", "bmifa", "hcfa", "acfa", "tsfa", "ssfa")
) {

  # ---- 0. Package dependency check ----------------------------------------
  if (!requireNamespace("anthro", quietly = TRUE))
    stop(
      "Package 'anthro' is required but is not installed.\n",
      "Install it with: install.packages(\"anthro\")",
      call. = FALSE
    )

  # ---- 1. Validate 'x' ----------------------------------------------------
  if (is.matrix(x)) {
    if (is.null(colnames(x)))
      stop("'x' is a matrix but has no column names. Please provide named columns.",
           call. = FALSE)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else if (!is.data.frame(x)) {
    stop("'x' must be a data.frame, tibble, or named matrix.", call. = FALSE)
  }

  if (nrow(x) == 0L)
    stop("'x' has zero rows. Please supply a non-empty dataset.", call. = FALSE)

  n <- nrow(x)

  # ---- 2. Validate column-name arguments ----------------------------------
  .chk_col <- function(arg_name, arg_val, df) {
    if (!is.character(arg_val) || length(arg_val) != 1L)
      stop("'", arg_name, "' must be a single character string.", call. = FALSE)
    if (!arg_val %in% names(df))
      stop("Column '", arg_val, "' specified in '", arg_name,
           "' was not found in 'x'.\n",
           "Available columns: ", paste(names(df), collapse = ", "),
           call. = FALSE)
  }

  .chk_col("sex_col", sex_col, x)
  .chk_col("age_col", age_col, x)

  optional_cols <- list(
    weight_col              = weight_col,
    height_col              = height_col,
    headcircumference_col   = headcircumference_col,
    armcircumference_col    = armcircumference_col,
    tricepskinfold_col      = tricepskinfold_col,
    subscapularskinfold_col = subscapularskinfold_col,
    oedema_col              = oedema_col
  )
  for (nm in names(optional_cols)) {
    v <- optional_cols[[nm]]
    if (!is.null(v)) .chk_col(nm, v, x)
  }

  # ---- 3. Match unit arguments --------------------------------------------
  age_unit                 <- match.arg(age_unit)
  weight_unit              <- match.arg(weight_unit)
  height_unit              <- match.arg(height_unit)
  headcircumference_unit   <- match.arg(headcircumference_unit)
  armcircumference_unit    <- match.arg(armcircumference_unit)
  tricepskinfold_unit      <- match.arg(tricepskinfold_unit)
  subscapularskinfold_unit <- match.arg(subscapularskinfold_unit)

  # ---- 4. Validate 'return' -----------------------------------------------
  valid_indices <- c("wfa", "hfa", "wfh", "bmifa", "hcfa", "acfa", "tsfa", "ssfa")
  return        <- unique(tolower(return))
  bad_idx       <- setdiff(return, valid_indices)
  if (length(bad_idx) > 0L)
    stop(
      "Unknown index code(s) in 'return': ",
      paste(bad_idx, collapse = ", "),
      ".\nValid codes are: ", paste(valid_indices, collapse = ", "),
      call. = FALSE
    )

  # ---- 5. Resolve sex column ----------------------------------------------
  # Both anthro and anthroplus require 1 = male, 2 = female.
  sex_int <- .anth_resolve_sex(x[[sex_col]], male_code, female_code)

  if (all(is.na(sex_int)))
    stop("All values in 'sex_col' resolved to NA. Check the column contents.",
         call. = FALSE)
  if (any(is.na(sex_int)))
    warning("Some values in 'sex_col' could not be resolved and are set to NA. ",
            "Corresponding rows will have NA z-scores.", call. = FALSE)

  # ---- 6. Resolve age -----------------------------------------------------
  age_raw <- x[[age_col]]
  if (!is.numeric(age_raw))
    stop("Column '", age_col, "' must be numeric.", call. = FALSE)
  if (any(is.na(age_raw)))
    warning("'age_col' contains ", sum(is.na(age_raw)),
            " NA value(s). Corresponding rows will have NA z-scores.",
            call. = FALSE)
  if (any(age_raw < 0, na.rm = TRUE))
    stop("'age_col' contains negative values. Age must be non-negative.",
         call. = FALSE)

  is_age_in_month <- age_unit != "day"
  # age_months is always in months (or raw days if age_unit == "day");
  # anthroplus always receives months, so this vector is used directly for it.
  age_months <- .anth_convert_age_to_months(age_raw, age_unit)

  # ---- 7. Resolve measurements --------------------------------------------
  weight_kg <- if (!is.null(weight_col)) {
    w <- x[[weight_col]]
    if (!is.numeric(w))
      stop("Column '", weight_col, "' must be numeric.", call. = FALSE)
    if (any(w <= 0, na.rm = TRUE))
      warning("'weight_col' contains non-positive values. These will produce NA z-scores.",
              call. = FALSE)
    .anth_convert_weight(w, weight_unit)
  } else NULL

  height_cm <- if (!is.null(height_col)) {
    h <- x[[height_col]]
    if (!is.numeric(h))
      stop("Column '", height_col, "' must be numeric.", call. = FALSE)
    if (any(h <= 0, na.rm = TRUE))
      warning("'height_col' contains non-positive values. These will produce NA z-scores.",
              call. = FALSE)
    .anth_convert_length(h, height_unit)
  } else NULL

  measure_vec <- if (!is.null(height_col)) {
    .anth_resolve_measure(x, height_measured, n)
  } else rep("H", n)

  headc_cm <- if (!is.null(headcircumference_col)) {
    hc <- x[[headcircumference_col]]
    if (!is.numeric(hc))
      stop("Column '", headcircumference_col, "' must be numeric.", call. = FALSE)
    .anth_convert_length(hc, headcircumference_unit)
  } else NULL

  armc_cm <- if (!is.null(armcircumference_col)) {
    ac <- x[[armcircumference_col]]
    if (!is.numeric(ac))
      stop("Column '", armcircumference_col, "' must be numeric.", call. = FALSE)
    .anth_convert_length(ac, armcircumference_unit)
  } else NULL

  triskin_mm <- if (!is.null(tricepskinfold_col)) {
    ts <- x[[tricepskinfold_col]]
    if (!is.numeric(ts))
      stop("Column '", tricepskinfold_col, "' must be numeric.", call. = FALSE)
    .anth_convert_skinfold(ts, tricepskinfold_unit)
  } else NULL

  subskin_mm <- if (!is.null(subscapularskinfold_col)) {
    ss <- x[[subscapularskinfold_col]]
    if (!is.numeric(ss))
      stop("Column '", subscapularskinfold_col, "' must be numeric.", call. = FALSE)
    .anth_convert_skinfold(ss, subscapularskinfold_unit)
  } else NULL

  oedema_vec <- if (!is.null(oedema_col)) {
    oe <- as.character(x[[oedema_col]])
    valid_n <- c("n", "N", "2")
    valid_y <- c("y", "Y", "1")
    bad_oe  <- !oe %in% c(valid_n, valid_y, NA)
    if (any(bad_oe, na.rm = TRUE))
      stop(
        "'oedema_col' contains invalid values: ",
        paste(unique(oe[bad_oe]), collapse = ", "),
        ".\nAccepted: ", paste(c(valid_n, valid_y), collapse = ", "),
        call. = FALSE
      )
    ifelse(oe %in% valid_y, "y",
    ifelse(oe %in% valid_n, "n", NA_character_))
  } else rep("n", n)

  # ---- 8. Check feasibility of each requested index -----------------------
  index_req <- list(
    wfa   = list(cols = list(weight_col),              msg = "weight_col"),
    hfa   = list(cols = list(height_col),              msg = "height_col"),
    wfh   = list(cols = list(weight_col, height_col),  msg = "weight_col and height_col"),
    bmifa = list(cols = list(weight_col, height_col),  msg = "weight_col and height_col"),
    hcfa  = list(cols = list(headcircumference_col),   msg = "headcircumference_col"),
    acfa  = list(cols = list(armcircumference_col),    msg = "armcircumference_col"),
    tsfa  = list(cols = list(tricepskinfold_col),      msg = "tricepskinfold_col"),
    ssfa  = list(cols = list(subscapularskinfold_col), msg = "subscapularskinfold_col")
  )

  feasible <- vapply(return, function(idx) {
    all(vapply(index_req[[idx]]$cols, Negate(is.null), logical(1L)))
  }, logical(1L))

  infeasible <- return[!feasible]
  for (idx in infeasible)
    warning("Index '", idx, "' was requested but requires ",
            index_req[[idx]]$msg, ", which was not supplied. Skipping '", idx, "'.",
            call. = FALSE)

  return_feasible <- return[feasible]
  if (length(return_feasible) == 0L)
    stop(
      "No feasible indices remain after checking required columns. ",
      "Please supply at least one measurement column.",
      call. = FALSE
    )

  # Warn early if anthroplus will be needed but is not installed
  needs_plus <- any(c("wfa", "hfa", "bmifa") %in% return_feasible) &&
                any(!is.na(age_months) & age_months >= 60 & age_months < 240)
  if (needs_plus && !requireNamespace("anthroplus", quietly = TRUE))
    warning(
      "Participants aged 5 to < 20 years are present and at least one of ",
      "'wfa', 'hfa', 'bmifa' was requested, but package 'anthroplus' is not ",
      "installed. z_score and 'Nutritional Status' will be NA for those rows.\n",
      "Install with: install.packages(\"anthroplus\")",
      call. = FALSE
    )

  # ---- 9. Call anthro_zscores once ----------------------------------------
  # anthro_zscores() calls is.numeric() on every argument and rejects NULL.
  # Supply numeric(0)/character(0) as the "not provided" sentinel.
  zscores_all <- anthro::anthro_zscores(
    sex             = sex_int,
    age             = age_months,
    is_age_in_month = is_age_in_month,
    weight          = if (!is.null(weight_kg))  weight_kg   else numeric(0),
    lenhei          = if (!is.null(height_cm))  height_cm   else numeric(0),
    measure         = if (!is.null(height_col)) measure_vec else character(0),
    headc           = if (!is.null(headc_cm))   headc_cm    else numeric(0),
    armc            = if (!is.null(armc_cm))    armc_cm     else numeric(0),
    triskin         = if (!is.null(triskin_mm)) triskin_mm  else numeric(0),
    subskin         = if (!is.null(subskin_mm)) subskin_mm  else numeric(0),
    oedema          = oedema_vec
  )

  # ---- 10. Column map: index → anthro output columns ----------------------
  col_map <- list(
    wfa   = c("clenhei", "cbmi", "fwei", "zwei",  "zweia", "sweig"),
    hfa   = c("clenhei",         "flen", "zlen",  "zlena", "slen"),
    wfh   = c("clenhei",         "fwfl", "zwfl",  "zwfla", "swfl"),
    bmifa = c("clenhei", "cbmi", "fbmi", "zbmi",  "zbmia", "sbmi"),
    hcfa  = c(                   "fhc",  "zhc",   "zhca",  "shc"),
    acfa  = c(                   "fac",  "zac",   "zaca",  "sac"),
    tsfa  = c(                   "fts",  "zts",   "ztsa",  "sts"),
    ssfa  = c(                   "fss",  "zss",   "zssa",  "sss")
  )

  common_cols <- intersect(c("age_in_days", "age_in_months"), names(zscores_all))

  # ---- 11. Build per-index data frames + Nutritional Status ---------------
  results <- vector("list", length(return_feasible))
  names(results) <- return_feasible

  for (idx in return_feasible) {
    keep   <- intersect(c(common_cols, col_map[[idx]]), names(zscores_all))
    df_idx <- cbind(x, zscores_all[, keep, drop = FALSE])

    df_idx <- .anth_add_status(
      df         = df_idx,
      index      = idx,
      age_months = age_months,
      sex_int    = sex_int,
      weight_kg  = if (!is.null(weight_kg)) weight_kg else rep(NA_real_, n),
      height_cm  = if (!is.null(height_cm)) height_cm else rep(NA_real_, n)
    )

    class(df_idx) <- c("run_anthroindex", "data.frame")
    results[[idx]] <- df_idx
  }

  # ---- 12. Return ---------------------------------------------------------
  if (length(results) == 1L) results[[1L]] else results
}