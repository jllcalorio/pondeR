#' Run Comprehensive Item-Level and Scale-Level Content Validity Analysis
#'
#' @description
#' Computes standard relevance-based content validity indices on a single
#' ratings matrix and returns a tidy summary. At the item level, the function
#' calculates the Item-level Content Validity Index (I-CVI), modified kappa,
#' and Aiken's V. At the scale level, it reports the Scale-level CVI by
#' averaging (S-CVI/Ave), Scale-level CVI by universal agreement (S-CVI/UA),
#' and the mean modified kappa across all items.
#'
#' @param x A numeric data frame or matrix where rows represent experts and
#'   columns represent items.
#' @param threshold An integer specifying the minimum rating value considered
#'   "relevant". Defaults to \code{3}.
#' @param low An integer specifying the minimum possible rating value on the
#'   scale. Passed to \code{aiken_v()}. Defaults to \code{1}.
#' @param high An integer specifying the maximum possible rating value on the
#'   scale. Passed to \code{aiken_v()}. Defaults to \code{4}.
#' @param na.rm Logical. If \code{TRUE}, ratings with missing values are
#'   removed prior to computation. Defaults to \code{FALSE}.
#' @param conf_level Confidence level for the score confidence interval of
#'   Aiken's V (Penfield & Giacobbi, 2004). Defaults to \code{0.95}.
#' @param agreement Logical. If \code{TRUE} (default), computes Gwet's AC1 and
#'   quadratic-weighted AC2 across all items (items as subjects, experts as
#'   raters) using the \code{irrCAC} package, with confidence intervals at
#'   \code{conf_level}. Skipped with a warning if \code{irrCAC} is not
#'   installed.
#'
#' @details
#' The function merges the unrounded numeric output of
#' \code{contentValidity::content_validity()} with the formatted column names
#' and interpretation labels from \code{contentValidity::apa_table()}, yielding
#' full-precision estimates alongside human-readable column headers.
#'
#' \strong{Modified-Kappa Cutoffs}
#'
#' Modified kappa values are interpreted using the criteria of
#' Cicchetti and Sparrow (1981), as adopted by Polit, Beck, and Owen (2007):
#'
#' \tabular{ll}{
#'   \strong{Range}      \tab \strong{Interpretation} \cr
#'   \eqn{\kappa > 0.74} \tab Excellent               \cr
#'   \eqn{0.60 \le \kappa \le 0.74} \tab Good          \cr
#'   \eqn{0.40 \le \kappa \le 0.59} \tab Fair          \cr
#'   \eqn{\kappa < 0.40} \tab Poor                     \cr
#' }
#'
#' \strong{Aiken's V Confidence Interval and Significance Test}
#'
#' Confidence intervals for Aiken's V use the score interval of Penfield and
#' Giacobbi (2004). The p-value is the exact one-tailed probability, from
#' Aiken (1985), of a rating sum at least as large as observed when each
#' rating is equally likely to fall in any scale category. With \eqn{n}
#' experts and \eqn{c} categories, the smallest attainable p-value is
#' \eqn{c^{-n}} (e.g., .0625 for 2 experts and 4 categories; .0156 for 3
#' experts), so very small panels cannot reach conventional significance.
#'
#' \strong{Inter-Expert Agreement}
#'
#' If \code{agreement = TRUE}, Gwet's AC1 (unweighted) and AC2 (quadratic
#' weights) are computed once across all items, with experts as raters, using
#' \code{irrCAC::gwet.ac1.raw()} (Gwet, 2008, 2014). They are reported at the
#' scale level because item-level coefficients based on a handful of raters
#' are degenerate.
#' 
#' @return A list of class \code{run_contentvalidity} with two elements:
#' \describe{
#'   \item{\code{$items}}{A data frame of item-level indices: I-CVI,
#'     modified kappa, Aiken's V, and Modified-Kappa Interpretation.}
#'   \item{\code{$scale}}{A named numeric vector of scale-level indices:
#'     S-CVI/Ave, S-CVI/UA, and mean modified kappa.}
#'   \item{\code{$agreement}}{A data frame with Gwet's AC1 and AC2 estimates
#'     and confidence limits, or \code{NULL} if \code{agreement = FALSE}.}
#' }
#'
#' @references
#' Aiken, L. R. (1985). Three coefficients for analyzing the reliability and
#' validity of ratings. \emph{Educational and Psychological Measurement},
#' \emph{45}(1), 131--142. \doi{10.1177/0013164485451012}
#'
#' Cicchetti, D. V., & Sparrow, S. A. (1981). Developing criteria for
#' establishing interrater reliability of specific items: Applications to
#' assessment of adaptive behavior. \emph{American Journal of Mental
#' Deficiency}, \emph{86}(2), 127--137.
#'
#' Gwet, K. L. (2008). Computing inter-rater reliability and its variance in
#' the presence of high agreement. \emph{British Journal of Mathematical and
#' Statistical Psychology}, \emph{61}(1), 29--48.
#' \doi{10.1348/000711006X126600}
#'
#' Gwet, K. L. (2014). \emph{Handbook of inter-rater reliability: The
#' definitive guide to measuring the extent of agreement among raters}
#' (4th ed.). Advanced Analytics, LLC.
#'
#' Gwet, K. L. (2019). \emph{irrCAC: Computing chance-corrected agreement
#' coefficients (CAC)} (R package).
#' \url{https://CRAN.R-project.org/package=irrCAC}
#' 
#' Lynn, M. R. (1986). Determination and quantification of content validity.
#' \emph{Nursing Research}, \emph{35}(6), 382--385.
#' \doi{10.1097/00006199-198611000-00017}
#' 
#' Penfield, R. D., & Giacobbi, Jr., P. R. (2004). Applying a Score Confidence 
#' Interval to Aiken’s Item Content-Relevance Index. 
#' \emph{Measurement in Physical Education and Exercise Science}, 
#' \emph{8}(4), 213–225. \doi{10.1207/s15327841mpee0804_3}
#' 
#' Polit, D. F., & Beck, C. T. (2006). The content validity index: Are you
#' sure you know what's being reported? Critique and recommendations.
#' \emph{Research in Nursing & Health}, \emph{29}(5), 489--497.
#' \doi{10.1002/nur.20147}
#'
#' Polit, D. F., Beck, C. T., & Owen, S. V. (2007). Is the CVI an acceptable
#' indicator of content validity? Appraisal and recommendations.
#' \emph{Research in Nursing & Health}, \emph{30}(4), 459--467.
#' \doi{10.1002/nur.20199}
#'
#' @seealso
#' \code{\link[contentValidity]{icvi}},
#' \code{\link[contentValidity]{scvi_ave}},
#' \code{\link[contentValidity]{scvi_ua}},
#' \code{\link[contentValidity]{mod_kappa}},
#' \code{\link[contentValidity]{aiken_v}}
#'
#' @examples
#' \dontrun{
#' ratings <- data.frame(
#'   item1 = c(4, 3, 4, 4, 3),
#'   item2 = c(2, 3, 2, 1, 2),
#'   item3 = c(4, 4, 3, 4, 4)
#' )
#' result <- run_contentvalidity(ratings)
#' result$items
#' result$scale
#' }
#'
#' @export
run_contentvalidity <- function(
    x,
    threshold = 3L,
    low       = 1L,
    high      = 4L,
    na.rm     = FALSE,
    conf_level = 0.95,
    agreement  = TRUE
) {

  # --- package check ----------------------------------------------------------
  if (!requireNamespace("contentValidity", quietly = TRUE)) {
    stop(
      "Package 'contentValidity' is required. ",
      "Install it with: install.packages('contentValidity')",
      call. = FALSE
    )
  }

  # --- input checks -----------------------------------------------------------
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("'x' must be a numeric data frame or matrix.", call. = FALSE)
  }

  x_mat <- as.matrix(x)

  if (!is.numeric(x_mat)) {
    stop("'x' must contain only numeric values.", call. = FALSE)
  }

  if (!na.rm && anyNA(x_mat)) {
    stop(
      "'x' contains missing values. Set 'na.rm = TRUE' to remove them, ",
      "or inspect your data before proceeding.",
      call. = FALSE
    )
  }

  if (!is.numeric(low)  || length(low)  != 1L ||
      !is.numeric(high) || length(high) != 1L) {
    stop("'low' and 'high' must each be a single numeric value.", call. = FALSE)
  }

  if (low >= high) {
    stop("'low' must be strictly less than 'high'.", call. = FALSE)
  }

  if (!is.numeric(threshold) || length(threshold) != 1L) {
    stop("'threshold' must be a single numeric value.", call. = FALSE)
  }

  if (threshold < low || threshold > high) {
    stop(
      "'threshold' must be within the rating scale range [", low, ", ", high, "].",
      call. = FALSE
    )
  }

  observed_vals <- x_mat[!is.na(x_mat)]

  if (any(observed_vals < low)) {
    stop(
      "Values below 'low' (", low, ") found in 'x'. ",
      "Verify that 'low' matches the minimum of your rating scale.",
      call. = FALSE
    )
  }

  if (any(observed_vals > high)) {
    stop(
      "Values above 'high' (", high, ") found in 'x'. ",
      "Verify that 'high' matches the maximum of your rating scale.",
      call. = FALSE
    )
  }

  if (!is.logical(na.rm) || length(na.rm) != 1L) {
    stop("'na.rm' must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(agreement) || length(agreement) != 1L) {
    stop("'agreement' must be a single logical value.", call. = FALSE)
  }

  # --- computation ------------------------------------------------------------
  cv_out  <- contentValidity::content_validity(
    ratings            = x,
    relevant_threshold = threshold,
    lo                 = low,
    hi                 = high,
    na.rm              = na.rm
  )

  apa_out <- contentValidity::apa_table(x = cv_out)

  # --- item-level assembly ----------------------------------------------------
  # apa_table() returns a flat data frame (not a list); it contains formatted
  # column names, an Item column, and an Interpretation column.
  # content_validity()$items holds the unrounded numeric values.
  # Goal: unrounded numerics + formatted col names + Interpretation (renamed).

  items_raw <- as.data.frame(cv_out$items)
  items_apa <- as.data.frame(apa_out)

  # Identify the Interpretation column
  # apa_table() names it "Kappa Interpretation" (contentValidity 0.2.0)
  interp_col <- grep("Interpretation", names(items_apa), value = TRUE)[1L]
  if (is.na(interp_col)) {
    stop("Could not find an interpretation column in apa_table() output.",
         call. = FALSE)
  }

  # Start from raw unrounded values; row names are item names
  items_out <- items_raw

  # Rename raw columns to display-friendly names
  raw_col_map <- c(
    icvi       = "I-CVI",
    mod_kappa  = "Modified Kappa",
    aiken_v    = "Aiken's V"
  )
  for (old in names(raw_col_map)) {
    if (old %in% names(items_out)) {
      names(items_out)[names(items_out) == old] <- raw_col_map[[old]]
    }
  }

  # Append interpretation under a less ambiguous name
  items_out[["Modified-Kappa Interpretation"]] <- items_apa[[interp_col]]

  # Prepend Item column from apa_out; drop any duplicate item-name column
  if ("Item" %in% names(items_apa)) {
    items_out <- cbind(Item = items_apa[["Item"]], items_out, row.names = NULL)
  }

  # Drop any residual duplicate item column (raw $items may carry item names)
  dup_cols <- which(
    names(items_out) %in% c("item") |
    (names(items_out) == "Item" & duplicated(names(items_out)))
  )
  if (length(dup_cols) > 0L) {
    items_out <- items_out[, -dup_cols, drop = FALSE]
  }

    # --- Aiken's V: score CI and exact significance test ------------------------
  k     <- high - low                      # number of scale steps (categories - 1)
  n_r   <- colSums(!is.na(x_mat))          # experts per item
  s_sum <- colSums(x_mat - low, na.rm = TRUE)
  nk    <- n_r * k
  v_obs <- s_sum / nk

  # Penfield & Giacobbi (2004) score interval
  z      <- stats::qnorm(1 - (1 - conf_level) / 2)
  centre <- 2 * nk * v_obs + z^2
  half   <- z * sqrt(4 * nk * v_obs * (1 - v_obs) + z^2)
  denom  <- 2 * (nk + z^2)
  v_lwr  <- pmax(0, (centre - half) / denom)
  v_upr  <- pmin(1, (centre + half) / denom)

  # Aiken (1985) exact one-tailed p: P(S >= observed) when each rating is
  # uniformly distributed over the (k + 1) categories
  .conv <- function(a, b) {
    out <- numeric(length(a) + length(b) - 1L)
    for (i in seq_along(a)) {
      idx <- i:(i + length(b) - 1L)
      out[idx] <- out[idx] + a[i] * b
    }
    out
  }
  .aiken_p <- function(s_obs, n, k) {
    if (n < 1L) return(NA_real_)
    dist <- Reduce(.conv, rep(list(rep(1 / (k + 1), k + 1)), n))
    sum(dist[(round(s_obs) + 1L):length(dist)])
  }
  v_p <- mapply(.aiken_p, s_sum, n_r, MoreArgs = list(k = k))

  # Safety check: recomputed V must match the package's V
  if ("Aiken's V" %in% names(items_out) &&
      !isTRUE(all.equal(unname(v_obs), unname(items_out[["Aiken's V"]])))) {
    warning("Recomputed Aiken's V differs from contentValidity output; ",
            "check item ordering.", call. = FALSE)
  }

  items_out[["Aiken's V Lower"]] <- unname(v_lwr)
  items_out[["Aiken's V Upper"]] <- unname(v_upr)
  items_out[["Aiken's V p"]]     <- unname(v_p)

  # remove the degenerate AC columns
  items_out[c("gwet_ac1", "gwet_ac2")] <- NULL
  
    # --- Inter-expert agreement across all items (Gwet's AC1 / AC2) -------------
  agreement_out <- NULL
  if (agreement) {
    if (!requireNamespace("irrCAC", quietly = TRUE)) {
      warning("Package 'irrCAC' is required for agreement coefficients; ",
              "skipping. Install it with: install.packages('irrCAC')",
              call. = FALSE)
    } else {
      r_mat <- t(x_mat)                    # items = rows, experts = columns
      cats  <- seq(low, high)
      .ci   <- function(ci) {
        as.numeric(strsplit(gsub("[()[:space:]]", "", ci), ",")[[1L]])
      }
      ac1 <- irrCAC::gwet.ac1.raw(r_mat, categ.labels = cats,
                                  conflev = conf_level)$est
      ac2 <- irrCAC::gwet.ac1.raw(r_mat, weights = "quadratic",
                                  categ.labels = cats,
                                  conflev = conf_level)$est
      ci1 <- .ci(ac1$conf.int)
      ci2 <- .ci(ac2$conf.int)
      agreement_out <- data.frame(
        Coefficient = c("Gwet's AC1", "Gwet's AC2 (quadratic)"),
        Estimate    = c(ac1$coeff.val, ac2$coeff.val),
        Lower       = c(ci1[1L], ci2[1L]),
        Upper       = c(ci1[2L], ci2[2L]),
        stringsAsFactors = FALSE
      )
    }
  }

  # --- output -----------------------------------------------------------------
  scale_out <- cv_out$scale
  scale_out <- scale_out[setdiff(names(scale_out), c("mean_ac1", "mean_ac2"))]
  scale_col_map <- c(
    scvi_ave   = "S-CVI/Ave",
    scvi_ua    = "S-CVI/UA",
    mean_kappa = "Average Kappa"
  )
  for (old in names(scale_col_map)) {
    if (old %in% names(scale_out)) {
      names(scale_out)[names(scale_out) == old] <- scale_col_map[[old]]
    }
  }

  out <- list(
    items     = items_out,
    scale     = scale_out,
    agreement = agreement_out
  )

  # summary() reads attr(object, "n_experts"), but nothing sets it
  attr(out, "n_experts") <- nrow(x_mat)

  class(out) <- c("run_contentvalidity", "list")
  out
}


# --- S3 methods ---------------------------------------------------------------

#' @export
print.run_contentvalidity <- function(x, digits = 3L, ...) {
  cat("Content Validity Analysis\n")
  cat(strrep("-", 40L), "\n\n")

  cat("Item-level indices:\n")
  items_print        <- x$items
  num_cols           <- vapply(items_print, is.numeric, logical(1L))
  items_print[num_cols] <- lapply(items_print[num_cols], round, digits = digits)
  print(items_print, row.names = TRUE)

  cat("\nScale-level indices:\n")
  scale_vals <- round(x$scale, digits = digits)
  print(scale_vals)

  if (!is.null(x$agreement)) {
    cat("\nInter-expert agreement (all items):\n")
    ag <- x$agreement
    num <- vapply(ag, is.numeric, logical(1L))
    ag[num] <- lapply(ag[num], round, digits = digits)
    print(ag, row.names = FALSE)
  }

  invisible(x)
}


#' @export
summary.run_contentvalidity <- function(object, digits = 3L, ...) {
  cat("Content Validity Analysis - Summary\n")
  cat(strrep("=", 42L), "\n\n")

  n_items   <- nrow(object$items)
  n_experts <- attr(object, "n_experts")

  cat("Items evaluated :", n_items, "\n")

  cat("\nItem-level indices (unrounded):\n")
  print(object$items, digits = digits)

  cat("\nScale-level indices:\n")
  sv <- round(object$scale, digits = digits)
  for (nm in names(sv)) {
    cat(sprintf("  %-20s %.*f\n", nm, digits, sv[[nm]]))
  }

  invisible(object)
}