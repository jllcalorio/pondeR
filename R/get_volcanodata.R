#' @title Get Volcano Plot Data
#'
#' @description
#' Extracts and combines results from \code{run_foldchange()} and \code{run_diff()}
#' into a single data frame. This function automatically aligns features/outcomes,
#' retrieves fold-change and significance thresholds from the provided objects,
#' and (optionally) keeps only "Up" or "Down" regulated features.
#'
#' @param ... Objects of class \code{run_foldchange} and/or \code{run_diff}.
#' @param up,down,pval Optional numeric overrides for thresholds. If \code{NULL}
#'   (default), thresholds are retrieved from the objects' parameters.
#' @param filter Logical. If \code{TRUE} (default), the data frame is filtered
#'   to include only "Up" or "Down" regulated features.
#' @param sort_es Row order by effect size: \code{"desc"} (default) puts the
#'   largest effects first, \code{"asc"} the smallest first, and \code{"none"}
#'   keeps the original order. Features are ordered by effect size magnitude
#'   (very large, large, medium, small, very small) and then by the absolute
#'   effect size estimate. Ignored when no \code{run_diff} effect sizes are
#'   available. Note that estimates from different effect size metrics (e.g.
#'   Cohen's d vs. rank-biserial correlation) are not strictly comparable.
#'
#' @return A \code{data.frame} containing merged data from both objects:
#'   fold change columns, effect size (95\% CI) and magnitude, and a single
#'   p-value column (\code{p_value_adj} when present, otherwise \code{p_value}).
#'   When the \code{run_diff} summary table is simplified, the p-value column
#'   keeps its superscript test markers (legend stored in the
#'   \code{"pvalue_legend"} attribute and shown by \code{\link{as_gt}}).
#'   Up/Down/NS classification is used only for \code{filter} and is not returned.
#'
#' @details
#' The function automatically detects the thresholds:
#' \itemize{
#'   \item \code{up} and \code{down} are retrieved from the \code{run_foldchange} parameters.
#'   \item Significance threshold (\code{pval}) is retrieved from \code{run_diff} parameters
#'         (defaulting to 0.05 if not found).
#'   \item When \code{run_diff} result has multiple outcomes (nrow > 1), a \code{p_value_adj}
#'         column is present; this function prefers it over \code{p_value} for the volcano plot.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming fc_res and diff_res are outputs from pondeR functions
#' volcano_data <- get_volcanodata(fc_res, diff_res, filter = TRUE)
#' }
#'
#' @author John Lennon L. Calorio
#' 
#' @seealso
#' \code{\link{run_foldchange}}, \code{\link{run_diff}}, \code{\link{plot_volcano}}, \code{\link[dimsprepr]{run_DIpreprocess}}
#' 
#' @export
get_volcanodata <- function(..., up = NULL, down = NULL, pval = NULL, filter = TRUE,
                           sort_es = c("desc", "asc", "none")) {

  sort_es <- match.arg(sort_es)
  dots <- list(...)
  fc_obj   <- NULL
  diff_obj <- NULL

  # 1. Identify the relevant objects (handle lists within dots)
  for (item in dots) {
    if (inherits(item, "run_foldchange") && is.null(fc_obj))   fc_obj   <- item
    if (inherits(item, "run_diff")       && is.null(diff_obj)) diff_obj <- item
  }

  # If diff_obj still NULL, check for multi-outcome list structure
  if (is.null(diff_obj)) {
    for (el in dots) {
      if (is.list(el) && !is.null(el$summary_table) &&
          is.data.frame(el$summary_table) && "p_value" %in% names(el$summary_table)) {
        diff_obj <- el
        break
      }
    }
  }

  if (is.null(fc_obj) && is.null(diff_obj)) {
    stop("At least one 'run_foldchange' or 'run_diff' object must be provided.")
  }

  # 2. Extract Data Tables and Parameters
  df_fc <- if (!is.null(fc_obj)) fc_obj$summary_table else NULL
  df_diff <- if (!is.null(diff_obj)) diff_obj$summary_table else NULL

  if (!is.null(fc_obj) && !is.null(fc_obj$params$log2) && !fc_obj$params$log2) {
    stop("get_volcanodata requires log2 fold changes in the input object.")
  }

  # Threshold resolution (passed arg > object param > default)
  up_val <- up %||% fc_obj$params$up %||% 1.5
  down_val <- down %||% fc_obj$params$down %||% 0.5
  pval_val <- pval %||% diff_obj$params$p_value %||% 
              diff_obj$params$pval %||% diff_obj$params$adj_pval %||% 0.05

  # 3. Construct base table and merge
  # prefer p_value_adj when present (multi-outcome, nrow > 1),
  # otherwise use p_value (single outcome or when no adjustment applied).
  # Only that column is kept in the output (superscript test markers retained).
  p_col <- if (!is.null(df_diff) && "p_value_adj" %in% names(df_diff)) "p_value_adj" else "p_value"
  if (!is.null(df_fc)) {
    res <- df_fc
    if (!is.null(df_diff)) {
      res[[p_col]] <- df_diff[[p_col]][match(res$feature, df_diff$outcome)]
    } else {
      warning("No 'run_diff' object provided; p-values set to 1.", call. = FALSE)
      res[[p_col]] <- 1
    }
  } else {
    warning("No 'run_foldchange' object provided; fold changes set to 1.", call. = FALSE)
    res <- df_diff
    names(res)[names(res) == "outcome"] <- "feature"
    res$fold_change <- 1
    res$log2_fc <- 0
    res[[setdiff(c("p_value", "p_value_adj"), p_col)]] <- NULL
  }

  # Add effect size (95% CI) + magnitude from run_diff
  es <- .get_diff_es(diff_obj)
  if (!is.null(es) && !"ES (95% CI)" %in% names(res) && !"effect_size_estimate" %in% names(res)) {
    idx <- match(res$feature, es$feature)
    res[["ES (95% CI)"]]  <- es[["ES (95% CI)"]][idx]
    res[["ES_magnitude"]] <- es[["ES_magnitude"]][idx]
  }
  res <- res[, c(setdiff(names(res), p_col), p_col), drop = FALSE]

  # Numeric p-values for classification only
  p_num <- .p_to_num(res[[p_col]])

  # 6. Classification logic (mimics plot_volcano)
  up_log2   <- log2(up_val)
  down_log2 <- log2(down_val)
  ok  <- !is.na(res$log2_fc) & !is.na(p_num) & p_num < pval_val
  reg <- ifelse(ok & res$log2_fc >= up_log2, "Up",
         ifelse(ok & res$log2_fc <= down_log2, "Down", "NS"))

  # 7. Filtering
  keep <- if (filter) reg %in% c("Up", "Down") else rep(TRUE, nrow(res))
  res  <- res[keep, , drop = FALSE]

  # 8. Sort by effect size magnitude, then |estimate|
  if (sort_es != "none" && !is.null(es) && nrow(res) > 1L) {
    idx <- match(res$feature, es$feature)
    mag <- .es_magnitude_rank(es$ES_magnitude[idx])
    ab  <- abs(es$estimate[idx])
    sgn <- if (sort_es == "desc") -1 else 1
    res <- res[order(sgn * mag, sgn * ab, na.last = TRUE), , drop = FALSE]
  }

  # Legend for the superscripts in the p-value column; shown by as_gt()
  # ponytail: attributes are lost by dplyr verbs / rbind, same as run_diff_summary
  leg <- attr(df_diff, "pvalue_legend", exact = TRUE)
  if (is.character(res[[p_col]]) && !is.null(leg)) attr(res, "pvalue_legend") <- leg

  return(res)
}

# Internal helper for default values
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Internal: effect size "est (lo, hi)" + magnitude per feature from a run_diff
# object (single) or its multi-outcome summary_table (raw or simplified).
# Returns data.frame(feature, `ES (95% CI)`, ES_magnitude, estimate) or NULL.
.get_diff_es <- function(diff_obj) {
  if (is.null(diff_obj)) return(NULL)
  st <- diff_obj$summary_table
  if (is.data.frame(st)) {
    if (all(c("ES (95% CI)", "ES_magnitude") %in% names(st))) {
      return(data.frame(feature = st$outcome, `ES (95% CI)` = st[["ES (95% CI)"]],
                        ES_magnitude = st$ES_magnitude,
                        # formatted "est (lo, hi)" string -> numeric estimate
                        estimate = suppressWarnings(as.numeric(sub(" \\(.*$", "", st[["ES (95% CI)"]]))),
                        check.names = FALSE, stringsAsFactors = FALSE))
    }
    if (all(c("effect_size_estimate", "effect_size_ci_low", "effect_size_ci_high",
              "effect_size_magnitude") %in% names(st))) {
      return(data.frame(feature = st$outcome,
                        `ES (95% CI)` = .fmt_es(st$effect_size_estimate,
                                                st$effect_size_ci_low, st$effect_size_ci_high),
                        ES_magnitude = st$effect_size_magnitude,
                        estimate = st$effect_size_estimate,
                        check.names = FALSE, stringsAsFactors = FALSE))
    }
    return(NULL)
  }
  es <- diff_obj$effect_size
  if (inherits(diff_obj, "run_diff") && !is.null(es)) {
    return(data.frame(feature = diff_obj$outcome,
                      `ES (95% CI)` = .fmt_es(es$estimate, es$ci_low, es$ci_high),
                      ES_magnitude = es$magnitude,
                      estimate = es$estimate,
                      check.names = FALSE, stringsAsFactors = FALSE))
  }
  NULL
}

.fmt_es <- function(est, lo, hi) {
  ifelse(is.na(est), NA_character_, sprintf("%.2f (%.2f, %.2f)", est, lo, hi))
}


# Internal: numeric p-values from run_diff(simplify = TRUE) strings such as
# ".047ᵃ" or "<.001". "<x" is only a bound, so use a value just below x.
.p_to_num <- function(p) {
  if (!is.character(p)) return(p)
  lt <- startsWith(p, "<")
  out <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", p)))
  out[lt %in% TRUE] <- out[lt %in% TRUE] * (1 - 1e-9)
  out
}


# Internal: ordinal rank of an effect size magnitude label (NA if unknown)
.es_magnitude_rank <- function(x) {
  lv <- c("tiny", "very small", "small", "medium", "large", "very large")
  r <- match(tolower(x), lv)
  # "tiny" and "very small" are the same tier
  ifelse(r == 1L, 2L, r)
}
