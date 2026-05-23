#' @title Get Volcano Plot Data
#'
#' @description
#' Extracts and combines results from \code{run_foldchange()} and \code{run_diff()}
#' into a single data frame. This function automatically aligns features/outcomes,
#' retrieves fold-change and significance thresholds from the provided objects,
#' and classifies each feature as "Up", "Down", or "NS" (Non-significant).
#'
#' @param ... Objects of class \code{run_foldchange} and/or \code{run_diff}.
#' @param up,down,pval Optional numeric overrides for thresholds. If \code{NULL}
#'   (default), thresholds are retrieved from the objects' parameters.
#' @param filter Logical. If \code{TRUE} (default), the data frame is filtered
#'   to include only "Up" or "Down" regulated features.
#'
#' @return A \code{data.frame} containing merged data from both objects,
#'   including columns for fold change, p-values, and a \code{.regulation}
#'   classification column.
#'
#' @details
#' The function automatically detects the thresholds:
#' \itemize{
#'   \item \code{up} and \code{down} are retrieved from the \code{run_foldchange} parameters.
#'   \item Significance threshold (\code{pval}) is retrieved from \code{run_diff} parameters
#'         (defaulting to 0.05 if not found).
#'   \item It automatically prefers \code{adj_p_value} over \code{p_value} if available
#'         in the \code{run_diff} summary table.
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
#' \code{\link{run_foldchange}}, \code{\link{run_diff}}, \code{\link{plot_volcano}}, \code{\link{run_DIpreprocess}}
#' 
#' @export
get_volcanodata <- function(..., up = NULL, down = NULL, pval = NULL, filter = TRUE) {

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
  if (!is.null(df_fc)) {
    res <- df_fc
    if (!is.null(df_diff)) {
      p_col <- if ("adj_p_value" %in% names(df_diff)) "adj_p_value" else "p_value"
      idx <- match(res$feature, df_diff$outcome)
      res$p_value_used <- df_diff[[p_col]][idx]
      res$p_type <- p_col
    } else {
      warning("No 'run_diff' object provided; p-values set to 1.", call. = FALSE)
      res$p_value_used <- 1
      res$p_type <- "none"
    }
  } else {
    warning("No 'run_foldchange' object provided; fold changes set to 1.", call. = FALSE)
    res <- df_diff
    names(res)[names(res) == "outcome"] <- "feature"
    res$fold_change <- 1
    res$log2_fc <- 0
    p_col <- if ("adj_p_value" %in% names(df_diff)) "adj_p_value" else "p_value"
    res$p_value_used <- df_diff[[p_col]]
    res$p_type <- p_col
  }

  # Handle zero p-values to match plot_volcano logic
  zero_p <- !is.na(res$p_value_used) & res$p_value_used == 0
  if (any(zero_p)) {
    res$p_value_used[zero_p] <- .Machine$double.eps
  }

  # 6. Classification logic (mimics plot_volcano)
  up_log2   <- log2(up_val)
  down_log2 <- log2(down_val)

  res$.regulation <- ifelse(
    !is.na(res$log2_fc) & !is.na(res$p_value_used) & 
      res$log2_fc >= up_log2 & res$p_value_used < pval_val, "Up",
    ifelse(
      !is.na(res$log2_fc) & !is.na(res$p_value_used) & 
        res$log2_fc <= down_log2 & res$p_value_used < pval_val, "Down",
      "NS"
    )
  )

  # 7. Filtering
  if (filter) {
    res <- res[res$.regulation %in% c("Up", "Down"), , drop = FALSE]
  }

  return(res)
}

# Internal helper for default values
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}