#' Summarise Group Sizes Relevant to Cross-Validation Fold Reduction
#'
#' @title Inspect Group Sample Sizes That May Trigger CV Fold Reduction in
#'   \code{run_regreg()}
#'
#' @description
#' Tabulates the number of samples per level of the dependent variable
#' (\code{pred}) after applying the same pre-processing steps that
#' \code{\link{run_regreg}} performs: removing rows with \code{NA} in
#' \code{pred}, and optionally subsetting rows to match those that would
#' survive the \code{not_penalized} NA removal. This lets the user see exactly
#' why \code{run_regreg} may warn that \code{cv_folds} was reduced — because
#' the smallest group has fewer samples than the requested number of folds —
#' and which group is the bottleneck.
#'
#' @details
#' The cross-validation fold reduction warning fires when
#' \code{min(table(y_train)) < cv_folds}. Because the training set is a
#' \code{train_percent} fraction of the data, the actual minimum group size
#' used in fold assignment is smaller than what appears in the full dataset.
#' This function reports both the full-data group sizes and (optionally) a
#' projected training-set minimum, so the user can anticipate the reduction
#' before running the model.
#'
#' An optional publication-ready summary table can be produced via
#' \code{\link{run_summarytable}} when \code{show_summary_table = TRUE}.
#'
#' When \code{metadata} is a named list (parallel to a \code{lapply} call),
#' each element is processed and results are stacked with a leading
#' \code{taxon} column.
#'
#' @param metadata A \code{data.frame}, \code{tibble}, \strong{or a named
#'   list} of such objects — identical to what was (or will be) passed to
#'   \code{run_regreg()}. When a list, names must be set.
#' @param pred A single character string naming the dependent variable column
#'   in \code{metadata}.
#' @param not_penalized Optional character vector. When supplied, rows with
#'   \code{NA} in any metadata-sourced \code{not_penalized} column are removed
#'   before counting, mirroring the exact sample set \code{run_regreg} would
#'   use. Default: \code{NULL} (no additional row removal).
#' @param is_numeric Optional character vector of column names listed in
#'   \code{not_penalized} that should be treated as numeric, bypassing the
#'   default factor/character detection. This mirrors the \code{is_numeric}
#'   argument of \code{\link{run_regreg}} and ensures that the NA-filtering
#'   step uses the same coercion logic, so group sizes are counted on exactly
#'   the same sample set. For example, if \code{"Serum albumin (g/L)"} is
#'   stored as character but is truly numeric, listing it here prevents it from
#'   being coerced to integer codes — which could introduce spurious NAs or
#'   hide real ones. Default: \code{NULL}.
#' @param train_percent A numeric scalar in \code{(0, 1)} matching the
#'   \code{train_percent} argument passed to \code{run_regreg()}. Used to
#'   compute the projected minimum training-group size. Default: \code{0.8}.
#' @param cv_folds Integer. The \code{cv_folds} value passed to
#'   \code{run_regreg()}. Used to flag groups whose projected training-set size
#'   is smaller than the requested folds. Default: \code{10L}.
#' @param show_summary_table Logical. If \code{TRUE} and \code{pred} is
#'   categorical, a \code{run_summarytable()} summary of \code{metadata}
#'   stratified by \code{pred} is printed as a side-effect. The numeric count
#'   data frame is still returned invisibly. Default: \code{FALSE}.
#' @param summarytable_args A named list of additional arguments forwarded to
#'   \code{\link{run_summarytable}} when \code{show_summary_table = TRUE}.
#'   Default: \code{list()}.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{\code{taxon}}{Only when \code{metadata} is a named list.}
#'   \item{\code{group}}{Level of the \code{pred} variable.}
#'   \item{\code{n_total}}{Number of samples in that group (after NA removal).}
#'   \item{\code{n_train_projected}}{Projected training-set size:
#'     \code{floor(n_total * train_percent)}.}
#'   \item{\code{folds_requested}}{The value of \code{cv_folds}.}
#'   \item{\code{fold_reduction_expected}}{Logical; \code{TRUE} when
#'     \code{n_train_projected < cv_folds}.}
#'   \item{\code{suggested_cv_folds}}{The fold count \code{run_regreg} would
#'     actually use: \code{max(3, n_train_projected)}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## ---- Single metadata data frame ------------------------------------------
#' gs <- get_groupsizes(
#'   metadata      = df_metadata_laboratory |>
#'                     dplyr::filter(`Disease Severity` != "Healthy"),
#'   pred          = "Disease Severity",
#'   not_penalized = vars_include_laboratory,
#'   is_numeric    = "Serum albumin (g/L)",
#'   train_percent = 0.6,
#'   cv_folds      = 10L
#' )
#' print(gs)
#'
#' ## ---- With a publication-ready summary table ------------------------------
#' get_groupsizes(
#'   metadata           = df_metadata_laboratory,
#'   pred               = "Disease Severity",
#'   not_penalized      = vars_include_laboratory,
#'   is_numeric         = "Serum albumin (g/L)",
#'   train_percent      = 0.6,
#'   cv_folds           = 10L,
#'   show_summary_table = TRUE
#' )
#'
#' ## ---- List of metadata (one per taxon) ------------------------------------
#' # If each taxon uses the same metadata, just pass it once (not a list).
#' # If each taxon has its own row-subsetted metadata, wrap them in a list:
#' gs_list <- get_groupsizes(
#'   metadata      = setNames(
#'     lapply(names(df_list_clr), function(nm) df_metadata_laboratory),
#'     names(df_list_clr)
#'   ),
#'   pred          = "Disease Severity",
#'   not_penalized = vars_include_laboratory,
#'   train_percent = 0.6,
#'   cv_folds      = 10L
#' )
#' # Taxa where fold reduction is expected
#' gs_list[gs_list$fold_reduction_expected, ]
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{run_regreg}}, \code{\link{get_removed_samples}},
#'   \code{\link{run_summarytable}}
#'
#' @export
get_groupsizes <- function(
    metadata,
    pred,
    not_penalized      = NULL,
    is_numeric         = NULL,
    train_percent      = 0.8,
    cv_folds           = 10L,
    show_summary_table = FALSE,
    summarytable_args  = list()
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.character(pred) || length(pred) != 1L || !nzchar(trimws(pred)))
    stop("'pred' must be a single non-empty character string.", call. = FALSE)

  if (!is.null(is_numeric)) {
    if (!is.character(is_numeric))
      stop(
        "'is_numeric' must be a character vector of column names to treat as ",
        "numeric, or NULL.",
        call. = FALSE
      )
    if (is.null(not_penalized))
      stop(
        "'is_numeric' was supplied but 'not_penalized' is NULL. ",
        "'is_numeric' only applies to columns listed in 'not_penalized'.",
        call. = FALSE
      )
    not_in_np <- setdiff(is_numeric, not_penalized)
    if (length(not_in_np) > 0L)
      warning(
        "The following column(s) in 'is_numeric' are not listed in ",
        "'not_penalized' and will be ignored: ",
        paste(not_in_np, collapse = ", "), ".",
        call. = FALSE
      )
  }

  if (!is.numeric(train_percent) || length(train_percent) != 1L ||
      !is.finite(train_percent) || train_percent <= 0 || train_percent >= 1)
    stop("'train_percent' must be a single numeric value strictly between 0 and 1.",
         call. = FALSE)

  if (!is.numeric(cv_folds) || length(cv_folds) != 1L ||
      cv_folds != round(cv_folds) || cv_folds < 3L)
    stop("'cv_folds' must be a single integer >= 3.", call. = FALSE)

  if (!is.logical(show_summary_table) || length(show_summary_table) != 1L)
    stop("'show_summary_table' must be TRUE or FALSE.", call. = FALSE)

  if (!is.list(summarytable_args))
    stop("'summarytable_args' must be a named list.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # Dispatch: list vs single
  # ---------------------------------------------------------------------------
  if (is.list(metadata) && !is.data.frame(metadata)) {

    if (is.null(names(metadata)) || any(!nzchar(names(metadata))))
      stop(
        "When 'metadata' is a list, every element must be named. ",
        "Use names(metadata) <- your_names to assign names.",
        call. = FALSE
      )

    rows <- lapply(names(metadata), function(nm) {
      out <- .ggs_single(metadata[[nm]], pred, not_penalized, is_numeric,
                         train_percent, cv_folds)
      cbind(taxon = nm, out, stringsAsFactors = FALSE)
    })

    combined <- do.call(rbind, rows)
    rownames(combined) <- NULL

    if (show_summary_table) {
      message(
        "'show_summary_table = TRUE' is not supported when 'metadata' is a ",
        "list. Call get_groupsizes() on a single metadata data frame instead."
      )
    }

    return(combined)
  }

  # Single metadata data frame
  result <- .ggs_single(metadata, pred, not_penalized, is_numeric,
                        train_percent, cv_folds)

  # ---------------------------------------------------------------------------
  # Optional summary table (single metadata only)
  # ---------------------------------------------------------------------------
  if (show_summary_table) {

    if (!is.character(metadata[[pred]]) && !is.factor(metadata[[pred]])) {
      message(
        "'show_summary_table = TRUE' is only supported for categorical 'pred' ",
        "variables. Skipping summary table."
      )
    } else if (!requireNamespace("pondeR", quietly = TRUE) &&
               !exists("run_summarytable", mode = "function")) {
      warning(
        "run_summarytable() was not found. ",
        "Make sure the pondeR package is loaded.",
        call. = FALSE
      )
    } else {

      # Filter metadata to only the rows that run_regreg would see
      meta_clean <- .ggs_filter_meta(metadata, pred, not_penalized)

      tbl_args <- c(
        list(df = meta_clean, split_by = pred),
        summarytable_args
      )

      cat("\n--- run_summarytable() output ---\n")
      tbl <- do.call(run_summarytable, tbl_args)
      print(tbl)
      cat("\n")
    }
  }

  result
}


# Internal worker: compute group sizes for one metadata data frame
.ggs_single <- function(metadata, pred, not_penalized, is_numeric,
                        train_percent, cv_folds) {

  if (!pred %in% colnames(metadata))
    stop("'pred' column \"", pred, "\" not found in 'metadata'.", call. = FALSE)

  meta_clean <- .ggs_filter_meta(metadata, pred, not_penalized)

  resp <- meta_clean[[pred]]

  # Apply is_numeric override: if pred itself is listed, coerce it
  if (!is.null(is_numeric) && pred %in% is_numeric) {
    resp <- suppressWarnings(as.numeric(as.character(resp)))
  }

  # Also coerce any not_penalized columns flagged in is_numeric, so that the
  # NA-filtering step in .ggs_filter_meta treated them correctly.
  # (Already handled upstream, but we recheck pred type here for group sizing.)

  # For numeric pred, just return total
  if (is.numeric(resp)) {
    n_valid <- sum(!is.na(resp))
    n_train <- floor(n_valid * train_percent)
    return(data.frame(
      group                   = "(numeric response)",
      n_total                 = n_valid,
      n_train_projected       = n_train,
      folds_requested         = as.integer(cv_folds),
      fold_reduction_expected = n_train < cv_folds,
      suggested_cv_folds      = max(3L, n_train),
      stringsAsFactors        = FALSE
    ))
  }

  # Categorical
  resp   <- resp[!is.na(resp)]
  counts <- sort(table(resp), decreasing = TRUE)
  lvls   <- names(counts)
  ns     <- as.integer(counts)

  n_train_proj    <- floor(ns * train_percent)
  will_reduce     <- n_train_proj < cv_folds
  suggested_folds <- pmax(3L, n_train_proj)

  data.frame(
    group                   = lvls,
    n_total                 = ns,
    n_train_projected       = n_train_proj,
    folds_requested         = as.integer(cv_folds),
    fold_reduction_expected = will_reduce,
    suggested_cv_folds      = suggested_folds,
    stringsAsFactors        = FALSE,
    row.names               = NULL
  )
}


# Internal: remove rows with NA in pred and in metadata-sourced not_penalized cols
.ggs_filter_meta <- function(metadata, pred, not_penalized, is_numeric = NULL) {

  # Drop NA in response
  keep     <- !is.na(metadata[[pred]])
  metadata <- metadata[keep, , drop = FALSE]

  if (is.null(not_penalized)) return(metadata)

  np        <- setdiff(not_penalized, pred)
  from_meta <- intersect(np, colnames(metadata))

  if (length(from_meta) == 0L) return(metadata)

  np_sub <- metadata[, from_meta, drop = FALSE]

  # Apply is_numeric coercion before the complete.cases() check,
  # mirroring what run_regreg does: columns declared numeric are coerced
  # with as.numeric() rather than as.factor(), which can change which rows
  # are NA (e.g. a column stored as character "1.5" is not NA after coercion).
  if (!is.null(is_numeric)) {
    for (col_nm in intersect(is_numeric, colnames(np_sub))) {
      np_sub[[col_nm]] <- suppressWarnings(
        as.numeric(as.character(np_sub[[col_nm]]))
      )
    }
  }

  complete <- complete.cases(np_sub)
  metadata[complete, , drop = FALSE]
}