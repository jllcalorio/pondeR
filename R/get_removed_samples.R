#' Identify Samples Removed Due to Missing Values in Unpenalized Covariates
#'
#' @title Get Samples Removed by \code{run_regreg()} Due to Missing Covariate Data
#'
#' @description
#' Replicates the NA-detection logic of \code{\link{run_regreg}} to identify
#' which samples would be (or were) removed because of missing values in one or
#' more columns supplied to \code{not_penalized}. Returns a data frame of those
#' samples together with the columns responsible, so the user can inspect,
#' impute, or exclude them deliberately before re-running the model.
#'
#' The function operates entirely on the original inputs — no fitted model
#' object is needed. This means it can be called before \code{run_regreg()} to
#' pre-screen data, or after to audit a run.
#'
#' @details
#' The function applies the same filtering sequence as \code{run_regreg()}:
#' \enumerate{
#'   \item The dependent variable column (\code{pred}) is excluded from
#'     \code{not_penalized} (it is never a covariate).
#'   \item Only columns in \code{not_penalized} that come from \code{metadata}
#'     (i.e., not already in \code{x}) are checked, because columns already in
#'     \code{x} are filtered earlier by the feature-NA check.
#'   \item Rows with at least one \code{NA} in those columns are flagged.
#' }
#'
#' When \code{x} is a named list (e.g., from a \code{lapply} call over
#' taxonomic levels), the function iterates over each list element and returns
#' a combined data frame with a leading \code{taxon} column.
#'
#' @param x A \code{data.frame}, \code{tibble}, \code{matrix}, \strong{or a
#'   named list} of such objects — matching what was passed (or would be
#'   passed) to \code{run_regreg()}. When a list is supplied, each element is
#'   processed independently.
#' @param metadata A \code{data.frame} or \code{tibble} of sample metadata,
#'   identical to the one passed to \code{run_regreg()}.
#' @param not_penalized A character vector of column names (from \code{metadata}
#'   and/or \code{x}) that were passed to the \code{not_penalized} argument of
#'   \code{run_regreg()}.
#' @param pred A single character string naming the dependent variable column
#'   in \code{metadata} (same as the \code{pred} argument of
#'   \code{run_regreg()}). Used only to exclude the response column from the
#'   covariate NA check.
#' @param include_metadata_cols Logical. If \code{TRUE} (default), all columns
#'   of \code{metadata} are appended to the returned data frame for context. If
#'   \code{FALSE}, only the columns specified in \code{not_penalized} (plus
#'   \code{taxon} when \code{x} is a list) are returned.
#' @param row_numbers Logical. If \code{TRUE}, a column \code{.row} containing
#'   the original row index of each flagged sample is prepended to the result.
#'   Default: \code{TRUE}.
#'
#' @return A \code{data.frame} of flagged samples. Columns include:
#' \describe{
#'   \item{\code{taxon}}{Only present when \code{x} is a named list. The name
#'     of the list element from which the sample originated.}
#'   \item{\code{.row}}{Original row index in \code{metadata} (when
#'     \code{row_numbers = TRUE}).}
#'   \item{\code{.missing_in}}{A semicolon-separated string naming the
#'     \code{not_penalized} columns that are \code{NA} for that sample.}
#'   \item{Metadata columns}{All columns of \code{metadata} (when
#'     \code{include_metadata_cols = TRUE}).}
#' }
#' Returns an empty \code{data.frame} (zero rows) with a message if no samples
#' are flagged.
#'
#' @examples
#' \dontrun{
#' ## ---- Single run ----------------------------------------------------------
#' removed <- get_removed_samples(
#'   x             = df_class_rel_abund_clr,
#'   metadata      = df_metadata_demog,
#'   not_penalized = names(df_metadata_demog),
#'   pred          = "Disease Severity"
#' )
#' print(removed)
#'
#' ## ---- lapply run (list of taxa) -------------------------------------------
#' removed_list <- get_removed_samples(
#'   x             = df_list_clr,
#'   metadata      = df_metadata_laboratory,
#'   not_penalized = vars_include_laboratory,
#'   pred          = "Disease Severity"
#' )
#' # See which taxa have the most removed samples
#' table(removed_list$taxon)
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{run_regreg}}, \code{\link{get_groupsizes}}
#'
#' @export
get_removed_samples <- function(
    x,
    metadata,
    not_penalized,
    pred,
    include_metadata_cols = TRUE,
    row_numbers           = TRUE
) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.data.frame(metadata) && !is.matrix(metadata))
    stop("'metadata' must be a data.frame or tibble.", call. = FALSE)

  if (!is.character(not_penalized) || length(not_penalized) == 0L)
    stop("'not_penalized' must be a non-empty character vector of column names.",
         call. = FALSE)

  if (!is.character(pred) || length(pred) != 1L || !nzchar(trimws(pred)))
    stop("'pred' must be a single non-empty character string.", call. = FALSE)

  if (!pred %in% colnames(metadata))
    stop("'pred' column \"", pred, "\" not found in 'metadata'.", call. = FALSE)

  if (!is.logical(include_metadata_cols) || length(include_metadata_cols) != 1L)
    stop("'include_metadata_cols' must be TRUE or FALSE.", call. = FALSE)

  if (!is.logical(row_numbers) || length(row_numbers) != 1L)
    stop("'row_numbers' must be TRUE or FALSE.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # Dispatch: list vs single
  # ---------------------------------------------------------------------------
  if (is.list(x) && !is.data.frame(x)) {

    if (is.null(names(x)) || any(!nzchar(names(x))))
      stop(
        "When 'x' is a list, every element must be named (e.g. by taxon). ",
        "Use names(x) <- your_names to assign names.",
        call. = FALSE
      )

    results <- lapply(names(x), function(nm) {
      out <- .grs_single(x[[nm]], metadata, not_penalized, pred,
                         include_metadata_cols, row_numbers)
      if (nrow(out) > 0L) cbind(taxon = nm, out, stringsAsFactors = FALSE)
      else NULL
    })

    combined <- do.call(rbind, Filter(Negate(is.null), results))

    if (is.null(combined) || nrow(combined) == 0L) {
      message("No samples with missing values found across any list element.")
      return(invisible(data.frame()))
    }
    rownames(combined) <- NULL
    return(combined)
  }

  # Single data frame / matrix
  .grs_single(x, metadata, not_penalized, pred,
              include_metadata_cols, row_numbers)
}


# Internal worker for one x / metadata pair
.grs_single <- function(x, metadata, not_penalized, pred,
                        include_metadata_cols, row_numbers) {

  x_cols <- if (is.matrix(x)) colnames(x) else colnames(x)

  # Mirror run_regreg: silently drop pred from not_penalized
  np <- setdiff(not_penalized, pred)

  # Only metadata-sourced columns matter here (x-sourced are handled earlier
  # in run_regreg's feature-NA check and are not relevant to this warning)
  from_meta <- intersect(np, colnames(metadata))
  from_meta <- setdiff(from_meta, x_cols)  # exclude cols that exist in x

  if (length(from_meta) == 0L) {
    message(
      "None of the 'not_penalized' columns originate exclusively from ",
      "'metadata'. No metadata-sourced NA check is applicable."
    )
    return(invisible(data.frame()))
  }

  np_sub   <- metadata[, from_meta, drop = FALSE]
  na_flags <- is.na(np_sub)               # logical matrix: rows x cols
  flagged  <- which(rowSums(na_flags) > 0L)

  if (length(flagged) == 0L) {
    message("No samples with missing values in the 'not_penalized' metadata columns.")
    return(invisible(data.frame()))
  }

  # Build the column listing which columns are NA per row
  missing_in <- apply(na_flags[flagged, , drop = FALSE], 1L, function(r) {
    paste(from_meta[r], collapse = "; ")
  })

  out <- if (include_metadata_cols) {
    metadata[flagged, , drop = FALSE]
  } else {
    metadata[flagged, from_meta, drop = FALSE]
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)

  # Prepend helper columns
  if (row_numbers) {
    out <- cbind(.row = flagged, .missing_in = missing_in, out,
                 stringsAsFactors = FALSE)
  } else {
    out <- cbind(.missing_in = missing_in, out, stringsAsFactors = FALSE)
  }

  rownames(out) <- NULL
  out
}