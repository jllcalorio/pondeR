#' @title Reduce Multiple Data Frames by Common or Unique Column/Row Names
#'
#' @description
#' Accepts two or more data frames (or tibbles or matrices with dimnames) and
#' returns a single data frame whose columns (or rows) are either the
#' \emph{intersection} (names present in \strong{all} supplied objects) or the
#' \emph{union} (names present in \strong{any} supplied object) of the column
#' or row names across all inputs.
#'
#' @param ... Two or more objects coercible to a data frame: \code{data.frame},
#'   \code{tbl_df} (tibble), or named \code{matrix}.  If exactly one object is
#'   supplied a warning is issued and that object is returned as-is (after
#'   coercion to \code{data.frame}).
#' @param check_cols Logical scalar (default \code{TRUE}).  When \code{TRUE}
#'   the reduction is performed over \strong{column} names; when \code{FALSE}
#'   it is performed over \strong{row} names.  All supplied objects must carry
#'   the relevant dimnames when \code{check_cols = FALSE}.
#' @param return_common Logical scalar (default \code{TRUE}).  When \code{TRUE}
#'   only names present in \strong{all} data frames are retained
#'   (intersection).  When \code{FALSE} names present in \strong{any} data
#'   frame are retained (union); missing entries are filled with \code{NA}.
#'
#' @details
#' \subsection{Single input}{
#'   Passing a single data frame raises a warning and returns that frame
#'   unchanged (coerced to \code{data.frame}).  No reduction is possible.
#' }
#'
#' \subsection{Column-wise reduction (\code{check_cols = TRUE})}{
#'   The set of column names to keep is computed as the intersection or union
#'   of \code{colnames()} across all inputs.  Each frame is then subset (or
#'   expanded with \code{NA} columns for the union case) to the target set, and
#'   the results are row-bound with \code{do.call(rbind, ...)}.
#' }
#'
#' \subsection{Row-wise reduction (\code{check_cols = FALSE})}{
#'   The set of row names to keep is computed as the intersection or union of
#'   \code{rownames()} across all inputs.  Each frame is then subset (or
#'   expanded with \code{NA} rows for the union case) to the target set, and
#'   the results are column-bound with \code{do.call(cbind, ...)}.  All inputs
#'   must have non-\code{NULL}, non-empty, and non-duplicated row names.
#' }
#'
#' \subsection{Duplicate names}{
#'   Duplicated column names within a single data frame, or duplicated row
#'   names when \code{check_cols = FALSE}, will trigger an error because the
#'   subsetting behaviour would be ambiguous.
#' }
#'
#' \subsection{Performance note}{
#'   Name intersection/union is computed once via \code{Reduce(intersect, ...)}
#'   or \code{Reduce(union, ...)} — an idiomatic base-R approach that avoids
#'   repeated pairwise comparisons and scales well with the number of inputs.
#' }
#'
#' @return A named list of class \code{"run_reduce"} with the following
#'   elements:
#'   \describe{
#'     \item{\code{data}}{A \code{data.frame} containing the reduced result.}
#'     \item{\code{summary}}{A named list with one element per input data frame
#'       (named \code{df1}, \code{df2}, \ldots).  Each element is a character
#'       vector of the column (or row) names that were \strong{removed} from
#'       that input during the reduction, or \code{character(0)} if nothing was
#'       removed.}
#'     \item{\code{call}}{The matched call, for reproducibility.}
#'     \item{\code{check_cols}}{Logical; mirrors the \code{check_cols}
#'       argument.}
#'     \item{\code{return_common}}{Logical; mirrors the \code{return_common}
#'       argument.}
#'     \item{\code{n_inputs}}{Integer; the number of data frames supplied.}
#'     \item{\code{target_names}}{Character vector of the names retained in
#'       the final result.}
#'   }
#'
#' @examples
#' ## --- Example data ---------------------------------------------------
#' df1 <- data.frame(a = 1:3, b = 4:6,  c = 7:9)
#' df2 <- data.frame(a = 10:12, b = 13:15, d = 16:18)
#' df3 <- data.frame(a = 19:21, b = 22:24, e = 25:27)
#'
#' ## 1. Return columns common to ALL data frames (default)
#' result_common <- run_reduce(df1, df2, df3)
#' result_common$data          # columns: a, b
#' result_common$summary       # columns dropped per frame
#'
#' ## 2. Return columns from ANY data frame (union, fill with NA)
#' result_union <- run_reduce(df1, df2, df3, return_common = FALSE)
#' result_union$data           # columns: a, b, c, d, e  (NAs where absent)
#'
#' ## 3. Row-wise reduction -------------------------------------------
#' m1 <- data.frame(x = c(1, 2, 3), row.names = c("r1", "r2", "r3"))
#' m2 <- data.frame(y = c(4, 5, 6), row.names = c("r1", "r3", "r4"))
#'
#' result_rows <- run_reduce(m1, m2, check_cols = FALSE)
#' result_rows$data            # rows: r1, r3 (common to both)
#'
#' result_rows_union <- run_reduce(m1, m2, check_cols = FALSE,
#'                                 return_common = FALSE)
#' result_rows_union$data      # rows: r1, r2, r3, r4 (NAs where absent)
#'
#' ## 4. Single data frame (warning, returned as-is) ------------------
#' result_single <- run_reduce(df1)
#'
#' @author John Lennon L. Calorio
#' @export
run_reduce <- function(..., check_cols = TRUE, return_common = TRUE) {

  # -----------------------------------------------------------------------
  # 0. Capture call and collect inputs
  # -----------------------------------------------------------------------
  mc <- match.call()
  dfs_raw <- list(...)
  n <- length(dfs_raw)

  # -----------------------------------------------------------------------
  # 1. Validate scalar logical parameters
  # -----------------------------------------------------------------------
  if (!is.logical(check_cols) || length(check_cols) != 1L || is.na(check_cols)) {
    stop(
      "'check_cols' must be a single non-NA logical value (TRUE or FALSE).",
      call. = FALSE
    )
  }

  if (!is.logical(return_common) || length(return_common) != 1L ||
      is.na(return_common)) {
    stop(
      "'return_common' must be a single non-NA logical value (TRUE or FALSE).",
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # 2. Require at least one input
  # -----------------------------------------------------------------------
  if (n == 0L) {
    stop(
      "No data frames were supplied to 'run_reduce()'. ",
      "Please provide at least one data frame.",
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # 3. Coerce inputs to data.frame and validate each
  # -----------------------------------------------------------------------
  .coerce_to_df <- function(obj, idx) {
    label <- paste0("Input ", idx)

    if (is.matrix(obj)) {
      if (is.null(dimnames(obj))) {
        stop(
          label, " is a matrix with no dimnames. ",
          "Assign column and/or row names before passing to 'run_reduce()'.",
          call. = FALSE
        )
      }
      obj <- as.data.frame(obj)
    } else if (inherits(obj, c("tbl_df", "tbl", "data.frame"))) {
      obj <- as.data.frame(obj)
    } else {
      stop(
        label, " is of class '", paste(class(obj), collapse = "', '"), "'. ",
        "'run_reduce()' accepts data.frames, tibbles, or named matrices only.",
        call. = FALSE
      )
    }
    obj
  }

  dfs <- vector("list", n)
  for (i in seq_len(n)) {
    dfs[[i]] <- .coerce_to_df(dfs_raw[[i]], i)
  }

  # -----------------------------------------------------------------------
  # 4. Single-input shortcut
  # -----------------------------------------------------------------------
  if (n == 1L) {
    warning(
      "Only one data frame was supplied. ",
      "Returning it unchanged (no reduction performed).",
      call. = FALSE
    )
    return(
      structure(
        list(
          data         = dfs[[1L]],
          summary      = list(df1 = character(0L)),
          call         = mc,
          check_cols   = check_cols,
          return_common = return_common,
          n_inputs     = 1L,
          target_names = if (check_cols) colnames(dfs[[1L]]) else rownames(dfs[[1L]])
        ),
        class = "run_reduce"
      )
    )
  }

  # -----------------------------------------------------------------------
  # 5. Extract and validate the relevant names
  # -----------------------------------------------------------------------
  if (check_cols) {

    name_list <- lapply(dfs, colnames)

    ## Check for NULL column names
    null_idx <- which(vapply(name_list, is.null, logical(1L)))
    if (length(null_idx) > 0L) {
      stop(
        "The following inputs have no column names: ",
        paste(paste0("input ", null_idx), collapse = ", "), ". ",
        "All data frames must have named columns when 'check_cols = TRUE'.",
        call. = FALSE
      )
    }

    ## Check for empty column name vectors
    empty_idx <- which(vapply(name_list, function(x) length(x) == 0L, logical(1L)))
    if (length(empty_idx) > 0L) {
      stop(
        "The following inputs have zero columns: ",
        paste(paste0("input ", empty_idx), collapse = ", "), ". ",
        "All data frames must have at least one column.",
        call. = FALSE
      )
    }

    ## Check for duplicated column names within a single frame
    dup_idx <- which(vapply(name_list, function(x) anyDuplicated(x) > 0L, logical(1L)))
    if (length(dup_idx) > 0L) {
      stop(
        "The following inputs contain duplicated column names: ",
        paste(paste0("input ", dup_idx), collapse = ", "), ". ",
        "Disambiguate column names before calling 'run_reduce()'.",
        call. = FALSE
      )
    }

  } else {

    name_list <- lapply(dfs, rownames)

    ## Check for NULL row names
    null_idx <- which(vapply(name_list, is.null, logical(1L)))
    if (length(null_idx) > 0L) {
      stop(
        "The following inputs have no row names: ",
        paste(paste0("input ", null_idx), collapse = ", "), ". ",
        "All data frames must have explicit row names when 'check_cols = FALSE'.",
        call. = FALSE
      )
    }

    ## Reject default integer-sequence row names (uninformative for row reduction)
    default_rn <- function(x) identical(rownames(x), as.character(seq_len(nrow(x))))
    default_idx <- which(vapply(dfs, default_rn, logical(1L)))
    if (length(default_idx) > 0L) {
      warning(
        "The following inputs appear to use default integer row names: ",
        paste(paste0("input ", default_idx), collapse = ", "), ". ",
        "Row-wise reduction on default row names is likely unintentional. ",
        "Set meaningful row names or use 'check_cols = TRUE'.",
        call. = FALSE
      )
    }

    ## Check for duplicated row names
    dup_idx <- which(vapply(name_list, function(x) anyDuplicated(x) > 0L, logical(1L)))
    if (length(dup_idx) > 0L) {
      stop(
        "The following inputs contain duplicated row names: ",
        paste(paste0("input ", dup_idx), collapse = ", "), ". ",
        "Disambiguate row names before calling 'run_reduce()'.",
        call. = FALSE
      )
    }

  }

  # -----------------------------------------------------------------------
  # 6. Compute target name set via Reduce (intersection or union)
  # -----------------------------------------------------------------------
  target_names <- if (return_common) {
    Reduce(intersect, name_list)
  } else {
    Reduce(union, name_list)
  }

  if (length(target_names) == 0L) {
    stop(
      "The intersection of ",
      if (check_cols) "column" else "row",
      " names across all supplied data frames is empty. ",
      "No common names exist. ",
      "Consider setting 'return_common = FALSE' to take the union instead.",
      call. = FALSE
    )
  }

  # -----------------------------------------------------------------------
  # 7. Build per-frame removal summary
  # -----------------------------------------------------------------------
  df_labels <- paste0("df", seq_len(n))

  removed_summary <- setNames(
    lapply(seq_len(n), function(i) {
      setdiff(name_list[[i]], target_names)
    }),
    df_labels
  )

  # -----------------------------------------------------------------------
  # 8. Subset / expand each frame to the target name set
  # -----------------------------------------------------------------------
  if (check_cols) {

    reduced_list <- lapply(dfs, function(df) {
      missing_cols <- setdiff(target_names, colnames(df))
      if (length(missing_cols) > 0L) {
        ## Union case: add NA columns for names absent in this frame
        df[missing_cols] <- NA
      }
      df[, target_names, drop = FALSE]
    })

    result_df <- do.call(rbind, reduced_list)
    rownames(result_df) <- NULL   # reset potentially duplicated row indices

  } else {

    reduced_list <- lapply(dfs, function(df) {
      missing_rows <- setdiff(target_names, rownames(df))
      if (length(missing_rows) > 0L) {
        ## Union case: add NA rows for names absent in this frame
        na_rows <- as.data.frame(
          matrix(NA_real_, nrow = length(missing_rows), ncol = ncol(df),
                 dimnames = list(missing_rows, colnames(df)))
        )
        df <- rbind(df, na_rows)
      }
      df[target_names, , drop = FALSE]
    })

    result_df <- do.call(cbind, reduced_list)

  }

  # -----------------------------------------------------------------------
  # 9. Assemble and return the "run_reduce" object
  # -----------------------------------------------------------------------
  structure(
    list(
      data          = result_df,
      summary       = removed_summary,
      call          = mc,
      check_cols    = check_cols,
      return_common = return_common,
      n_inputs      = n,
      target_names  = target_names
    ),
    class = "run_reduce"
  )
}


# =========================================================================
# S3 print method
# =========================================================================

#' @title Print Method for \code{run_reduce} Objects
#'
#' @description Provides a concise, informative summary when a
#'   \code{run_reduce} object is printed to the console.
#'
#' @param x An object of class \code{"run_reduce"}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @author John Lennon L. Calorio
#' @export
print.run_reduce <- function(x, ...) {

  mode_label  <- if (x$check_cols) "column" else "row"
  set_label   <- if (x$return_common) "intersection (common)" else "union (all unique)"

  cat("── run_reduce result ─────────────────────────────────────────────\n")
  cat(sprintf("  Inputs      : %d data frame(s)\n", x$n_inputs))
  cat(sprintf("  Dimension   : %s names\n", mode_label))
  cat(sprintf("  Name set    : %s\n", set_label))
  cat(sprintf("  Target names: %s\n",
              paste(x$target_names, collapse = ", ")))
  cat(sprintf("  Output dim  : %d row(s) x %d col(s)\n",
              nrow(x$data), ncol(x$data)))

  any_removed <- any(vapply(x$summary, function(s) length(s) > 0L, logical(1L)))

  if (any_removed) {
    cat("\n  Names removed per input:\n")
    for (nm in names(x$summary)) {
      dropped <- x$summary[[nm]]
      if (length(dropped) > 0L) {
        cat(sprintf("    %-6s : %s\n", nm, paste(dropped, collapse = ", ")))
      } else {
        cat(sprintf("    %-6s : (none)\n", nm))
      }
    }
  } else {
    cat("\n  Names removed per input: (none in any input)\n")
  }

  cat("──────────────────────────────────────────────────────────────────\n")
  cat("\n$data:\n")
  print(x$data)

  invisible(x)
}