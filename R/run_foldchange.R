# =============================================================================
#  run_foldchange
# =============================================================================

#' @title Fold Change Analysis Across Groups
#'
#' @description
#' Computes pairwise fold changes (and optionally log2 fold changes) for each
#' feature (column) in \code{x} across ordered group levels. Group means are
#' used as the basis of comparison. When \code{log2 = TRUE} the data are
#' shifted so that all values are strictly positive before log-transformation.
#' The function is designed to integrate seamlessly with \code{run_diff()} and
#' \code{plot_volcano()}.
#'
#' @param x A \code{data.frame}, \code{tibble}, or named \code{matrix} whose
#'   columns are numeric features (variables) to be analysed. Column names may
#'   contain special characters.
#' @param metadata A \code{data.frame}, \code{tibble}, or named \code{matrix}
#'   with the same number of rows as \code{x} and containing sample
#'   annotations. Must include the column specified by \code{group}.
#' @param group A single character string naming the column in \code{metadata}
#'   that contains the grouping factor. Must have at least two distinct
#'   non-\code{NA} levels after optional filtering via \code{filter}.
#' @param arrange An optional character vector listing \emph{all} levels of
#'   \code{group} (after filtering) in the desired comparison order. Fold
#'   changes are computed as consecutive ratios:
#'   \code{arrange[1] / arrange[2]}, \code{arrange[2] / arrange[3]}, etc.
#'   If \code{NULL} (default), levels are sorted alphanumerically.
#' @param filter An optional character vector of group levels to \emph{exclude}
#'   before analysis. Rows whose \code{group} value appears in \code{filter}
#'   are dropped. Default \code{NULL} (no filtering).
#' @param select An optional character vector of column names in \code{x} to
#'   include in the analysis. Default \code{NULL} uses all columns of \code{x}.
#' @param sort Logical. If \code{TRUE} (default), each fold-change column in
#'   the output is sorted from highest to lowest absolute fold change.
#' @param log2 Logical. If \code{TRUE} (default), log2 fold changes are also
#'   computed. Data are shifted prior to log-transformation so that the minimum
#'   value equals \code{eps} (see below), guaranteeing all values are strictly
#'   positive.
#' @param eps A small positive numeric value. When \code{log2 = TRUE}, data are
#'   shifted so that the global minimum of \code{x} (after filtering and
#'   selection) becomes \code{eps}. Default \code{1e-8}.
#'
#' @details
#' \strong{Fold change computation}\cr
#' Group means are computed per feature. For \eqn{k} ordered levels
#' \eqn{g_1, g_2, \ldots, g_k}, the fold change for comparison
#' \eqn{i} is \eqn{\bar{x}_{g_i} / \bar{x}_{g_{i+1}}}, yielding
#' \eqn{k - 1} pairwise comparisons.
#'
#' \strong{Log2 shifting}\cr
#' When \code{log2 = TRUE}, let \eqn{m} be the global minimum of the
#' (filtered, selected) data matrix. A shift of \eqn{\delta = \text{eps} - m}
#' is added to every value when \eqn{m \le 0}; otherwise no shift is applied
#' (\eqn{\delta = 0}). Log2 fold changes are then computed from the shifted
#' group means.
#'
#' \strong{Integration with \code{run_diff()} and \code{plot_volcano()}}\cr
#' The returned list exposes \code{fc_table}, \code{log2fc_table}, and
#' \code{shifted_data} in a consistent format that \code{plot_volcano()} can
#' consume directly.
#'
#' @return A named list of class \code{"run_foldchange"} containing:
#' \describe{
#'   \item{\code{fc_table}}{A \code{data.frame} of pairwise fold changes.
#'     Rows are features; columns are labelled
#'     \code{"<level_i>_vs_<level_{i+1}>"}. When \code{sort = TRUE} each
#'     column is independently sorted (descending).}
#'   \item{\code{log2fc_table}}{A \code{data.frame} of log2 fold changes in
#'     the same structure as \code{fc_table}. \code{NULL} when
#'     \code{log2 = FALSE}.}
#'   \item{\code{shifted_data}}{The (possibly shifted) numeric data matrix
#'     used for log2 computation. \code{NULL} when \code{log2 = FALSE}.}
#'   \item{\code{min_value}}{The global minimum of the analysis data
#'     \emph{before} shifting (numeric scalar).}
#'   \item{\code{shift}}{The numeric shift \eqn{\delta} applied.
#'     \code{0} when no shift was needed or \code{log2 = FALSE}.}
#'   \item{\code{group_means}}{A \code{data.frame} of per-group, per-feature
#'     means (unshifted) with one row per group level.}
#'   \item{\code{comparisons}}{A character vector of comparison labels, e.g.
#'     \code{c("A_vs_B", "B_vs_C")}.}
#'   \item{\code{params}}{A list of all resolved parameter values for
#'     reproducibility.}
#' }
#'
#' @examples
#' ## -----------------------------------------------------------------------
#' ## Example 1 — basic usage with iris
#' ## -----------------------------------------------------------------------
#' data(iris)
#' x_iris  <- iris[, 1:4]
#' meta_iris <- data.frame(Species = iris$Species)
#'
#' res <- run_foldchange(
#'   x        = x_iris,
#'   metadata = meta_iris,
#'   group    = "Species"
#' )
#' print(res$fc_table)
#' print(res$log2fc_table)
#'
#' ## -----------------------------------------------------------------------
#' ## Example 2 — custom order, filter one level, select two features
#' ## -----------------------------------------------------------------------
#' res2 <- run_foldchange(
#'   x        = x_iris,
#'   metadata = meta_iris,
#'   group    = "Species",
#'   arrange  = c("virginica", "versicolor", "setosa"),
#'   filter   = NULL,
#'   select   = c("Sepal.Length", "Petal.Length")
#' )
#' print(res2$comparisons)
#' print(res2$fc_table)
#'
#' ## -----------------------------------------------------------------------
#' ## Example 3 — data with special characters in column names
#' ## -----------------------------------------------------------------------
#' set.seed(1)
#' n   <- 60
#' xsc <- data.frame(
#'   `LPC 18:2`  = rnorm(n, 5, 1),
#'   `PA O-28:1` = rnorm(n, 3, 0.5),
#'   check.names = FALSE
#' )
#' meta_sc <- data.frame(group = rep(c("Healthy", "Disease"), each = 30))
#'
#' res3 <- run_foldchange(x = xsc, metadata = meta_sc, group = "group")
#' print(res3$fc_table)
#'
#' @author John Lennon L. Calorio
#' @export
run_foldchange <- function(
    x,
    metadata,
    group,
    arrange  = NULL,
    filter   = NULL,
    select   = NULL,
    sort     = TRUE,
    log2     = TRUE,
    eps      = 1e-8
) {

  # ---------------------------------------------------------------------------
  # 1.  Coerce x and metadata to data.frame
  # ---------------------------------------------------------------------------
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop(
      "`x` must be a data.frame, tibble, or named matrix. ",
      "Received class: ", paste(class(x), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.matrix(x)) {
    if (is.null(colnames(x)))
      stop("`x` is a matrix but has no column names.", call. = FALSE)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(metadata) && !is.matrix(metadata)) {
    stop(
      "`metadata` must be a data.frame, tibble, or named matrix. ",
      "Received class: ", paste(class(metadata), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.matrix(metadata)) metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  else                      metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (nrow(x) == 0L) stop("`x` has zero rows.",        call. = FALSE)
  if (ncol(x) == 0L) stop("`x` has zero columns.",     call. = FALSE)
  if (nrow(metadata) != nrow(x))
    stop(
      "`metadata` must have the same number of rows as `x`. ",
      "nrow(x) = ", nrow(x), "; nrow(metadata) = ", nrow(metadata), ".",
      call. = FALSE
    )

  # ---------------------------------------------------------------------------
  # 2.  Validate group
  # ---------------------------------------------------------------------------
  if (!is.character(group) || length(group) != 1L || is.na(group))
    stop("`group` must be a single non-NA character string.", call. = FALSE)
  if (!group %in% colnames(metadata))
    stop(
      "`group` = \"", group, "\" is not a column in `metadata`. ",
      "Available columns: ", paste(colnames(metadata), collapse = ", "), ".",
      call. = FALSE
    )

  # ---------------------------------------------------------------------------
  # 3.  Validate scalar logical / numeric arguments
  # ---------------------------------------------------------------------------
  .chk_bool <- function(v, nm) {
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop("`", nm, "` must be a single non-NA logical (TRUE or FALSE).", call. = FALSE)
  }
  .chk_bool(sort, "sort")
  .chk_bool(log2, "log2")

  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0)
    stop("`eps` must be a single positive numeric value.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 4.  Apply filter (exclude group levels)
  # ---------------------------------------------------------------------------
  grp_raw <- as.character(metadata[[group]])

  if (!is.null(filter)) {
    if (!is.character(filter))
      stop("`filter` must be a character vector of group levels to exclude.", call. = FALSE)
    unknown_filt <- setdiff(filter, unique(grp_raw))
    if (length(unknown_filt) > 0L)
      warning(
        "The following `filter` values are not present in `group` and will be ignored: ",
        paste(unknown_filt, collapse = ", "), ".",
        call. = FALSE
      )
    keep_rows <- !grp_raw %in% filter
    if (sum(keep_rows) == 0L)
      stop("After applying `filter`, no rows remain.", call. = FALSE)
    x        <- x[keep_rows, , drop = FALSE]
    metadata <- metadata[keep_rows, , drop = FALSE]
    grp_raw  <- grp_raw[keep_rows]
  }

  # ---------------------------------------------------------------------------
  # 5.  Apply select (subset columns of x)
  # ---------------------------------------------------------------------------
  if (!is.null(select)) {
    if (!is.character(select))
      stop("`select` must be a character vector of column names.", call. = FALSE)
    missing_sel <- setdiff(select, colnames(x))
    if (length(missing_sel) > 0L)
      stop(
        "The following `select` names are not columns in `x`: ",
        paste(missing_sel, collapse = ", "), ".",
        call. = FALSE
      )
    x <- x[, select, drop = FALSE]
  }

  # Enforce all columns numeric
  non_num <- colnames(x)[!vapply(x, is.numeric, logical(1L))]
  if (length(non_num) > 0L)
    stop(
      "All columns in `x` (after `select`) must be numeric. ",
      "Non-numeric column(s): ", paste(non_num, collapse = ", "), ".",
      call. = FALSE
    )

  # ---------------------------------------------------------------------------
  # 6.  Resolve group levels and arrange
  # ---------------------------------------------------------------------------
  uniq_levels <- sort(unique(grp_raw[!is.na(grp_raw)]))

  if (length(uniq_levels) < 2L)
    stop(
      "After filtering, `group` contains fewer than 2 distinct levels. ",
      "At least 2 levels are required for fold change computation.",
      call. = FALSE
    )

  if (!is.null(arrange)) {
    if (!is.character(arrange))
      stop("`arrange` must be a character vector.", call. = FALSE)
    missing_arr <- setdiff(arrange, uniq_levels)
    if (length(missing_arr) > 0L)
      stop(
        "The following `arrange` values are not present in `group` (after filtering): ",
        paste(missing_arr, collapse = ", "), ".",
        call. = FALSE
      )
    extra_arr <- setdiff(uniq_levels, arrange)
    if (length(extra_arr) > 0L)
      stop(
        "`arrange` must list all group levels present after filtering. ",
        "Missing from `arrange`: ", paste(extra_arr, collapse = ", "), ".",
        call. = FALSE
      )
    lvl_order <- arrange
  } else {
    lvl_order <- uniq_levels   # alphanumeric ascending (already sorted)
  }

  grp_factor <- factor(grp_raw, levels = lvl_order)

  # ---------------------------------------------------------------------------
  # 7.  Compute per-group column means (unshifted)
  # ---------------------------------------------------------------------------
  feat_cols <- colnames(x)
  n_feat    <- length(feat_cols)
  n_lvl     <- length(lvl_order)

  # Pre-allocate: rows = levels, cols = features
  means_mat <- matrix(NA_real_, nrow = n_lvl, ncol = n_feat,
                      dimnames = list(lvl_order, feat_cols))

  for (lvl in lvl_order) {
    idx <- which(grp_factor == lvl)
    if (length(idx) == 0L) next
    means_mat[lvl, ] <- colMeans(x[idx, , drop = FALSE], na.rm = TRUE)
  }

  group_means_df           <- as.data.frame(means_mat)
  group_means_df[[group]]  <- lvl_order
  group_means_df           <- group_means_df[, c(group, feat_cols), drop = FALSE]

  # ---------------------------------------------------------------------------
  # 8.  Compute fold changes (consecutive ratios of means)
  # ---------------------------------------------------------------------------
  n_comp       <- n_lvl - 1L
  comp_labels  <- paste0(lvl_order[seq_len(n_comp)], "_vs_",
                         lvl_order[seq_len(n_comp) + 1L])

  # FC matrix: rows = features, cols = comparisons
  fc_mat <- matrix(NA_real_, nrow = n_feat, ncol = n_comp,
                   dimnames = list(feat_cols, comp_labels))

  for (i in seq_len(n_comp)) {
    num <- means_mat[lvl_order[i],     ]
    den <- means_mat[lvl_order[i + 1L], ]
    fc_mat[, i] <- num / den
  }

  # ---------------------------------------------------------------------------
  # 9.  Log2 fold change (with shifting)
  # ---------------------------------------------------------------------------
  log2fc_df    <- NULL
  shifted_data <- NULL
  min_val      <- min(as.matrix(x), na.rm = TRUE)
  shift_delta  <- 0

  if (log2) {
    if (min_val <= 0) {
      shift_delta  <- eps - min_val        # ensures min becomes eps
      x_shifted    <- x + shift_delta
    } else {
      x_shifted    <- x
    }
    shifted_data <- x_shifted

    # Shifted group means
    means_shift <- matrix(NA_real_, nrow = n_lvl, ncol = n_feat,
                          dimnames = list(lvl_order, feat_cols))
    for (lvl in lvl_order) {
      idx <- which(grp_factor == lvl)
      if (length(idx) == 0L) next
      means_shift[lvl, ] <- colMeans(x_shifted[idx, , drop = FALSE], na.rm = TRUE)
    }

    log2fc_mat <- matrix(NA_real_, nrow = n_feat, ncol = n_comp,
                         dimnames = list(feat_cols, comp_labels))
    for (i in seq_len(n_comp)) {
      num <- means_shift[lvl_order[i],     ]
      den <- means_shift[lvl_order[i + 1L], ]
      log2fc_mat[, i] <- base::log2(num / den)
    }
  }

  # ---------------------------------------------------------------------------
  # 10.  Convert matrices to data frames; optionally sort each comparison col
  # ---------------------------------------------------------------------------
  .mat_to_df <- function(mat) {
    df           <- as.data.frame(mat, stringsAsFactors = FALSE)
    df$feature   <- rownames(mat)
    rownames(df) <- NULL
    df[, c("feature", comp_labels), drop = FALSE]
  }

  fc_df     <- .mat_to_df(fc_mat)
  log2fc_df <- if (log2) .mat_to_df(log2fc_mat) else NULL

  if (sort && n_comp >= 1L) {
    # Sort each comparison independently by descending absolute value
    fc_df <- do.call(cbind, c(
      list(data.frame(feature = fc_df$feature, stringsAsFactors = FALSE)),
      lapply(comp_labels, function(cl) {
        ord <- order(abs(fc_df[[cl]]), decreasing = TRUE, na.last = TRUE)
        df_out           <- data.frame(fc_df$feature[ord], fc_df[[cl]][ord],
                                       stringsAsFactors = FALSE)
        colnames(df_out) <- c(paste0("feature_", cl), cl)
        df_out
      })
    ))
    # Rebuild a tidy sorted frame: one feature column per comparison (keep wide)
    # Simpler: just return sorted-by-first-comparison for the feature column
    ord_primary <- order(abs(fc_mat[, 1L]), decreasing = TRUE, na.last = TRUE)
    fc_df       <- .mat_to_df(fc_mat[ord_primary, , drop = FALSE])

    if (log2) {
      log2fc_df <- .mat_to_df(log2fc_mat[ord_primary, , drop = FALSE])
    }
  }

  # ---------------------------------------------------------------------------
  # 11.  Assemble and return
  # ---------------------------------------------------------------------------
  structure(
    list(
      fc_table     = fc_df,
      log2fc_table = log2fc_df,
      shifted_data = shifted_data,
      min_value    = min_val,
      shift        = shift_delta,
      group_means  = group_means_df,
      comparisons  = comp_labels,
      params       = list(
        group    = group,
        arrange  = lvl_order,
        filter   = filter,
        select   = if (!is.null(select)) select else feat_cols,
        sort     = sort,
        log2     = log2,
        eps      = eps,
        n_feat   = n_feat,
        n_groups = n_lvl
      )
    ),
    class = "run_foldchange"
  )
}


# =============================================================================
#  S3 print method
# =============================================================================

#' @title Print Method for \code{run_foldchange} Objects
#' @description Compact console summary of a \code{run_foldchange} result.
#' @param x An object of class \code{"run_foldchange"}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @author John Lennon L. Calorio
#' @export
print.run_foldchange <- function(x, ...) {
  cat("── run_foldchange results ──────────────────────────────────\n")
  cat(sprintf("  Features   : %d\n",  x$params$n_feat))
  cat(sprintf("  Groups     : %d  (%s)\n",
              x$params$n_groups, paste(x$params$arrange, collapse = " > ")))
  cat(sprintf("  Comparisons: %s\n",  paste(x$comparisons, collapse = ", ")))
  cat(sprintf("  log2 FC    : %s\n",  x$params$log2))
  if (x$params$log2 && x$shift != 0)
    cat(sprintf("  Shift (δ)  : %.3e  (min before shift: %.3e)\n",
                x$shift, x$min_value))
  cat("\n  Fold-change table (first 6 rows):\n")
  print(utils::head(x$fc_table), row.names = FALSE)
  if (!is.null(x$log2fc_table)) {
    cat("\n  log2 FC table (first 6 rows):\n")
    print(utils::head(x$log2fc_table), row.names = FALSE)
  }
  invisible(x)
}