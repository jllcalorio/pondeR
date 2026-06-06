#' @title Randomize a Data Frame, Matrix, or Vector
#'
#' @description
#' Randomly samples rows from a data frame, matrix, or elements from a vector,
#' with optional group-stratified sampling. Serves as a structured wrapper
#' around \code{\link[dplyr]{slice_sample}}, extending it with flexible
#' per-group control (joint or sequential/independent stratification),
#' reproducible seeding, and strict input validation suitable for biological
#' and clinical research workflows.
#'
#' @param x A data frame, tibble, matrix, or atomic vector. Column names with
#'   special characters are supported. Matrices are coerced to data frames
#'   internally and the result is recast to a matrix on return.
#' @param n A single positive integer specifying the number of rows (or
#'   elements, if \code{x} is a vector) to return. Defaults to \code{NULL}.
#'   Mutually exclusive with \code{p}; if both are supplied, \code{n} takes
#'   priority and a warning is issued. When \code{per_group} is supplied
#'   without an explicit \code{n}/\code{p} inside the list, this top-level
#'   value is used as the fallback.
#' @param p A single numeric value in \code{[0, 1]} specifying the proportion
#'   of rows (or elements) to return. Defaults to \code{NULL} (equivalent to
#'   \code{p = 1}, i.e., all rows shuffled). Mutually exclusive with \code{n}.
#'   When \code{per_group} is supplied without an explicit \code{n}/\code{p}
#'   inside the list, this top-level value is used as the fallback.
#' @param per A character string specifying the randomization strategy for 
#'   data frames, matrices, or tibbles. 
#'   \itemize{
#'     \item \strong{"all"} (default): Will randomize all of the items in the data frame. Any value can be placed anywhere in the randomized data frame.
#'     \item \strong{"row"}: Will randomize the items row-wise. The randomization is done for each row independently.
#'     \item \strong{"col"}: Will randomize the items column-wise. The randomization is done for each column independently.
#'   }
#' @param per_group An optional list of per-group sampling specs controlling
#'   stratified sampling. Each element must follow one of two forms:
#'   \itemize{
#'     \item \strong{Single column}: \code{list("ColumnName", n = <integer>)}
#'       or \code{list("ColumnName", p = <number>)}
#'     \item \strong{Multiple columns (joint or independent)}: \cr
#'       \code{list(c("Col1", "Col2"), n = <integer>)} or \cr
#'       \code{list(c("Col1", "Col2"), p = <number>)} or \cr
#'       \code{list(c("Col1", "Col2"), n = c(2, 3))} (per-column values,
#'       for \code{joint = FALSE} + \code{joint_behaviour = "independent"} only)
#'   }
#'   The \code{n}/\code{p} inside each list element is optional. If omitted,
#'   the top-level \code{n} or \code{p} argument is used as the fallback
#'   (top-level \code{n} takes priority over top-level \code{p} if both are
#'   provided). If both \code{n} and \code{p} are supplied inside a list
#'   element, \code{n} takes priority and a warning is issued.
#'
#'   Requires \code{x} to be a data frame or matrix. Grouping columns must be
#'   categorical (character, factor, or logical). See \code{joint} and
#'   \code{joint_behaviour} for how multiple columns are handled.
#' @param joint A single logical value. Defaults to \code{TRUE}.
#'   \itemize{
#'     \item \code{TRUE}: all columns listed in a \code{per_group} entry are
#'       treated \emph{simultaneously} — rows are grouped by every unique
#'       combination of levels across those columns (via
#'       \code{\link[base]{interaction}}), and \code{n}/\code{p} is enforced
#'       for each combination. For example, with \code{cyl} and \code{gear},
#'       a temporary group key such as \code{"6:4"} or \code{"8:3"} is created
#'       internally and discarded after sampling.
#'     \item \code{FALSE}: the \code{joint_behaviour} argument controls whether
#'       columns are handled in \emph{cascading} or \emph{independent} fashion.
#'   }
#'   Ignored when only a single column is listed in a \code{per_group} entry.
#' @param joint_behaviour A character string, either \code{"cascading"} (or
#'   \code{"c"}) or \code{"independent"} (or \code{"i"}). Only applies when
#'   \code{joint = FALSE} and a \code{per_group} entry lists multiple columns.
#'   Defaults to \code{"cascading"}.
#'   \itemize{
#'     \item \code{"cascading"}: columns are applied sequentially in the order
#'       listed. The data are first grouped by column 1 and sampled; then the
#'       \emph{reduced result} is grouped by column 2 and sampled again, and
#'       so on. Each step's \code{n}/\code{p} is applied to whatever rows
#'       remain after the previous step, which can aggressively reduce the
#'       output when \code{n} is small.
#'     \item \code{"independent"}: each column is sampled independently from
#'       the \emph{original} \code{x}. The results are combined by taking the
#'       \emph{union} of selected row indices, so a row qualifies if it was
#'       chosen under \emph{any} of the grouping columns. When per-column
#'       \code{n}/\code{p} vectors are supplied (e.g.,
#'       \code{n = c(2, 3)}), each value corresponds positionally to each
#'       column. Vectors must not exceed the number of columns in length.
#'   }
#' @param replace A single logical value. If \code{FALSE} (default), sampling
#'   is performed without replacement. If \code{TRUE}, rows or elements may be
#'   repeated in the output.
#' @param seed A single integer passed to \code{\link[base]{set.seed}} for
#'   reproducibility. Defaults to \code{123L}. Supply \code{NA} to skip
#'   seeding (non-reproducible runs).
#'
#' @details
#' ## Sampling modes
#'
#' \code{run_randomize()} operates in three mutually exclusive top-level modes,
#' resolved in the following order of priority:
#'
#' \enumerate{
#'   \item \strong{Per-group stratified sampling} (when \code{per_group} is not
#'     \code{NULL}): Rows are sampled within groups defined by the column(s) in
#'     each \code{per_group} entry. The \code{joint} and \code{joint_behaviour}
#'     arguments control how multiple columns are handled (see below).
#'   \item \strong{Fixed-count sampling} (\code{n} supplied, \code{per_group =
#'     NULL}): Exactly \code{n} rows are returned. Setting \code{n} equal to
#'     \code{nrow(x)} returns all rows in a new random order.
#'   \item \strong{Proportional sampling} (\code{p} supplied, or both \code{n}
#'     and \code{p} omitted): A fraction \code{p} of rows is returned. Omitting
#'     both defaults to \code{p = 1} (full shuffle). \code{p = 0} returns a
#'     zero-row object of the same class.
#' }
#'
#' ## Per-group syntax
#'
#' Each element of \code{per_group} bundles one or more grouping columns with
#' an optional sampling quantity:
#' \preformatted{
#' # Single column, explicit n
#' per_group = list(
#'   list("Species", n = 10)
#' )
#'
#' # Two columns, joint (default), explicit n shared across all combinations
#' per_group = list(
#'   list(c("cyl", "gear"), n = 2)
#' )
#'
#' # Two columns, independent, per-column n vector
#' per_group = list(
#'   list(c("cyl", "gear"), n = c(2, 3))
#' )
#'
#' # Omit n/p to fall back on the top-level n or p argument
#' per_group = list(
#'   list("Species")
#' )
#' }
#'
#' ## Behaviour of \code{joint} and \code{joint_behaviour}
#'
#' \tabular{lll}{
#'   \strong{joint} \tab \strong{joint_behaviour} \tab \strong{Effect} \cr
#'   \code{TRUE}  \tab (ignored)         \tab Groups by all columns simultaneously via \code{interaction()}; one \code{n}/\code{p} for every combination \cr
#'   \code{FALSE} \tab \code{"cascading"} \tab Applies each column's \code{n}/\code{p} sequentially on the progressively reduced data \cr
#'   \code{FALSE} \tab \code{"independent"} \tab Samples each column independently from the original data; returns the union of selected rows \cr
#' }
#'
#' ## Output row order
#'
#' When \code{per_group} is used, the final output is randomly shuffled so that
#' rows from the same group are not clustered together. For example, if grouping
#' by \code{Species}, the returned rows will be intermixed across species rather
#' than appearing as contiguous blocks of \emph{setosa}, then \emph{versicolor},
#' and so on. This mimics a fully randomized layout and avoids order effects in
#' downstream analyses.
#'
#' @return
#' An object of the same class as \code{x}:
#' \itemize{
#'   \item \strong{Data frame / tibble}: a data frame or tibble with sampled
#'     rows, preserving all columns and their classes. Row names are dropped
#'     (consistent with \code{\link[dplyr]{slice_sample}} behaviour).
#'   \item \strong{Matrix}: a matrix with sampled rows, preserving column names
#'     and storage mode.
#'   \item \strong{Vector}: an atomic vector of sampled elements, preserving
#'     names if present.
#' }
#'
#' @examples
#' ## ------------------------------------------------------------------
#' ## 1. Shuffle all rows (default: p = 1)
#' ## ------------------------------------------------------------------
#' run_randomize(mtcars)
#'
#' ## ------------------------------------------------------------------
#' ## 2. Sample a fixed number of rows
#' ## ------------------------------------------------------------------
#' run_randomize(mtcars, n = 5)
#'
#' ## ------------------------------------------------------------------
#' ## 3. Sample a proportion of rows
#' ## ------------------------------------------------------------------
#' run_randomize(mtcars, p = 0.25)
#'
#' ## ------------------------------------------------------------------
#' ## 4. Sampling with replacement
#' ## ------------------------------------------------------------------
#' run_randomize(mtcars, n = 10, replace = TRUE)
#'
#' ## ------------------------------------------------------------------
#' ## 5. Randomize a vector
#' ## ------------------------------------------------------------------
#' run_randomize(1:20, n = 8, seed = 42)
#'
#' run_randomize(c(a = 1, b = 2, c = 3, d = 4), p = 0.5, seed = 7)
#'
#' ## ------------------------------------------------------------------
#' ## 6. Per-group: single column, explicit n
#' ## ------------------------------------------------------------------
#' run_randomize(
#'   iris,
#'   per_group = list(list("Species", n = 10)),
#'   seed = 2024
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 7. Per-group: single column, explicit p
#' ## ------------------------------------------------------------------
#' run_randomize(
#'   iris,
#'   per_group = list(list("Species", p = 0.6)),
#'   seed = 2024
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 8. Per-group: single column, n/p omitted — falls back to top-level n
#' ## ------------------------------------------------------------------
#' run_randomize(
#'   iris,
#'   n = 5,
#'   per_group = list(list("Species")),
#'   seed = 2024
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 9. Per-group: two columns, joint = TRUE (default)
#' ##    Samples n = 2 from every unique cyl x gear combination.
#' ##    Expected combinations and row counts in mtcars:
#' ##      4:3 (1), 4:4 (8), 4:5 (2), 6:3 (2), 6:4 (4),
#' ##      6:5 (1), 8:3 (12), 8:5 (2)
#' ##    => returns 1+2+2+2+2+1+2+2 = 14 rows
#' ##       (4:3 and 6:5 are undersized; all their rows are returned)
#' ## ------------------------------------------------------------------
#' mtcars2 <- mtcars
#' mtcars2$cyl  <- factor(mtcars2$cyl)
#' mtcars2$gear <- factor(mtcars2$gear)
#'
#' run_randomize(
#'   mtcars2,
#'   per_group = list(list(c("cyl", "gear"), n = 2)),
#'   joint = TRUE,
#'   seed = 99
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 10. Per-group: two columns, joint = FALSE, cascading (default)
#' ##     Step 1 — group by cyl, take n = 2 per level => 6 rows
#' ##     Step 2 — group by gear, take n = 3 per level from those 6 rows
#' ##     Final row count depends on how many gear levels survived step 1
#' ## ------------------------------------------------------------------
#' run_randomize(
#'   mtcars2,
#'   per_group = list(list(c("cyl", "gear"), n = c(2, 3))),
#'   joint = FALSE,
#'   joint_behaviour = "cascading",
#'   seed = 99
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 11. Per-group: two columns, joint = FALSE, independent
#' ##     cyl sampled with n = 2 per level from original data  => row set A
#' ##     gear sampled with n = 3 per level from original data => row set B
#' ##     Returns union of A and B (rows qualifying under either column)
#' ## ------------------------------------------------------------------
#' run_randomize(
#'   mtcars2,
#'   per_group = list(list(c("cyl", "gear"), n = c(2, 3))),
#'   joint = FALSE,
#'   joint_behaviour = "independent",
#'   seed = 99
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 12. Per-group: column name with special characters
#' ## ------------------------------------------------------------------
#' df <- data.frame(
#'   `Age Group` = c("Young", "Young", "Middle", "Middle", "Old", "Old"),
#'   Score       = c(80, 85, 70, 75, 90, 95),
#'   check.names = FALSE
#' )
#'
#' run_randomize(
#'   df,
#'   per_group = list(list("Age Group", n = 1)),
#'   seed = 1L
#' )
#'
#' ## ------------------------------------------------------------------
#' ## 13. If both n and p are supplied, n takes priority (with warning)
#' ## ------------------------------------------------------------------
#' run_randomize(iris, n = 10, p = 0.5, seed = 1L)
#'
#' @author John Lennon L. Calorio
#' @export
run_randomize <- function(
    x,
    n               = NULL,
    p               = NULL,
    per             = c("all", "row", "col"),
    per_group       = NULL,
    joint           = TRUE,
    joint_behaviour = "cascading",
    replace         = FALSE,
    seed            = 123L
) {

  # ── 0. Dependency check ──────────────────────────────────────────────────────
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "'dplyr' is required by run_randomize() but is not installed.\n",
      "Install it with: install.packages(\"dplyr\")",
      call. = FALSE
    )
  }

  # ── 1. Validate `x` ──────────────────────────────────────────────────────────
  if (missing(x) || is.null(x)) {
    stop(
      "`x` must be supplied and cannot be NULL.\n",
      "Provide a data frame, matrix, or atomic vector.",
      call. = FALSE
    )
  }

  is_vec <- is.atomic(x) && !is.matrix(x)
  is_mat <- is.matrix(x)
  is_df  <- is.data.frame(x)

  if (!is_vec && !is_mat && !is_df) {
    stop(
      "`x` must be a data frame, tibble, matrix, or atomic vector.\n",
      "Received an object of class: ", paste(class(x), collapse = ", "),
      call. = FALSE
    )
  }

  if (is_vec && length(x) == 0L) {
    warning("`x` is an empty vector. Returning `x` unchanged.", call. = FALSE)
    return(x)
  }

  if ((is_df || is_mat) && nrow(x) == 0L) {
    warning("`x` has zero rows. Returning `x` unchanged.", call. = FALSE)
    return(x)
  }

  # Validate 'per'
  per <- match.arg(per)
  
  mat_colnames <- NULL
  mat_mode     <- NULL
  if (is_mat) {
    mat_colnames <- colnames(x)
    mat_mode     <- storage.mode(x)
    x            <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  # ── 2. Validate `joint` ──────────────────────────────────────────────────────
  if (!is.logical(joint) || length(joint) != 1L || is.na(joint)) {
    stop(
      "`joint` must be a single non-NA logical value (TRUE or FALSE).",
      call. = FALSE
    )
  }

  # ── 3. Validate `joint_behaviour` ────────────────────────────────────────────
  valid_jb <- c("cascading", "c", "independent", "i")
  if (!is.character(joint_behaviour) || length(joint_behaviour) != 1L ||
      !joint_behaviour %in% valid_jb) {
    stop(
      "`joint_behaviour` must be one of: \"cascading\" (\"c\") or ",
      "\"independent\" (\"i\").\n",
      "Received: ", deparse(joint_behaviour),
      call. = FALSE
    )
  }
  joint_behaviour <- if (joint_behaviour %in% c("cascading", "c")) {
    "cascading"
  } else {
    "independent"
  }

  # ── 4. Validate `replace` ────────────────────────────────────────────────────
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop(
      "`replace` must be a single non-NA logical value (TRUE or FALSE).",
      call. = FALSE
    )
  }

  # ── 5. Validate `seed` ───────────────────────────────────────────────────────
  if (length(seed) != 1L || (!is.na(seed) && !is.numeric(seed))) {
    stop(
      "`seed` must be a single integer or NA (to skip seeding).\n",
      "Example: seed = 42L",
      call. = FALSE
    )
  }
  if (!is.na(seed)) {
    if (seed != as.integer(seed)) {
      warning(
        "`seed` was coerced to integer: ", as.integer(seed), ".",
        call. = FALSE
      )
    }
    set.seed(as.integer(seed))
  }

  # ── 6. Resolve top-level n / p ───────────────────────────────────────────────
  n_rows          <- if (is_vec) length(x) else nrow(x)
  n_supplied      <- !is.null(n)
  p_supplied      <- !is.null(p)
  p_user_supplied <- p_supplied   # track explicit user intent

  if (!n_supplied && !p_supplied) {
    p          <- 1
    p_supplied <- TRUE
  }

  if (n_supplied && p_supplied) {
    warning(
      "Both `n` and `p` were supplied. `n` takes priority; `p` is ignored.\n",
      "To use `p`, omit the `n` argument.",
      call. = FALSE
    )
    p_supplied      <- FALSE
    p_user_supplied <- FALSE
    p               <- NULL
  }

  if (n_supplied) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        n < 1 || n != as.integer(n)) {
      stop(
        "`n` must be a single positive integer (e.g., n = 5).\n",
        "Received: ", deparse(n),
        call. = FALSE
      )
    }
    n <- as.integer(n)
    if (!replace && n > n_rows) {
      stop(
        "`n` (", n, ") exceeds the number of rows/elements in `x` (",
        n_rows, ").\n",
        "Either reduce `n`, or set `replace = TRUE` to allow repeated sampling.",
        call. = FALSE
      )
    }
  }

  if (p_supplied) {
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
      stop(
        "`p` must be a single numeric value between 0 and 1 (inclusive).\n",
        "Received: ", deparse(p),
        call. = FALSE
      )
    }
  }

  # ── 7. Validate and parse `per_group` ────────────────────────────────────────
  pg_specs <- NULL

  if (!is.null(per_group)) {

    if (is_vec) {
      stop(
        "`per_group` is not supported when `x` is a vector.\n",
        "Supply a data frame or matrix to use stratified sampling.",
        call. = FALSE
      )
    }

    if (!is.list(per_group) || length(per_group) == 0L) {
      stop(
        "`per_group` must be a non-empty list.\n",
        "Example: per_group = list(list(\"Species\", n = 10))",
        call. = FALSE
      )
    }

    # Fallback quantity from top-level n/p
    fallback_type  <- if (n_supplied) "n" else if (p_user_supplied) "p" else "p"
    fallback_value <- if (n_supplied) n   else if (p_user_supplied) p   else 1

    pg_specs <- .parse_per_group(per_group, x, fallback_type, fallback_value)

    # Validate grouping columns are categorical
    for (spec in pg_specs) {
      for (col in spec$cols) {
        col_vec <- x[[col]]
        if (!is.character(col_vec) && !is.factor(col_vec) &&
            !is.logical(col_vec)) {
          stop(
            "Grouping column '", col, "' must be categorical ",
            "(character, factor, or logical).\n",
            "Column class: ", paste(class(col_vec), collapse = ", "), ".\n",
            "Convert it first, e.g.: x[[\"", col, "\"]] <- as.factor(x[[\"",
            col, "\"]])",
            call. = FALSE
          )
        }
      }
    }
  }

  # ── 7.5 Dimension-wise randomization ─────────────────────────────────────────
  if (!is_vec) {
    m <- as.matrix(x)
    if (per == "row") {
      # Randomize the items row-wise
      for (i in seq_len(nrow(m))) m[i, ] <- m[i, sample.int(ncol(m))]
    } else if (per == "col") {
      # Randomize the items column-wise
      for (j in seq_len(ncol(m))) m[, j] <- m[sample.int(nrow(m)), j]
    } else if (per == "all") {
      # Randomize all of the items in the data frame
      m[] <- sample(m)
    }
    # Assign back to preserve structure/class
    x[] <- as.data.frame(m, stringsAsFactors = FALSE)
  }

  # ── 8. Execute sampling ───────────────────────────────────────────────────────
  if (!is.null(pg_specs)) {

    result <- x
    for (spec in pg_specs) {
      n_cols <- length(spec$cols)
      if (n_cols == 1L || joint) {
        result <- .sample_joint(result, spec, replace)
      } else if (joint_behaviour == "cascading") {
        result <- .sample_cascading(result, spec, replace)
      } else {
        result <- .sample_independent(x, spec, replace)
      }
    }

    # Shuffle so same-group rows are not contiguous
    result <- dplyr::slice_sample(result, prop = 1, replace = FALSE)

  } else if (is_vec) {
    size   <- if (n_supplied) n else max(0L, round(p * n_rows))
    idx    <- sample.int(n_rows, size = size, replace = replace)
    result <- x[idx]

  } else {
    result <- if (n_supplied) {
      dplyr::slice_sample(x, n = n, replace = replace)
    } else {
      dplyr::slice_sample(x, prop = p, replace = replace)
    }
  }

  # ── 9. Restore matrix class ───────────────────────────────────────────────────
  if (is_mat) {
    result               <- as.matrix(result)
    storage.mode(result) <- mat_mode
    colnames(result)     <- mat_colnames
  }

  result
}


# ── Internal helpers ────────────────────────────────────────────────────────────

#' Parse per_group into a structured spec list
#'
#' Each spec has: cols (character vector), types (character vector of "n"/"p"),
#' values (numeric vector, positionally aligned with cols).
#' @noRd
.parse_per_group <- function(per_group, x, fallback_type, fallback_value) {

  df_cols <- names(x)

  lapply(seq_along(per_group), function(i) {

    el       <- per_group[[i]]
    el_names <- names(el)

    if (!is.list(el)) {
      stop(
        "`per_group` element ", i, " is not a list.\n",
        "Example: list(\"Species\", n = 10)  or  list(c(\"cyl\", \"gear\"), n = 2)",
        call. = FALSE
      )
    }

    # ── Column name(s): first unnamed element ──────────────────────────────────
    unnamed_idx <- which(is.null(el_names) | el_names == "")
    cols <- if (length(unnamed_idx) >= 1L) el[[unnamed_idx[1L]]] else NULL

    if (is.null(cols) || !is.character(cols) || length(cols) == 0L ||
        any(nchar(cols) == 0L)) {
      stop(
        "`per_group` element ", i, " is missing a valid column name.\n",
        "Supply column name(s) as the first (unnamed) element:\n",
        "  list(\"ColumnName\", n = 3)\n",
        "  list(c(\"Col1\", \"Col2\"), n = 2)",
        call. = FALSE
      )
    }

    bad_cols <- cols[!cols %in% df_cols]
    if (length(bad_cols) > 0L) {
      stop(
        "Column(s) not found in `x` (per_group element ", i, "): ",
        paste(paste0("'", bad_cols, "'"), collapse = ", "), ".\n",
        "Available columns: ", paste(df_cols, collapse = ", "),
        call. = FALSE
      )
    }

    n_cols   <- length(cols)
    has_n    <- "n" %in% el_names
    has_p    <- "p" %in% el_names

    # ── No n/p supplied: fall back to top-level ────────────────────────────────
    if (!has_n && !has_p) {
      return(list(
        cols   = cols,
        types  = rep(fallback_type,  n_cols),
        values = rep(fallback_value, n_cols)
      ))
    }

    # ── Both n and p supplied: warn, n wins ────────────────────────────────────
    if (has_n && has_p) {
      warning(
        "`per_group` element ", i, " supplies both `n` and `p`. ",
        "`n` takes priority; `p` is ignored.",
        call. = FALSE
      )
      has_p <- FALSE
    }

    param_name <- if (has_n) "n" else "p"
    raw_val    <- el[[param_name]]

    # ── Validate value: scalar or vector ──────────────────────────────────────
    if (!is.numeric(raw_val) || any(is.na(raw_val))) {
      stop(
        "`", param_name, "` in `per_group` element ", i,
        " must be a numeric vector with no NA values.\n",
        "Received: ", deparse(raw_val),
        call. = FALSE
      )
    }

    if (length(raw_val) > n_cols) {
      stop(
        "Length of `", param_name, "` in `per_group` element ", i,
        " (", length(raw_val), ") exceeds the number of columns supplied (",
        n_cols, ").\n",
        "Either reduce the length of `", param_name, "` to at most ", n_cols,
        ", or add more column names to match.",
        call. = FALSE
      )
    }

    # Recycle scalar to match number of columns
    values <- if (length(raw_val) == 1L) rep(raw_val, n_cols) else raw_val

    if (param_name == "n") {
      bad <- which(values < 1 | values != as.integer(values))
      if (length(bad) > 0L) {
        stop(
          "`n` value(s) at position(s) ", paste(bad, collapse = ", "),
          " in `per_group` element ", i,
          " must be positive integer(s).\n",
          "Received: ", deparse(values[bad]),
          call. = FALSE
        )
      }
      values <- as.integer(values)
    } else {
      bad <- which(values < 0 | values > 1)
      if (length(bad) > 0L) {
        stop(
          "`p` value(s) at position(s) ", paste(bad, collapse = ", "),
          " in `per_group` element ", i,
          " must be in [0, 1].\n",
          "Received: ", deparse(values[bad]),
          call. = FALSE
        )
      }
    }

    list(
      cols   = cols,
      types  = rep(param_name, n_cols),
      values = values
    )
  })
}


#' Warn when groups are undersized (without replacement, fixed n)
#' @noRd
.warn_undersized <- function(data, col, n_req, replace) {
  if (replace || is.na(n_req)) return(invisible(NULL))
  grp_sizes  <- tapply(seq_len(nrow(data)), data[[col]], length)
  small_grps <- names(grp_sizes)[grp_sizes < n_req]
  if (length(small_grps) > 0L) {
    warning(
      "Some groups in column '", col, "' have fewer rows than n = ", n_req, ":\n",
      paste0("  - '", small_grps, "' (", grp_sizes[small_grps], " row(s))",
             collapse = "\n"), "\n",
      "All rows from these groups will be returned.\n",
      "Set `replace = TRUE` to allow oversampling.",
      call. = FALSE
    )
  }
}


#' Sample one spec using interaction() — joint or single-column
#' @noRd
.sample_joint <- function(data, spec, replace) {

  cols   <- spec$cols
  type   <- spec$types[1L]
  value  <- spec$values[1L]

  if (length(cols) == 1L) {
    col_sym <- as.name(cols)
    if (type == "n") .warn_undersized(data, cols, value, replace)
    grouped <- dplyr::group_by(data, !!col_sym)
  } else {
    # Build interaction key
    interact_col <- ".pg_joint_key_"
    data[[interact_col]] <- do.call(
      interaction,
      c(lapply(cols, function(col) data[[col]]), list(drop = TRUE, sep = ":"))
    )
    if (type == "n") {
      grp_sizes  <- tapply(seq_len(nrow(data)), data[[interact_col]], length)
      small_grps <- names(grp_sizes)[grp_sizes < value]
      if (length(small_grps) > 0L && !replace) {
        warning(
          "Some joint groups have fewer rows than n = ", value, ":\n",
          paste0("  - '", small_grps, "' (", grp_sizes[small_grps], " row(s))",
                 collapse = "\n"), "\n",
          "All rows from these groups will be returned.\n",
          "Set `replace = TRUE` to allow oversampling.",
          call. = FALSE
        )
      }
    }
    key_sym <- as.name(interact_col)
    grouped <- dplyr::group_by(data, !!key_sym)
  }

  result <- if (type == "n") {
    dplyr::slice_sample(grouped, n = value, replace = replace)
  } else {
    dplyr::slice_sample(grouped, prop = value, replace = replace)
  }

  result <- dplyr::ungroup(result)
  if (length(cols) > 1L) result[[".pg_joint_key_"]] <- NULL
  result
}


#' Cascading stratified sampling: apply each column sequentially
#' @noRd
.sample_cascading <- function(data, spec, replace) {

  result <- data

  for (j in seq_along(spec$cols)) {
    col   <- spec$cols[j]
    type  <- spec$types[j]
    value <- spec$values[j]

    # Drop unused factor levels from prior step
    if (is.factor(result[[col]])) result[[col]] <- droplevels(result[[col]])

    if (type == "n") .warn_undersized(result, col, value, replace)

    col_sym <- as.name(col)
    grouped <- dplyr::group_by(result, !!col_sym)
    result  <- if (type == "n") {
      dplyr::slice_sample(grouped, n = value, replace = replace)
    } else {
      dplyr::slice_sample(grouped, prop = value, replace = replace)
    }
    result <- dplyr::ungroup(result)
  }

  result
}


#' Independent stratified sampling: sample each column from original data,
#' return the union of selected rows.
#' @noRd
.sample_independent <- function(data, spec, replace) {

  n_rows      <- nrow(data)
  selected    <- logical(n_rows)   # TRUE = row selected by at least one column

  for (j in seq_along(spec$cols)) {
    col   <- spec$cols[j]
    type  <- spec$types[j]
    value <- spec$values[j]

    if (type == "n") .warn_undersized(data, col, value, replace)

    col_sym   <- as.name(col)
    grouped   <- dplyr::group_by(data, !!col_sym)
    sampled   <- if (type == "n") {
      dplyr::slice_sample(grouped, n = value, replace = replace)
    } else {
      dplyr::slice_sample(grouped, prop = value, replace = replace)
    }
    sampled   <- dplyr::ungroup(sampled)

    # Mark selected rows by matching row positions
    row_idx   <- which(
      duplicated(rbind(data, sampled), fromLast = TRUE)[seq_len(n_rows)]
    )
    selected[row_idx] <- TRUE
  }

  data[selected, , drop = FALSE]
}