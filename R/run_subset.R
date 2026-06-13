#' Balanced Subsetting of Rows or Columns
#'
#' @description
#' Subsets a data frame, tibble, or matrix either row-wise or column-wise,
#' with optional balancing across categorical groupings. When balancing is
#' requested via `match`, the function ensures that every category within
#' each specified column is represented by the same number of observations.
#' An optional `n` argument acts as a row ceiling: the per-category count
#' is derived so that the chosen set is as close to `n` as possible while
#' remaining perfectly balanced within each `match` column.
#'
#' When `match` contains more than one column, the `independent` argument
#' controls how they are combined:
#' \itemize{
#'   \item \strong{`independent = TRUE` (default):} Each `match` column is
#'     balanced separately. A row is selected if it satisfies the balance
#'     requirement of **any** `match` column (union). This maximises the
#'     number of chosen rows.
#'   \item \strong{`independent = FALSE`:} All `match` columns are crossed
#'     into a single Cartesian product of categories (e.g.
#'     `"1-2 years × Male"`), and every unique combination is represented
#'     equally. This guarantees perfect joint balance.
#' }
#'
#' @param x A data frame, tibble, or matrix to be subsetted.
#' @param n An integer specifying the desired number of rows (when
#'   `by = "row"`) or columns (when `by = "col"`) in `$chosen`.
#' @param match `NULL` (default) or a character vector of one or more
#'   column names in `x` containing categorical variables to balance on.
#' @param independent Logical; `TRUE` (default). When `match` has multiple columns:
#'   If `TRUE`, each column is balanced independently and their results are unioned.
#'   If `FALSE`, columns are crossed into a single interaction and balanced jointly.
#' @param balance A string; either `"within_match"` (default) or `"within_group"`.
#'   Controls the direction of balancing when both `match` and `group` are provided.
#'   \itemize{
#'     \item \strong{`"within_match"`}: For every category in `match`, the levels of `group`
#'       will have an equal number of observations. For example, if balancing by Age Group,
#'       each Age Group will have an equal number of Diseased and Healthy patients.
#'     \item \strong{`"within_group"`}: For every category in `group`, the levels of `match`
#'       will have an equal number of observations. For example, within the Diseased
#'       group, there will be an equal number of patients from each Age Group.
#'   }
#' @param group `NULL` (default) or a character vector of column names
#'   that define additional grouping strata. 
#' @param by A string; either `"row"` (default) or `"col"`. 
#' @param random Logical; `TRUE` (default) to select rows or columns randomly.
#' @param top_to_bottom Logical; `TRUE` (default). For deterministic row pulls.
#' @param left_to_right Logical; `TRUE` (default). For deterministic col pulls.
#' @param seed An integer; `123` (default). The random seed. `NULL` disables seeding.
#'
#' @importFrom("utils", "tail")
#' 
#' @author John Lennon L. Calorio
#' 
#' @return A named list of class `"run_subset"`.
#' @export
run_subset <- function(
    x,
    n              = NULL,
    match          = NULL,
    independent    = TRUE,
    balance        = "within_match",
    group          = NULL,
    by             = "row",
    random         = TRUE,
    top_to_bottom  = TRUE,
    left_to_right  = TRUE,
    seed           = 123
) {

  # ── 0. Input validation ─────────────────────────────────────────────────────
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("`x` must be a data frame, tibble, or matrix.", call. = FALSE)
  }

  was_matrix <- is.matrix(x)
  if (was_matrix) {
    x <- as.data.frame(x)
    message("Note: `x` was a matrix and has been coerced to a data frame.")
  }

  if (nrow(x) == 0L) stop("`x` has 0 rows. Nothing to subset.", call. = FALSE)
  if (ncol(x) == 0L) stop("`x` has 0 columns. Nothing to subset.", call. = FALSE)

  by <- match.arg(by, choices = c("row", "col"))
  balance <- match.arg(balance, choices = c("within_match", "within_group"))

  .assert_scalar_logical <- function(val, nm) {
    if (!is.logical(val) || length(val) != 1L || is.na(val))
      stop(sprintf("`%s` must be a single TRUE or FALSE.", nm), call. = FALSE)
  }
  .assert_scalar_logical(random, "random")
  .assert_scalar_logical(top_to_bottom, "top_to_bottom")
  .assert_scalar_logical(left_to_right, "left_to_right")
  .assert_scalar_logical(independent, "independent")

  if (!is.null(match)) {
    missing_cols <- setdiff(match, names(x))
    if (length(missing_cols) > 0) stop("Missing `match` columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.null(group)) {
    missing_cols <- setdiff(group, names(x))
    if (length(missing_cols) > 0) stop("Missing `group` columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # ── 1. Seed Management ──────────────────────────────────────────────────────
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) stop("`seed` must be a valid integer.", call. = FALSE)
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(as.integer(seed))
  }

  # ── 2. Row subsetting ───────────────────────────────────────────────────────
  if (by == "row") {
    all_idx <- seq_len(nrow(x))

    # Helper: Mathematical Engine for Balancing
    .get_balanced_idx <- function(target_n, m_vec, g_vec) {
      unique_m <- unique(m_vec)
      unique_g <- unique(g_vec)
      n_m <- length(unique_m)
      n_g <- length(unique_g)

      if (n_m == 0 || n_g == 0) return(integer(0))
      selected <- integer(0)

      if (balance == "within_match") {
        # 1. Find the bottleneck (maximum balanced capacity) for each match level
        max_k_per_m <- setNames(numeric(n_m), unique_m)
        for (mv in unique_m) {
          avail <- sapply(unique_g, function(gv) sum(m_vec == mv & g_vec == gv))
          max_k_per_m[mv] <- min(avail)
        }

        # 2. Proportionally scale down ONLY if target_n requires it
        target_k_per_m <- max_k_per_m
        if (!is.null(target_n)) {
          total_possible_rows <- sum(max_k_per_m) * n_g
          if (target_n < total_possible_rows) {
            target_pairs <- floor(target_n / n_g)
            props <- max_k_per_m / sum(max_k_per_m)
            target_k_per_m <- floor(props * target_pairs)
            
            # Distribute remainder to categories with highest fractional loss
            remainder <- target_pairs - sum(target_k_per_m)
            if (remainder > 0) {
               frac <- (props * target_pairs) - target_k_per_m
               add_idx <- names(sort(frac, decreasing = TRUE))[seq_len(remainder)]
               target_k_per_m[add_idx] <- target_k_per_m[add_idx] + 1
            }
          }
        }

        # 3. Pull the rows
        for (mv in unique_m) {
          k <- target_k_per_m[mv]
          if (k > 0) {
            for (gv in unique_g) {
              pool <- all_idx[m_vec == mv & g_vec == gv]
              chosen <- if (random) sample(pool, k) else if (top_to_bottom) head(pool, k) else tail(pool, k)
              selected <- c(selected, chosen)
            }
          }
        }
        
      } else {
        # balance == "within_group"
        # 1. Find the bottleneck (maximum balanced capacity) for each group level
        max_k_per_g <- setNames(numeric(n_g), unique_g)
        for (gv in unique_g) {
          avail <- sapply(unique_m, function(mv) sum(m_vec == mv & g_vec == gv))
          max_k_per_g[gv] <- min(avail)
        }

        # 2. Proportionally scale down ONLY if target_n requires it
        target_k_per_g <- max_k_per_g
        if (!is.null(target_n)) {
          total_possible_rows <- sum(max_k_per_g) * n_m
          if (target_n < total_possible_rows) {
            target_pairs <- floor(target_n / n_m)
            props <- max_k_per_g / sum(max_k_per_g)
            target_k_per_g <- floor(props * target_pairs)
            
            # Distribute remainder
            remainder <- target_pairs - sum(target_k_per_g)
            if (remainder > 0) {
               frac <- (props * target_pairs) - target_k_per_g
               add_idx <- names(sort(frac, decreasing = TRUE))[seq_len(remainder)]
               target_k_per_g[add_idx] <- target_k_per_g[add_idx] + 1
            }
          }
        }

        # 3. Pull the rows
        for (gv in unique_g) {
          k <- target_k_per_g[gv]
          if (k > 0) {
            for (mv in unique_m) {
              pool <- all_idx[m_vec == mv & g_vec == gv]
              chosen <- if (random) sample(pool, k) else if (top_to_bottom) head(pool, k) else tail(pool, k)
              selected <- c(selected, chosen)
            }
          }
        }
      }
      return(selected)
    }

    # Group Vector Prep
    g_vec <- if (!is.null(group)) do.call(paste, c(x[group], sep = "|||")) else rep("all", nrow(x))

    if (!is.null(match)) {
      if (!independent && length(match) > 1) {
        # Joint Cross-product
        m_vec <- do.call(paste, c(x[match], sep = " x "))
        chosen_idx <- .get_balanced_idx(n, m_vec, g_vec)
      } else {
        # Independent draws unioned together
        sub_n <- if (!is.null(n)) floor(n / length(match)) else NULL
        idx_list <- lapply(match, function(m_col) {
          .get_balanced_idx(sub_n, as.character(x[[m_col]]), g_vec)
        })
        chosen_idx <- unique(unlist(idx_list))
      }
    } else {
      # No Match: Proportional / Random Stratification
      if (!is.null(group)) {
        unique_g <- unique(g_vec)
        counts <- table(g_vec)
        target_counts <- counts
        
        if (!is.null(n)) {
          props <- counts / sum(counts)
          target_counts <- floor(props * n)
          remainder <- n - sum(target_counts)
          if (remainder > 0) {
            frac <- (props * n) - target_counts
            add_idx <- names(sort(frac, decreasing = TRUE))[seq_len(remainder)]
            target_counts[add_idx] <- target_counts[add_idx] + 1
          }
        }

        chosen_idx <- integer(0)
        for (gv in unique_g) {
          pool <- all_idx[g_vec == gv]
          k <- min(target_counts[gv], length(pool))
          if (k > 0) chosen_idx <- c(chosen_idx, if (random) sample(pool, k) else if (top_to_bottom) head(pool, k) else tail(pool, k))
        }
      } else {
        # Simple random sample
        k <- if (!is.null(n)) min(n, nrow(x)) else nrow(x)
        chosen_idx <- if (random) sample(all_idx, k) else if (top_to_bottom) head(all_idx, k) else tail(all_idx, k)
      }
    }

    if (length(chosen_idx) > 1 && random) chosen_idx <- sample(chosen_idx) # Final shuffle
    
    # Mathematical Ceiling Check
    if (!is.null(n) && length(chosen_idx) < n) {
      warning(sprintf(
        "Could not reach target n = %d due to categorical balancing constraints, division math, or bottlenecks. Returned maximum balanced size: %d.",
        n, length(chosen_idx)
      ), call. = FALSE)
    }

    not_chosen_idx <- setdiff(all_idx, chosen_idx)
    result_chosen <- x[chosen_idx, , drop = FALSE]
    result_not_chosen <- if (length(not_chosen_idx) > 0L) x[not_chosen_idx, , drop = FALSE] else NULL
    
    # Generate actual balance summary
    balance_summary <- if (!is.null(match)) {
      list(match_table = table(result_chosen[match]), group_table = if (!is.null(group)) table(result_chosen[group]) else NULL)
    } else NULL

  # ── 3. Column subsetting ────────────────────────────────────────────────────
  } else {
    all_col_idx <- seq_len(ncol(x))
    k <- if (!is.null(n)) min(n, length(all_col_idx)) else length(all_col_idx)
    
    chosen_col_idx <- if (random) {
      sort(sample(all_col_idx, size = k))
    } else {
      if (left_to_right) head(all_col_idx, k) else tail(all_col_idx, k)
    }

    not_chosen_col_idx <- setdiff(all_col_idx, chosen_col_idx)
    result_chosen <- x[, chosen_col_idx, drop = FALSE]
    result_not_chosen <- if (length(not_chosen_col_idx) > 0L) x[, not_chosen_col_idx, drop = FALSE] else NULL
    balance_summary <- NULL
  }

  # ── 4. Parameter Output ─────────────────────────────────────────────────────
  parameters <- list(
    by = by, n_requested = n, 
    n_chosen = if (by == "row") nrow(result_chosen) else ncol(result_chosen),
    match = match, independent = independent, group = group, random = random, 
    top_to_bottom = top_to_bottom, left_to_right = left_to_right, seed = seed,
    was_matrix = was_matrix, balance_summary = balance_summary
  )

  out <- list(chosen = result_chosen, not_chosen = result_not_chosen, parameters = parameters)
  class(out) <- "run_subset"
  return(out)
}