# =============================================================================
# run_correl.R
# Part of the pondeR package
# =============================================================================

# Internal helpers -------------------------------------------------------------

#' @keywords internal
.detect_binary <- function(v) {
  uv <- unique(stats::na.omit(v))
  length(uv) == 2L
}

#' @keywords internal
.detect_ordinal <- function(v) {
  is.ordered(v) || (is.factor(v) && !.detect_binary(as.integer(v)))
}

#' @keywords internal
.is_continuous <- function(v) {
  is.numeric(v) && !.detect_binary(v) && !is.integer(v)
}

#' @keywords internal
.count_ties <- function(v) {
  tb  <- table(v)
  sum(tb[tb > 1L])
}

#' @keywords internal
.auto_method <- function(vx, vy,
                          norm_result_x, norm_result_y,
                          var_violated,
                          alpha_normality = 0.05) {

  bin_x  <- .detect_binary(vx)
  bin_y  <- .detect_binary(vy)
  ord_x  <- .detect_ordinal(vx)
  ord_y  <- .detect_ordinal(vy)
  cont_x <- .is_continuous(vx)
  cont_y <- .is_continuous(vy)
  n      <- sum(!is.na(vx) & !is.na(vy))
  ties   <- (.count_ties(stats::na.omit(vx)) + .count_ties(stats::na.omit(vy))) / (2 * n)

  norm_x_viol <- isTRUE(norm_result_x$violated)
  norm_y_viol <- isTRUE(norm_result_y$violated)

  if (bin_x && bin_y)  return("phi")
  if (xor(bin_x, bin_y) && (cont_x || cont_y)) return("pointbiserial")
  if ((ord_x && ord_y)) {
    if (requireNamespace("polycor", quietly = TRUE)) return("poly")
    return("spearman")
  }
  if (xor(ord_x, ord_y) && (cont_x || cont_y)) {
    if (requireNamespace("polycor", quietly = TRUE)) return("poly")
    return("spearman")
  }
  if ((!is.ordered(vx) && is.factor(vx) && nlevels(vx) > 2L) ||
      (!is.ordered(vy) && is.factor(vy) && nlevels(vy) > 2L)) return("cramerv")
  if (cont_x && cont_y) {
    if (norm_x_viol || norm_y_viol) {
      if (n < 30 || ties > 0.2) return("kendall")
      return("spearman")
    }
    if (var_violated) return("spearman")
    return("pearson")
  }
  return("spearman")
}

#' @keywords internal
.check_normality_correl <- function(data_vector,
                                     label          = NULL,
                                     method         = "auto",
                                     alpha_normality = 0.05) {
  n <- length(stats::na.omit(data_vector))
  if (n < 3L)
    return(list(violated = TRUE,  p_value = NA_real_,
                method   = "Insufficient data", label = label))

  use_method <- if (method == "auto") {
    if (n < 50L) "shapiro" else "lilliefors"
  } else method

  if (use_method == "shapiro") {
    res   <- stats::shapiro.test(stats::na.omit(data_vector))
    p_val <- res$p.value
    return(list(violated = p_val < alpha_normality, p_value = p_val,
                method   = "Shapiro-Wilk test", label = label))
  }

  if (!requireNamespace("nortest", quietly = TRUE)) {
    warning("Package 'nortest' not available; falling back to Shapiro-Wilk.",
            call. = FALSE)
    res   <- stats::shapiro.test(stats::na.omit(data_vector))
    p_val <- res$p.value
    return(list(violated = p_val < alpha_normality, p_value = p_val,
                method   = "Shapiro-Wilk test (fallback)", label = label))
  }

  res   <- nortest::lillie.test(stats::na.omit(data_vector))
  p_val <- res$p.value
  list(violated = p_val < alpha_normality, p_value = p_val,
       method   = "Lilliefors (K-S) test", label = label)
}

#' @keywords internal
.correl_pair <- function(vx, vy, xname, yname, method,
                          sig_threshold, normality_method, alpha_normality,
                          cor_args) {

  complete_idx <- stats::complete.cases(vx, vy)
  vx_c <- vx[complete_idx]
  vy_c <- vy[complete_idx]
  n    <- length(vx_c)

  if (n < 3L) {
    return(list(
      estimate   = NA_real_, p_value = NA_real_,
      ci_lower   = NA_real_, ci_upper = NA_real_,
      method_used = method,  n = n,
      sig        = FALSE,    xname = xname, yname = yname
    ))
  }

  norm_x <- .check_normality_correl(vx_c, xname, normality_method,
                                     alpha_normality)
  norm_y <- .check_normality_correl(vy_c, yname, normality_method,
                                     alpha_normality)

  var_violated <- FALSE
  if (is.numeric(vx_c) && is.numeric(vy_c)) {
    var_violated <- tryCatch({
      if (requireNamespace("car", quietly = TRUE)) {
        grp <- rep(c("x", "y"), c(n, n))
        car::leveneTest(c(vx_c, vy_c), factor(grp))$`Pr(>F)`[1L] < alpha_normality
      } else {
        stats::bartlett.test(list(vx_c, vy_c))$p.value < alpha_normality
      }
    }, error = function(e) FALSE)
  }

  resolved_method <- if (method == "auto") {
    .auto_method(vx_c, vy_c, norm_x, norm_y, var_violated, alpha_normality)
  } else method

  result <- tryCatch({
    switch(resolved_method,
      pearson = {
        args <- c(list(x = vx_c, y = vy_c, method = "pearson"), cor_args)
        do.call(stats::cor.test, args)
      },
      spearman = {
        # cor.test does not compute CIs for Spearman; use Fisher z on rho
        args   <- c(list(x = vx_c, y = vy_c, method = "spearman",
                         exact = FALSE), cor_args)
        ct     <- do.call(stats::cor.test, args)
        n_obs  <- length(vx_c)
        rho    <- unname(ct$estimate)
        # Fisher z-transformation CI for Spearman rho
        if (!is.na(rho) && abs(rho) < 1 && n_obs > 3) {
          z     <- atanh(rho)
          se_z  <- 1 / sqrt(n_obs - 3)
          z_crit <- stats::qnorm(1 - (1 - 0.95) / 2)
          ci_lo <- tanh(z - z_crit * se_z)
          ci_hi <- tanh(z + z_crit * se_z)
          ct$conf.int <- c(ci_lo, ci_hi)
        } else {
          ct$conf.int <- c(NA_real_, NA_real_)
        }
        ct
      },
      kendall = {
        # cor.test does not compute CIs for Kendall; use Fisher z on tau
        args  <- c(list(x = vx_c, y = vy_c, method = "kendall",
                        exact = FALSE), cor_args)
        ct    <- do.call(stats::cor.test, args)
        n_obs <- length(vx_c)
        tau   <- unname(ct$estimate)
        if (!is.na(tau) && abs(tau) < 1 && n_obs > 3) {
          z     <- atanh(tau)
          se_z  <- 1 / sqrt(n_obs - 3)
          z_crit <- stats::qnorm(1 - (1 - 0.95) / 2)
          ci_lo <- tanh(z - z_crit * se_z)
          ci_hi <- tanh(z + z_crit * se_z)
          ct$conf.int <- c(ci_lo, ci_hi)
        } else {
          ct$conf.int <- c(NA_real_, NA_real_)
        }
        ct
      },
      pointbiserial = {
        if (.detect_binary(vx_c)) {
          vx_c <- as.integer(factor(vx_c)) - 1L
        } else {
          vy_c <- as.integer(factor(vy_c)) - 1L
        }
        args <- c(list(x = vx_c, y = vy_c, method = "pearson"), cor_args)
        do.call(stats::cor.test, args)
      },
      phi = {
        vx_c <- as.integer(factor(vx_c)) - 1L
        vy_c <- as.integer(factor(vy_c)) - 1L
        args <- c(list(x = vx_c, y = vy_c, method = "pearson"), cor_args)
        do.call(stats::cor.test, args)
      },
      cramerv = {
        tbl <- table(vx_c, vy_c)
        cs  <- suppressWarnings(stats::chisq.test(tbl))
        n_t <- sum(tbl)
        k   <- min(nrow(tbl), ncol(tbl))
        v   <- sqrt(cs$statistic / (n_t * (k - 1L)))
        list(estimate  = unname(v),
             p.value   = cs$p.value,
             conf.int  = c(NA_real_, NA_real_),
             method    = "Cramer's V")
      },
      poly = {
        if (!requireNamespace("polycor", quietly = TRUE))
          stop("Package 'polycor' is required for method = 'poly'. ",
               "Install it with: install.packages('polycor')", call. = FALSE)
        both_ord <- .detect_ordinal(vx_c) && .detect_ordinal(vy_c)
        if (both_ord) {
          res <- polycor::polychor(vx_c, vy_c, std.err = TRUE)
          list(estimate = res$rho, p.value = NA_real_,
               conf.int = c(NA_real_, NA_real_), method = "Polychoric")
        } else {
          cont_v <- if (is.numeric(vx_c) && !.detect_ordinal(vx_c)) vx_c else vy_c
          ord_v  <- if (.detect_ordinal(vx_c)) vx_c else vy_c
          res    <- polycor::polyserial(cont_v, ord_v, std.err = TRUE)
          list(estimate = res$rho, p.value = NA_real_,
               conf.int = c(NA_real_, NA_real_), method = "Polyserial")
        }
      },
      stop("Unknown method: '", resolved_method, "'.", call. = FALSE)
    )
  }, error = function(e) {
    warning(sprintf("Correlation failed for '%s' vs '%s': %s",
                    xname, yname, conditionMessage(e)), call. = FALSE)
    list(estimate = NA_real_, p.value = NA_real_,
         conf.int = c(NA_real_, NA_real_), method = resolved_method)
  })

  est   <- if (is.list(result)) result$estimate   else result$estimate
  pval  <- if (is.list(result)) result$p.value     else result$p.value
  ci    <- if (is.list(result)) result$conf.int    else result$conf.int
  mused <- if (is.list(result)) result$method      else result$method

  if (is.null(ci) || all(is.na(ci))) ci <- c(NA_real_, NA_real_)

  list(
    estimate    = unname(est),
    p_value     = unname(pval),
    ci_lower    = ci[1L],
    ci_upper    = ci[2L],
    method_used = resolved_method,
    n           = n,
    sig         = !is.na(pval) && pval < sig_threshold,
    xname       = xname,
    yname       = yname
  )
}

# Main function ----------------------------------------------------------------

#' @title Perform Correlation Analysis
#'
#' @description
#' Computes pairwise or targeted correlation analysis between numeric, ordinal,
#' or binary variables in a data frame, tibble, or matrix. Supports automatic
#' method selection (Pearson, Spearman, Kendall, point-biserial, phi,
#' Cramer's V, polychoric, and polyserial) based on data characteristics.
#' Returns tidy output tables of correlation estimates, masked correlations,
#' confidence intervals, and p-values.
#'
#' @param x A data frame, tibble, or matrix whose columns are the variables
#'   to correlate. Column names may contain special characters.
#' @param y A single character string naming a column in \code{x} to use as
#'   the primary variable. When both \code{y} and \code{z} are \code{NULL},
#'   all pairwise correlations among numeric columns of \code{x} are computed.
#'   When only \code{y} is supplied, correlations of \code{y} against all other
#'   columns are returned.
#' @param z A single character string naming a second column in \code{x}.
#'   When supplied together with \code{y}, only the correlation between
#'   \code{y} and \code{z} is computed.
#' @param metadata A data frame with the same number of rows as \code{x}
#'   containing sample-level metadata. Appended to the returned tidy tables
#'   when provided but not used in correlation computation.
#' @param remove A named list of one-sided formulas of the form
#'   \code{list("colname" ~ c("level1", "level2"))} specifying levels of a
#'   column in \code{x} to exclude before analysis. The left-hand side must be
#'   a single column name; the right-hand side may be a character vector of
#'   one or more levels.
#' @param method A single character string specifying the correlation method.
#'   One of \code{"auto"} (default), \code{"pearson"}, \code{"spearman"},
#'   \code{"kendall"}, \code{"pointbiserial"}, \code{"phi"}, \code{"cramerv"},
#'   or \code{"poly"}. When \code{"auto"}, the method is selected per variable
#'   pair based on data type, normality, and sample size.
#' @param force_test A list of two-element lists, each specifying a variable
#'   pair and the correlation method to use for that pair, overriding the
#'   global \code{method} argument. Each element must be of the form
#'   \code{list(c("var1", "var2"), "method")}. For example:
#'   \preformatted{force_test = list(
#'     list(c("var1", "var2"), "pearson"),
#'     list(c("var1", "var3"), "spearman")
#'   )}
#'   Valid methods are all values accepted by \code{method} except
#'   \code{"auto"}. Column names are order-insensitive within each pair.
#'   If a pair appears more than once, the last entry takes precedence with
#'   a warning. Pairs not listed here are analyzed using the global
#'   \code{method}. Default \code{NULL} (no overrides).
#' @param sig_threshold A single numeric value in \code{(0, 1)} specifying
#'   the significance threshold for p-values. Correlations with
#'   \code{p >= sig_threshold} are masked in the masked table. Default
#'   \code{0.05}.
#' @param normality_method A single character string, one of \code{"auto"}
#'   (default), \code{"shapiro"}, or \code{"lilliefors"}, controlling the
#'   normality test used in automatic method selection. \code{"auto"} applies
#'   Shapiro-Wilk for \eqn{n < 50} and Lilliefors otherwise.
#' @param alpha_normality A single numeric value in \code{(0, 1)} specifying
#'   the significance threshold used internally for normality and variance
#'   tests when \code{method = "auto"}. Default \code{0.05}.
#' @param verbose Logical. When \code{TRUE}, diagnostic messages about
#'   assumption checks and method selection are printed. Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[stats]{cor.test}}.
#'
#' @details
#' \strong{Automatic method selection} (\code{method = "auto"}):
#'
#' \tabular{ll}{
#'   \strong{Method} \tab \strong{Conditions} \cr
#'   Pearson         \tab Both continuous, normal, homoscedastic \cr
#'   Spearman        \tab Ordinal or continuous, non-normal or heteroscedastic,
#'                        \eqn{n \ge 30}, few ties \cr
#'   Kendall         \tab Ordinal or continuous, non-normal, \eqn{n < 30} or
#'                        many ties \cr
#'   Point-biserial  \tab One continuous, one binary variable \cr
#'   Phi             \tab Both binary \cr
#'   Cramer's V      \tab Both nominal (multi-category) factors \cr
#'   Polychoric      \tab Both ordered factors (requires \pkg{polycor}) \cr
#'   Polyserial      \tab One ordered factor, one continuous
#'                        (requires \pkg{polycor}) \cr
#' }
#'
#' Normality is assessed via Shapiro-Wilk (\eqn{n < 50}) or Lilliefors
#' (\eqn{n \ge 50}; requires \pkg{nortest}). Homoscedasticity is assessed via
#' Levene's test (requires \pkg{car}) with Bartlett's test as fallback.
#'
#' Confidence intervals for Cramer's V and polychoric/polyserial correlations
#' are not available and are returned as \code{NA}.
#'
#' @return A named list of class \code{"pondeR_correl"} with the following
#'   elements:
#'   \describe{
#'     \item{\code{correlations}}{A tidy data frame of pairwise correlation
#'       estimates with columns \code{var1}, \code{var2}, \code{estimate},
#'       \code{method_used}, and \code{n}.}
#'     \item{\code{correlations_masked}}{Same as \code{correlations} but
#'       \code{estimate} is \code{NA} for pairs where
#'       \code{p_value >= sig_threshold}.}
#'     \item{\code{p_values}}{A tidy data frame with columns \code{var1},
#'       \code{var2}, and \code{p_value}.}
#'     \item{\code{confidence_intervals}}{A tidy data frame with columns
#'       \code{var1}, \code{var2}, \code{ci_lower}, and \code{ci_upper}.}
#'     \item{\code{sig_threshold}}{The significance threshold used.}
#'     \item{\code{method}}{The method argument as supplied.}
#'     \item{\code{call}}{The matched call.}
#'   }
#'
#' @examples
#' ## Basic pairwise correlations (all numeric columns)
#' run_correl(mtcars)
#'
#' ## Correlate one variable against all others
#' run_correl(mtcars, y = "mpg")
#'
#' ## Single pair
#' run_correl(mtcars, y = "mpg", z = "wt")
#'
#' ## Specify method explicitly
#' run_correl(mtcars, method = "spearman")
#'
#' ## Exclude specific levels before analysis
#' df <- mtcars
#' df$cyl <- as.numeric(factor(df$cyl))
#' run_correl(df, remove = list("cyl" ~ c("4", "6")))
#'
#' ## Access individual output tables
#' res <- run_correl(mtcars)
#' res$correlations
#' res$p_values
#' res$correlations_masked
#'
#' @author John Lennon L. Calorio
#' @export
run_correl <- function(x,
                        y                = NULL,
                        z                = NULL,
                        metadata         = NULL,
                        remove           = NULL,
                        method           = "auto",
                        force_test       = NULL,
                        sig_threshold    = 0.05,
                        normality_method = "auto",
                        alpha_normality  = 0.05,
                        verbose          = FALSE,
                        ...) {

  mc <- match.call()

  # --- Input validation -------------------------------------------------------

  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("`x` must be a data frame, tibble, or matrix.",
         "\n  Supplied: ", class(x)[1L], call. = FALSE)
  }
  if (is.matrix(x)) x <- as.data.frame(x)
  if (ncol(x) < 1L) stop("`x` must have at least one column.", call. = FALSE)
  if (nrow(x) < 3L) stop("`x` must have at least 3 rows.", call. = FALSE)

  if (!is.null(y)) {
    if (!is.character(y) || length(y) != 1L)
      stop("`y` must be a single character string.", call. = FALSE)
    if (!y %in% names(x))
      stop("`y` ('", y, "') is not a column in `x`.", call. = FALSE)
  }
  if (!is.null(z)) {
    if (!is.character(z) || length(z) != 1L)
      stop("`z` must be a single character string.", call. = FALSE)
    if (!z %in% names(x))
      stop("`z` ('", z, "') is not a column in `x`.", call. = FALSE)
    if (is.null(y))
      stop("`z` cannot be specified without `y`.", call. = FALSE)
    if (identical(y, z))
      stop("`y` and `z` must be different column names.", call. = FALSE)
  }
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata))
      stop("`metadata` must be a data frame.", call. = FALSE)
    if (nrow(metadata) != nrow(x))
      stop("`metadata` must have the same number of rows as `x` (",
           nrow(x), " rows). Supplied: ", nrow(metadata), " rows.", call. = FALSE)
  }

  valid_methods <- c("auto", "pearson", "spearman", "kendall",
                     "pointbiserial", "phi", "cramerv", "poly")
  if (!is.character(method) || length(method) != 1L ||
      !method %in% valid_methods) {
    stop("`method` must be one of: ",
         paste0('"', valid_methods, '"', collapse = ", "), ".", call. = FALSE)
  }
  if (!is.numeric(sig_threshold) || length(sig_threshold) != 1L ||
      sig_threshold <= 0 || sig_threshold >= 1) {
    stop("`sig_threshold` must be a single numeric value between 0 and 1 ",
         "(exclusive). Supplied: ", sig_threshold, call. = FALSE)
  }
  valid_norm <- c("auto", "shapiro", "lilliefors")
  if (!normality_method %in% valid_norm)
    stop("`normality_method` must be one of: ",
         paste0('"', valid_norm, '"', collapse = ", "), ".", call. = FALSE)
  if (!is.numeric(alpha_normality) || length(alpha_normality) != 1L ||
      alpha_normality <= 0 || alpha_normality >= 1)
    stop("`alpha_normality` must be a single numeric value between 0 and 1 ",
         "(exclusive).", call. = FALSE)

  # --- Validate `force_test` --------------------------------------------------

  .ft_lookup <- list()  # will map "var1\rvar2" -> method string

  if (!is.null(force_test)) {
    if (!is.list(force_test))
      stop("`force_test` must be a list of two-element lists, e.g.:\n",
           '  list(list(c("var1", "var2"), "pearson"), list(c("var1", "var3"), "spearman"))',
           call. = FALSE)

    for (i in seq_along(force_test)) {
      item <- force_test[[i]]

      if (!is.list(item) || length(item) != 2L)
        stop("`force_test[[", i, "]]` must be a two-element list: ",
             'list(c("var1", "var2"), "method").', call. = FALSE)

      pair   <- item[[1L]]
      ft_mtd <- item[[2L]]

      if (!is.character(pair) || length(pair) != 2L)
        stop("The first element of `force_test[[", i, "]]` must be a ",
             "character vector of length 2 naming two columns.", call. = FALSE)
      if (!is.character(ft_mtd) || length(ft_mtd) != 1L ||
          !ft_mtd %in% valid_methods || ft_mtd == "auto")
        stop("The second element of `force_test[[", i, "]]` must be one of: ",
             paste0('"', setdiff(valid_methods, "auto"), '"', collapse = ", "),
             ".\n  Supplied in `force_test[[", i, "]]`: '", ft_mtd, "'.",
             call. = FALSE)
      if (!pair[1L] %in% names(x))
        stop("Column '", pair[1L], "' in `force_test[[", i, "]]` ",
             "not found in `x`.", call. = FALSE)
      if (!pair[2L] %in% names(x))
        stop("Column '", pair[2L], "' in `force_test[[", i, "]]` ",
             "not found in `x`.", call. = FALSE)

      # Store under a canonical sorted key so order doesn't matter
      key <- paste(sort(pair), collapse = "\r")
      if (key %in% names(.ft_lookup))
        warning(sprintf(
          "Duplicate `force_test` entry for pair ('%s', '%s'). ",
          pair[1L], pair[2L]
        ), "The last entry will be used.", call. = FALSE)
      .ft_lookup[[key]] <- ft_mtd
    }
  }  
  
  # --- Apply `remove` filters -------------------------------------------------

  if (!is.null(remove)) {
    if (!is.list(remove))
      stop("`remove` must be a list of formulas, e.g. ",
           'list("col" ~ c("level1", "level2")).', call. = FALSE)
    for (i in seq_along(remove)) {
      item <- remove[[i]]
      if (!inherits(item, "formula"))
        stop("`remove[[", i, "]]` must be a formula of the form ",
             '"colname" ~ c("level1", "level2").', call. = FALSE)
      lhs <- as.character(item[[2L]])
      rhs <- eval(item[[3L]])
      if (length(lhs) != 1L)
        stop("The left-hand side of `remove[[", i, "]]` must be a single ",
             "column name.", call. = FALSE)
      if (!lhs %in% names(x))
        stop("Column '", lhs, "' in `remove[[", i, "]]` not found in `x`.",
             call. = FALSE)
      keep <- !(as.character(x[[lhs]]) %in% as.character(rhs))
      x    <- x[keep, , drop = FALSE]
      if (!is.null(metadata)) metadata <- metadata[keep, , drop = FALSE]
      if (verbose)
        message(sprintf("Removed %d rows where `%s` %%in%% {%s}.",
                        sum(!keep), lhs,
                        paste0('"', rhs, '"', collapse = ", ")))
    }
    if (nrow(x) < 3L)
      stop("After applying `remove`, fewer than 3 rows remain. ",
           "Correlation analysis cannot proceed.", call. = FALSE)
  }

  # --- Determine variable pairs -----------------------------------------------

  extra_args <- list(...)
  col_names  <- names(x)

  pairs <- if (!is.null(y) && !is.null(z)) {
    list(c(y, z))
  } else if (!is.null(y)) {
    others <- setdiff(col_names, y)
    lapply(others, function(v) c(y, v))
  } else {
    nc    <- length(col_names)
    combs <- utils::combn(col_names, 2L, simplify = FALSE)
    combs
  }

  if (length(pairs) == 0L)
    stop("No variable pairs to correlate. ",
         "Ensure `x` has at least two columns.", call. = FALSE)

  # --- Compute correlations ---------------------------------------------------

  results <- vector("list", length(pairs))

  # for (i in seq_along(pairs)) {
  #   vnames     <- pairs[[i]]
  #   xname      <- vnames[1L]
  #   yname      <- vnames[2L]
  #   vx         <- x[[xname]]
  #   vy         <- x[[yname]]

  #   results[[i]] <- .correl_pair(
  #     vx               = vx,
  #     vy               = vy,
  #     xname            = xname,
  #     yname            = yname,
  #     method           = method,
  #     sig_threshold    = sig_threshold,
  #     normality_method = normality_method,
  #     alpha_normality  = alpha_normality,
  #     cor_args         = extra_args
  #   )

  #   if (verbose)
  #     message(sprintf("[%d/%d] %s vs %s — method: %s, r = %.4f, p = %.4f",
  #                     i, length(pairs), xname, yname,
  #                     results[[i]]$method_used,
  #                     results[[i]]$estimate %||% NA,
  #                     results[[i]]$p_value   %||% NA))
  # }
  for (i in seq_along(pairs)) {
    vnames     <- pairs[[i]]
    xname      <- vnames[1L]
    yname      <- vnames[2L]
    vx         <- x[[xname]]
    vy         <- x[[yname]]

    # Resolve per-pair forced method (if any), else fall back to `method`
    ft_key      <- paste(sort(c(xname, yname)), collapse = "\r")
    pair_method <- if (!is.null(.ft_lookup[[ft_key]])) {
      if (verbose)
        message(sprintf("[%d/%d] %s vs %s — force_test override: '%s'",
                        i, length(pairs), xname, yname, .ft_lookup[[ft_key]]))
      .ft_lookup[[ft_key]]
    } else {
      method
    }

    results[[i]] <- .correl_pair(
      vx               = vx,
      vy               = vy,
      xname            = xname,
      yname            = yname,
      method           = pair_method,
      sig_threshold    = sig_threshold,
      normality_method = normality_method,
      alpha_normality  = alpha_normality,
      cor_args         = extra_args
    )

    if (verbose)
      message(sprintf("[%d/%d] %s vs %s — method: %s, r = %.4f, p = %.4f",
                      i, length(pairs), xname, yname,
                      results[[i]]$method_used,
                      results[[i]]$estimate %||% NA,
                      results[[i]]$p_value   %||% NA))
  }

  # --- Assemble tidy tables ---------------------------------------------------

  var1   <- vapply(results, `[[`, character(1L), "xname")
  var2   <- vapply(results, `[[`, character(1L), "yname")
  est    <- vapply(results, `[[`, numeric(1L),   "estimate")
  pval   <- vapply(results, `[[`, numeric(1L),   "p_value")
  ci_lo  <- vapply(results, `[[`, numeric(1L),   "ci_lower")
  ci_hi  <- vapply(results, `[[`, numeric(1L),   "ci_upper")
  mused  <- vapply(results, `[[`, character(1L), "method_used")
  n_used <- vapply(results, `[[`, integer(1L),   "n")
  sig    <- vapply(results, `[[`, logical(1L),   "sig")

  correlations <- data.frame(
    var1        = var1,
    var2        = var2,
    estimate    = est,
    method_used = mused,
    n           = n_used,
    stringsAsFactors = FALSE
  )

  est_masked             <- est
  est_masked[!sig]       <- NA_real_

  correlations_masked <- data.frame(
    var1        = var1,
    var2        = var2,
    estimate    = est_masked,
    method_used = mused,
    n           = n_used,
    stringsAsFactors = FALSE
  )

  p_values <- data.frame(
    var1    = var1,
    var2    = var2,
    p_value = pval,
    stringsAsFactors = FALSE
  )

  confidence_intervals <- data.frame(
    var1     = var1,
    var2     = var2,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    stringsAsFactors = FALSE
  )

  out <- list(
    correlations         = correlations,
    correlations_masked  = correlations_masked,
    p_values             = p_values,
    confidence_intervals = confidence_intervals,
    sig_threshold        = sig_threshold,
    method               = method,
    call                 = mc
  )
  class(out) <- c("pondeR_correl", "list")
  out
}

# Null-coalescing operator (internal) -----------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b