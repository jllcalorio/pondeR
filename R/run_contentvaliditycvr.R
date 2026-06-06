#' Run Lawshe's Content Validity Ratio (CVR) Analysis
#'
#' @description
#' Computes Lawshe's (1975) Content Validity Ratio (CVR) for one or more items
#' rated by an expert panel. Each expert classifies an item as "essential",
#' "useful but not essential", or "not necessary"; CVR captures the proportion
#' of experts endorsing "essential" relative to chance. The function returns a
#' tidy data frame of CVR values, critical values, and significance
#' interpretations for each item.
#'
#' @param x A numeric data frame or matrix where rows represent experts and
#'   columns represent items.
#' @param essential A numeric vector indicating the rating value(s) that
#'   classify an item as "essential". Defaults to \code{1}. A vector may be
#'   supplied (e.g., \code{c(1, 2)}) when multiple adjacent categories should
#'   count as essential.
#' @param n An integer specifying the panel size used to compute the minimum
#'   CVR considered statistically significant via \code{cvr_critical()}.
#'   Defaults to \code{dim(x)[1]} (the number of rows in \code{x}).
#' @param alpha A numeric value for the one-tailed significance level used in
#'   \code{cvr_critical()}. Defaults to \code{0.05}.
#' @param na.rm Logical. If \code{TRUE}, missing ratings are excluded when
#'   counting experts per item. Defaults to \code{FALSE}.
#'
#' @details
#' CVR is defined as:
#' \deqn{\text{CVR} = \frac{n_e - n/2}{n/2}}
#' where \eqn{n_e} is the number of experts rating the item as "essential" and
#' \eqn{n} is the total number of experts. CVR ranges from \eqn{-1} (no expert
#' endorses "essential") to \eqn{+1} (all experts endorse "essential"). A value
#' of \eqn{0} indicates that exactly half the panel rated the item "essential".
#'
#' Statistical significance of each CVR is evaluated against the critical value
#' returned by \code{contentValidity::cvr_critical(n, alpha)}, which implements
#' the recalculated critical values of Wilson, Pan, and Schumsky (2012).
#'
#' @return A data frame of class \code{run_contentvaliditycvr} with one row
#' per item and the following columns:
#' \describe{
#'   \item{\code{CVR}}{Lawshe's Content Validity Ratio.}
#'   \item{\code{Critical Value}}{Minimum CVR for statistical significance at
#'     the specified \code{alpha} and panel size \code{n}.}
#'   \item{\code{Interpretation}}{Character; \code{"Significant"} if
#'     \code{CVR >= Critical Value}, otherwise \code{"Not significant"}.}
#' }
#'
#' @references
#' Lawshe, C. H. (1975). A quantitative approach to content validity.
#' \emph{Personnel Psychology}, \emph{28}(4), 563--575.
#' \doi{10.1111/j.1744-6570.1975.tb01393.x}
#'
#' Wilson, F. R., Pan, W., & Schumsky, D. A. (2012). Recalculation of the
#' critical values for Lawshe's content validity ratio. \emph{Measurement and
#' Evaluation in Counseling and Development}, \emph{45}(3), 197--210.
#' \doi{10.1177/0748175612440286}
#'
#' @examples
#' \dontrun{
#' # Lawshe's (1975) original 3-point scale:
#' # 1 = essential, 2 = useful but not essential, 3 = not necessary
#' ratings <- data.frame(
#'   item1 = c(1, 1, 1, 2, 1),
#'   item2 = c(3, 3, 1, 3, 3),
#'   item3 = c(1, 1, 1, 1, 1)
#' )
#' result <- run_contentvaliditycvr(ratings, essential = 1)
#' result
#'
#' # Treat both "essential" and "useful but not essential" as endorsing the item
#' result2 <- run_contentvaliditycvr(ratings, essential = c(1, 2))
#' result2
#' }
#'
#' @seealso
#' \code{\link[contentValidity]{cvr}} for the underlying CVR computation,
#' \code{\link[contentValidity]{cvr_critical}} for the critical value lookup.
#'
#' @export
run_contentvaliditycvr <- function(
    x,
    essential = 1,
    n         = dim(x)[1L],
    alpha     = 0.05,
    na.rm     = FALSE
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
      "'x' contains missing values. Set 'na.rm = TRUE' to exclude them, ",
      "or inspect your data before proceeding.",
      call. = FALSE
    )
  }

  if (!is.numeric(essential) || length(essential) < 1L) {
    stop(
      "'essential' must be a non-empty numeric vector of rating value(s).",
      call. = FALSE
    )
  }

  observed_vals <- x_mat[!is.na(x_mat)]
  if (!any(observed_vals %in% essential)) {
    stop(
      "None of the 'essential' value(s) (", paste(essential, collapse = ", "),
      ") were found in 'x'. ",
      "Verify that 'essential' matches the coding used in your data.",
      call. = FALSE
    )
  }

  if (!is.numeric(n) || length(n) != 1L || n < 1L) {
    stop("'n' must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value strictly between 0 and 1.",
         call. = FALSE)
  }

  if (!is.logical(na.rm) || length(na.rm) != 1L) {
    stop("'na.rm' must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  # --- computation ------------------------------------------------------------
  cvr_vals     <- contentValidity::cvr(ratings = x, essential = essential, na.rm = na.rm)
  cvr_critical <- contentValidity::cvr_critical(n_experts = n, alpha = alpha)

  # cvr() returns a named numeric vector; coerce to data frame
  items_out <- data.frame(
    CVR              = as.numeric(cvr_vals),
    `Critical Value` = cvr_critical,
    check.names      = FALSE,
    row.names        = names(cvr_vals)
  )

  items_out[["Interpretation"]] <- ifelse(
    items_out[["CVR"]] >= items_out[["Critical Value"]],
    "Significant",
    "Not significant"
  )

  # --- output -----------------------------------------------------------------
  class(items_out) <- c("run_contentvaliditycvr", "data.frame")
  attr(items_out, "n")     <- n
  attr(items_out, "alpha") <- alpha
  items_out
}


# --- S3 methods ---------------------------------------------------------------

#' @export
print.run_contentvaliditycvr <- function(x, digits = 3L, ...) {
  cat("Content Validity Ratio (CVR) Analysis\n")
  cat(sprintf(
    "Panel size: %d  |  alpha = %.3f  |  Critical value: %.3f\n",
    attr(x, "n"),
    attr(x, "alpha"),
    x[["Critical Value"]][1L]
  ))
  cat(strrep("-", 50L), "\n\n")

  print_df           <- as.data.frame(x)
  num_cols           <- vapply(print_df, is.numeric, logical(1L))
  print_df[num_cols] <- lapply(print_df[num_cols], round, digits = digits)
  print(print_df, row.names = TRUE)

  invisible(x)
}


#' @export
summary.run_contentvaliditycvr <- function(object, digits = 3L, ...) {
  cat("Content Validity Ratio (CVR) Analysis — Summary\n")
  cat(strrep("=", 50L), "\n\n")
  cat("Panel size  :", attr(object, "n"),    "\n")
  cat("Alpha       :", attr(object, "alpha"), "\n")
  cat("Critical CVR:", round(object[["Critical Value"]][1L], digits), "\n\n")

  n_sig     <- sum(object[["Interpretation"]] == "Significant",     na.rm = TRUE)
  n_not_sig <- sum(object[["Interpretation"]] == "Not significant", na.rm = TRUE)
  cat(sprintf("Items significant    : %d / %d\n", n_sig, nrow(object)))
  cat(sprintf("Items not significant: %d / %d\n\n", n_not_sig, nrow(object)))

  print_df           <- as.data.frame(object)
  num_cols           <- vapply(print_df, is.numeric, logical(1L))
  print_df[num_cols] <- lapply(print_df[num_cols], round, digits = digits)
  print(print_df, row.names = TRUE)

  invisible(object)
}