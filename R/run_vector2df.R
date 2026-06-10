#' Convert a Named Vector to a Two-Column Data Frame
#'
#' Transforms a named vector into a tidy two-column data frame, where one
#' column holds the names and the other holds the corresponding values.
#' Optionally filters rows by a value threshold or range. Useful for quickly
#' tabulating named outputs such as VIP scores, loadings, or any named
#' numeric/character vector.
#'
#' @param x A named vector (numeric, integer, character, or logical). Must have
#'   names; an error is raised if \code{names(x)} is \code{NULL}.
#' @param name A single character string giving the name of the column that
#'   will contain the vector names. Defaults to \code{"Variable"}.
#' @param value A single character string giving the name of the column that
#'   will contain the vector values. Defaults to \code{"Value"}.
#' @param sort_by A single character string controlling row ordering. One of
#'   \code{"none"} (original order), \code{"name"} (alphabetical by name
#'   column), \code{"value_asc"} (ascending by value), or
#'   \code{"value_desc"} (descending by value). Defaults to \code{"none"}.
#' @param min A single numeric value used as a lower threshold for filtering
#'   the value column. Behaviour depends on \code{range_type}:
#'   \itemize{
#'     \item \code{"inner"}: rows with \code{value >= min} are kept (when
#'       \code{max} is \code{NULL}), or rows satisfying
#'       \code{value >= min & value <= max} (when both are supplied).
#'     \item \code{"outer"}: rows with \code{value <= min} are kept (when
#'       \code{max} is \code{NULL}), or rows satisfying
#'       \code{value <= min | value >= max} (when both are supplied).
#'   }
#'   Defaults to \code{NULL} (no lower filtering).
#' @param max A single numeric value used as an upper threshold for filtering
#'   the value column. Behaviour depends on \code{range_type}:
#'   \itemize{
#'     \item \code{"inner"}: rows with \code{value <= max} are kept (when
#'       \code{min} is \code{NULL}), or rows satisfying
#'       \code{value >= min & value <= max} (when both are supplied).
#'     \item \code{"outer"}: rows with \code{value >= max} are kept (when
#'       \code{min} is \code{NULL}), or rows satisfying
#'       \code{value <= min | value >= max} (when both are supplied).
#'   }
#'   Defaults to \code{NULL} (no upper filtering).
#' @param range_type A single character string, either \code{"inner"} or
#'   \code{"outer"}, that controls how \code{min} and \code{max} are applied
#'   when both are supplied. \code{"inner"} keeps values \emph{between} the
#'   thresholds (inclusive); \code{"outer"} keeps values \emph{outside} the
#'   thresholds (inclusive). Ignored when only one of \code{min}/\code{max}
#'   is supplied. Defaults to \code{"inner"}.
#' @param stringsAsFactors Logical. If \code{TRUE}, the name column is
#'   returned as a \code{factor}. Defaults to \code{FALSE}.
#'
#' @return A \code{data.frame} with two columns named according to \code{name}
#'   and \code{value}, and one row per (retained) element of \code{x}. Row
#'   names are reset to sequential integers.
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' vip <- c(`623.0149@9.699` = 0.3616770, `769.4391@9.699` = 1.0151602,
#'          `1147.7305@9.699` = 1.2068827, `397.3033@9.7028` = 0.5754963)
#' run_vector2df(vip, name = "Feature", value = "VIP",
#'               min = 1.0, sort_by = "value_desc")
#'
#' @export
run_vector2df <- function(x,
                          name             = "Variable",
                          value            = "Value",
                          sort_by          = c("none", "name", "value_asc", "value_desc"),
                          min              = NULL,
                          max              = NULL,
                          range_type       = c("inner", "outer"),
                          stringsAsFactors = FALSE) {

  ## --- input checks ----------------------------------------------------------
  if (is.null(names(x))) {
    stop("`x` must be a named vector (names(x) is NULL).", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || nchar(trimws(name)) == 0L) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(value) || length(value) != 1L || nchar(trimws(value)) == 0L) {
    stop("`value` must be a single non-empty character string.", call. = FALSE)
  }
  if (name == value) {
    stop("`name` and `value` must be different strings.", call. = FALSE)
  }
  if (!is.null(min) && (!is.numeric(min) || length(min) != 1L)) {
    stop("`min` must be a single numeric value or NULL.", call. = FALSE)
  }
  if (!is.null(max) && (!is.numeric(max) || length(max) != 1L)) {
    stop("`max` must be a single numeric value or NULL.", call. = FALSE)
  }
  if (!is.null(min) && !is.null(max) && min > max) {
    stop("`min` must be less than or equal to `max`.", call. = FALSE)
  }

  sort_by    <- match.arg(sort_by)
  range_type <- match.arg(range_type)

  ## --- build data frame ------------------------------------------------------
  out <- data.frame(
    names(x),
    unname(x),
    stringsAsFactors = stringsAsFactors,
    row.names        = NULL
  )
  colnames(out) <- c(name, value)

  ## --- filtering -------------------------------------------------------------
  v <- out[[value]]

  keep <- if (!is.null(min) && !is.null(max)) {
    if (range_type == "inner") v >= min & v <= max else v <= min | v >= max
  } else if (!is.null(min)) {
    if (range_type == "inner") v >= min else v <= min
  } else if (!is.null(max)) {
    if (range_type == "inner") v <= max else v >= max
  } else {
    rep(TRUE, nrow(out))
  }

  out <- out[keep, ]

  ## --- sorting ---------------------------------------------------------------
  out <- switch(
    sort_by,
    none       = out,
    name       = out[order(out[[name]]), ],
    value_asc  = out[order(out[[value]]), ],
    value_desc = out[order(out[[value]], decreasing = TRUE), ]
  )

  rownames(out) <- NULL
  out
}