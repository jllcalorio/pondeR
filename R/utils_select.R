#' Select Numeric Variables from a Data Frame
#'
#' @description
#' This function selects variables from a data frame based on their numeric
#' types (integer, double, logical, or complex). It also allows forcing
#' specific columns to be treated as numeric even if they are stored as
#' other types.
#'
#' @param x A data frame or matrix.
#' @param force_num Character vector of column names in \code{x} to be
#'   forcefully included as numeric variables.
#' @param include_int Logical. If \code{TRUE} (default), includes integer
#'   type variables.
#' @param include_dbl Logical. If \code{TRUE} (default), includes double
#'   type variables.
#' @param include_lgl Logical. If \code{TRUE} (default), includes logical
#'   type variables.
#' @param include_comp Logical. If \code{TRUE} (default), includes complex
#'   type variables.
#'
#' @return A data frame containing only the selected numeric variables.
#' @author John Lennon L. Calorio
#' @export
#'
#' @examples
#' df <- data.frame(a = 1:3, b = c(1.1, 2.2, 3.3), 
#' c = c("x", "y", "z"),
#' d = c(TRUE, FALSE, TRUE))
#' select_numvars(df)
select_numvars <- function(x,
                           force_num    = NULL,
                           include_int  = TRUE,
                           include_dbl  = TRUE,
                           include_lgl  = TRUE,
                           include_comp = TRUE) {
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) stop("'x' must be a data frame or matrix.")

  check_num <- function(col, name) {
    if (!is.null(force_num) && name %in% force_num) return(TRUE)
    if (include_int  && is.integer(col)) return(TRUE)
    if (include_dbl  && is.double(col) && !is.factor(col)) return(TRUE)
    if (include_lgl  && is.logical(col)) return(TRUE)
    if (include_comp && is.complex(col)) return(TRUE)
    return(FALSE)
  }

  keep <- vapply(names(x), function(n) check_num(x[[n]], n), logical(1))
  x[, keep, drop = FALSE]
}

#' Select Categorical Variables from a Data Frame
#'
#' @description
#' This function selects variables from a data frame that are categorical
#' (factors, ordered factors, characters, or logicals). It also allows
#' forcing specific columns to be treated as categorical even if they are
#' stored as numeric types.
#'
#' @param x A data frame or matrix.
#' @param force_cat Character vector of column names in \code{x} to be
#'   forcefully included as categorical variables.
#' @param include_chr Logical. If \code{TRUE} (default), includes character
#'   type variables.
#' @param include_lgl Logical. If \code{FALSE} (default), excludes logical
#'   type variables.
#'
#' @return A data frame containing only the selected categorical variables.
#' @author John Lennon L. Calorio
#' @export
#'
#' @examples
#' df <- data.frame(a = factor(c("A", "B")), b = 1:2, 
#' c = c("x", "y"),
#' d = c(TRUE, FALSE))
#' select_catvars(df)
select_catvars <- function(x,
                           force_cat   = NULL,
                           include_chr = TRUE,
                           include_lgl = FALSE) {
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) stop("'x' must be a data frame or matrix.")

  check_cat <- function(col, name) {
    if (!is.null(force_cat) && name %in% force_cat) return(TRUE)
    if (is.factor(col)) return(TRUE)
    if (include_chr && is.character(col)) return(TRUE)
    if (include_lgl && is.logical(col))   return(TRUE)
    return(FALSE)
  }

  keep <- vapply(names(x), function(n) check_cat(x[[n]], n), logical(1))
  x[, keep, drop = FALSE]
}