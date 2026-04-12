# =============================================================================
# pondeR - utility checkers
# is_tabular | is_list | has_names | has_na | has_zero | has_negative | has_inf
# =============================================================================


# -----------------------------------------------------------------------------
# is_tabular
# -----------------------------------------------------------------------------

#' @title Check if an Object is a Tabular Data Structure
#'
#' @description
#' Tests whether an object belongs to a recognized tabular data class -
#' \code{data.frame}, \code{tbl_df} (tibble), \code{data.table}, or
#' \code{matrix} - and confirms that it has exactly two dimensions. Returns a
#' single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test.
#' @param warn Logical scalar. If \code{TRUE} (default), informative warnings
#'   are emitted when \code{x} fails the check, describing the specific reason
#'   for failure. Set to \code{FALSE} to suppress all warnings and obtain a
#'   silent \code{FALSE}.
#'
#' @details
#' An object is considered tabular if it satisfies **both** of the following
#' conditions:
#'
#' \enumerate{
#'   \item Its class is one of \code{"data.frame"}, \code{"tbl_df"},
#'         \code{"tbl"} (generic tibble superclass), \code{"data.table"}, or
#'         \code{"matrix"}.
#'   \item It has exactly two dimensions (i.e., \code{length(dim(x)) == 2}).
#' }
#'
#' The class check uses \code{inherits()}, so subclasses of \code{data.frame}
#' (e.g., tibbles, which inherit from \code{"data.frame"}) are naturally
#' detected. A \code{matrix} always has two dimensions by construction, but the
#' dimension check is applied uniformly for consistency and to guard against
#' unusual edge cases such as zero-dimensional matrices produced by certain
#' external packages.
#'
#' \code{NULL} and objects with no class attribute are handled gracefully:
#' they always return \code{FALSE}.
#'
#' Note that \code{is_tabular()} does **not** check whether the object is
#' non-empty. A zero-row or zero-column tabular object will still return
#' \code{TRUE}, as it retains the correct class and dimensions. Callers that
#' require non-empty data should perform that check separately.
#'
#' @return A single logical value: \code{TRUE} if \code{x} is a recognized
#'   tabular structure with two dimensions, \code{FALSE} otherwise. The return
#'   value is never \code{NA}.
#'
#' @examples
#' # --- TRUE cases ---
#'
#' # Base data.frame
#' is_tabular(mtcars)
#'
#' # Matrix
#' is_tabular(matrix(1:12, nrow = 3))
#'
#' # tibble (if tibble is installed)
#' if (requireNamespace("tibble", quietly = TRUE)) {
#'   is_tabular(tibble::tibble(x = 1:3, y = letters[1:3]))
#' }
#'
#' # data.table (if data.table is installed)
#' if (requireNamespace("data.table", quietly = TRUE)) {
#'   is_tabular(data.table::data.table(a = 1:5, b = rnorm(5)))
#' }
#'
#' # Zero-row data.frame - still tabular
#' is_tabular(mtcars[0, ])
#'
#' # --- FALSE cases ---
#'
#' # Atomic vector
#' is_tabular(1:10)
#'
#' # List
#' is_tabular(list(a = 1, b = 2))
#'
#' # 3-D array (more than two dimensions)
#' is_tabular(array(1:24, dim = c(2, 3, 4)))
#'
#' # NULL
#' is_tabular(NULL)
#'
#' # Suppress warnings for programmatic use
#' is_tabular("not a table", warn = FALSE)
#'
#' @author John Lennon L. Calorio
#' @export
is_tabular <- function(x, warn = TRUE) {

  # --- Input validation for `warn` itself -----------------------------------
  .check_warn(warn)

  # --- NULL guard -----------------------------------------------------------
  if (is.null(x)) {
    if (warn) {
      warning(
        "'x' is NULL. A tabular object (data.frame, tibble, data.table, or ",
        "matrix) was expected.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  # --- Recognized tabular classes -------------------------------------------
  tabular_classes <- c("data.frame", "tbl_df", "tbl", "data.table", "matrix")
  has_tabular_class <- inherits(x, tabular_classes)

  if (!has_tabular_class) {
    if (warn) {
      actual_class <- paste(class(x), collapse = ", ")
      warning(
        "'x' has class <", actual_class, ">, which is not a recognized ",
        "tabular type. Expected one of: data.frame, tibble (tbl_df), ",
        "data.table, or matrix.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  # --- Dimension check ------------------------------------------------------
  dims <- dim(x)

  if (is.null(dims)) {
    # Should be rare for the classes above, but guard defensively.
    if (warn) {
      warning(
        "'x' inherits from a tabular class but has no dim() attribute. ",
        "This is unexpected; the object may be malformed.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  if (length(dims) != 2L) {
    if (warn) {
      warning(
        "'x' has ", length(dims), " dimension(s) (dim = ",
        paste(dims, collapse = " x "), "), but a tabular object must have ",
        "exactly 2 dimensions.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  TRUE
}

# -----------------------------------------------------------------------------
# is_list
# -----------------------------------------------------------------------------

#' @title Check if an Object is a Plain List
#'
#' @description
#' Tests whether an object is a base R \code{list} (but not a \code{data.frame},
#' tibble, or \code{data.table}, which are technically lists under the hood).
#' Returns a single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test.
#' @param warn Logical scalar. If \code{TRUE} (default), an informative warning
#'   is emitted when \code{x} fails the check. Set to \code{FALSE} to suppress
#'   all warnings and obtain a silent \code{FALSE}.
#'
#' @details
#' \code{is_list()} returns \code{TRUE} only when \code{x} is a plain
#' \code{list} - i.e., \code{is.list(x)} is \code{TRUE} \emph{and} \code{x}
#' does not inherit from any tabular class (\code{data.frame}, \code{tbl_df},
#' \code{tbl}, \code{data.table}, or \code{matrix}). This distinction matters
#' because \code{data.frame} and \code{data.table} objects are internally lists
#' but should not be treated as such in the context of \pkg{pondeR} workflows.
#'
#' \code{NULL} always returns \code{FALSE}.
#'
#' @return A single logical value: \code{TRUE} if \code{x} is a plain list,
#'   \code{FALSE} otherwise. Never \code{NA}.
#'
#' @examples
#' # TRUE cases
#' is_list(list(a = 1, b = "hello"))
#' is_list(list(1:3, TRUE, NULL))
#'
#' # FALSE cases - tabular objects are excluded
#' is_list(mtcars)
#' is_list(matrix(1:4, 2))
#'
#' # FALSE - atomic vector
#' is_list(1:10)
#'
#' # FALSE - NULL
#' is_list(NULL)
#'
#' # Silent FALSE for programmatic use
#' is_list(42, warn = FALSE)
#'
#' @author John Lennon L. Calorio
#' @export
is_list <- function(x, warn = TRUE) {

  .check_warn(warn)

  if (is.null(x)) {
    if (warn) warning("'x' is NULL. A plain list was expected.", call. = FALSE)
    return(FALSE)
  }

  # Tabular objects are internally lists - exclude them explicitly.
  if (is_tabular(x, warn = FALSE)) {
    if (warn) {
      warning(
        "'x' has class <", paste(class(x), collapse = ", "), ">, which is a ",
        "tabular object, not a plain list. Use is_tabular() to check for ",
        "data.frame, tibble, data.table, or matrix objects.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  if (!is.list(x)) {
    if (warn) {
      warning(
        "'x' has class <", paste(class(x), collapse = ", "), ">, which is not ",
        "a list. A plain list object was expected.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  TRUE
}


# -----------------------------------------------------------------------------
# has_names
# -----------------------------------------------------------------------------

#' @title Check if an Object has Names
#'
#' @description
#' Tests whether an object has names. For tabular objects (\code{data.frame},
#' tibble, \code{data.table}, \code{matrix}), column names are checked. For
#' plain lists, element names are checked. For atomic vectors, element names
#' are checked. Returns a single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test. Supported types: \code{data.frame}, tibble,
#'   \code{data.table}, \code{matrix}, plain \code{list}, or atomic vector.
#' @param warn Logical scalar. If \code{TRUE} (default), informative warnings
#'   are emitted describing why \code{x} fails the check.
#'
#' @details
#' An object is considered named if \code{names()} (or \code{colnames()} for
#' tabular objects) returns a non-\code{NULL} character vector where every
#' element is a non-empty, non-\code{NA} string. A partially named object
#' (some names are \code{""} or \code{NA}) returns \code{FALSE}, and a warning
#' is emitted identifying the problematic positions.
#'
#' \code{NULL} always returns \code{FALSE}. Unsupported object types (e.g.,
#' environments, S4 objects) also return \code{FALSE} with a warning.
#'
#' @return A single logical value: \code{TRUE} if all names are present and
#'   non-empty, \code{FALSE} otherwise. Never \code{NA}.
#'
#' @examples
#' # Tabular - TRUE
#' has_names(mtcars)
#'
#' # Matrix - TRUE
#' has_names(matrix(1:4, 2, dimnames = list(NULL, c("A", "B"))))
#'
#' # Matrix without column names - FALSE
#' has_names(matrix(1:4, 2))
#'
#' # Named list - TRUE
#' has_names(list(a = 1, b = 2))
#'
#' # Unnamed list - FALSE
#' has_names(list(1, 2, 3))
#'
#' # Named vector - TRUE
#' has_names(c(x = 1, y = 2))
#'
#' # Partially named list - FALSE
#' has_names(list(a = 1, 2, b = 3))
#'
#' # NULL - FALSE
#' has_names(NULL)
#'
#' @author John Lennon L. Calorio
#' @export
has_names <- function(x, warn = TRUE) {

  .check_warn(warn)

  if (is.null(x)) {
    if (warn) warning("'x' is NULL.", call. = FALSE)
    return(FALSE)
  }

  # Route name extraction by type.
  if (is_tabular(x, warn = FALSE)) {
    nms <- colnames(x)
    context <- "column names"
  } else if (is_list(x, warn = FALSE) || is.atomic(x)) {
    nms <- names(x)
    context <- "names"
  } else {
    if (warn) {
      warning(
        "'x' has class <", paste(class(x), collapse = ", "), ">, which is not ",
        "a supported type. Supply a data.frame, tibble, data.table, matrix, ",
        "list, or atomic vector.",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  # No names attribute at all.
  if (is.null(nms)) {
    if (warn) {
      warning("'x' has no ", context, " attribute.", call. = FALSE)
    }
    return(FALSE)
  }

  # Detect empty-string or NA names.
  bad_idx <- which(is.na(nms) | nms == "")
  if (length(bad_idx) > 0L) {
    if (warn) {
      warning(
        "'x' has missing or empty ", context, " at position(s): ",
        paste(bad_idx, collapse = ", "), ".",
        call. = FALSE
      )
    }
    return(FALSE)
  }

  TRUE
}


# -----------------------------------------------------------------------------
# Internal helper - .check_warn
# Used by has_na / has_zero / has_negative / has_inf to validate 'warn'.
# Not exported.
# -----------------------------------------------------------------------------

.check_warn <- function(warn) {
  if (!is.logical(warn) || length(warn) != 1L || is.na(warn)) {
    stop("'warn' must be a single non-NA logical value (TRUE or FALSE).",
         call. = FALSE)
  }
}


# -----------------------------------------------------------------------------
# Internal helper - .coerce_to_atomic
# Extracts a flat atomic vector from supported input types.
# Returns NULL (invisibly) if x is unsupported, also emits warning if warn=TRUE.
# Not exported.
# -----------------------------------------------------------------------------

.coerce_to_atomic <- function(x, fn_name, warn) {
  if (is.null(x)) {
    if (warn) warning("'x' is NULL in ", fn_name, "().", call. = FALSE)
    return(NULL)
  }
  if (is_tabular(x, warn = FALSE)) {
    # Unlist column-wise; preserves NAs / Infs faithfully.
    return(unlist(x, use.names = FALSE))
  }
  if (is_list(x, warn = FALSE)) {
    return(unlist(x, use.names = FALSE))
  }
  if (is.atomic(x)) {
    return(as.vector(x))
  }
  if (warn) {
    warning(
      "'x' has class <", paste(class(x), collapse = ", "), ">. ",
      fn_name, "() supports data.frame, tibble, data.table, matrix, ",
      "list, or atomic vector.",
      call. = FALSE
    )
  }
  NULL
}


# -----------------------------------------------------------------------------
# has_na
# -----------------------------------------------------------------------------

#' @title Check if an Object Contains Missing Values
#'
#' @description
#' Tests whether an object contains at least one \code{NA} or \code{NaN} value.
#' Accepts tabular objects (\code{data.frame}, tibble, \code{data.table},
#' \code{matrix}), plain lists, or atomic vectors (including individual table
#' columns). Returns a single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test.
#' @param warn Logical scalar. If \code{TRUE} (default), a warning is emitted
#'   when \code{x} is \code{NULL} or of an unsupported type.
#'
#' @details
#' \code{NaN} is treated as a missing value (consistent with
#' \code{is.na(NaN) == TRUE} in base R). Zero (\code{0}) is \emph{not}
#' considered missing. For tabular inputs the entire table is scanned; for
#' vectors or list inputs all elements are scanned after unlisting.
#'
#' @return \code{TRUE} if at least one \code{NA} or \code{NaN} is found,
#'   \code{FALSE} otherwise. Never \code{NA}.
#'
#' @examples
#' has_na(c(1, 2, NA, 4))          # TRUE
#' has_na(c(1, 2, NaN, 4))         # TRUE
#' has_na(c(1, 2, 0, 4))           # FALSE - 0 is not NA
#' has_na(c(1, 2, 3))              # FALSE
#'
#' df <- data.frame(x = c(1, NA), y = c(3, 4))
#' has_na(df)                       # TRUE
#'
#' has_na(mtcars)                   # FALSE
#'
#' @author John Lennon L. Calorio
#' @export
has_na <- function(x, warn = TRUE) {
  .check_warn(warn)
  v <- .coerce_to_atomic(x, "has_na", warn)
  if (is.null(v)) return(FALSE)
  anyNA(v)
}


# -----------------------------------------------------------------------------
# has_zero
# -----------------------------------------------------------------------------

#' @title Check if an Object Contains Zero Values
#'
#' @description
#' Tests whether an object contains at least one element equal to \code{0}.
#' Accepts tabular objects, plain lists, or atomic vectors. Returns a single
#' logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test.
#' @param warn Logical scalar. If \code{TRUE} (default), a warning is emitted
#'   when \code{x} is \code{NULL} or of an unsupported type.
#'
#' @details
#' Only numeric (or coercible-to-numeric) elements are tested. Non-numeric
#' elements (e.g., character, logical) are silently skipped. \code{NA} values
#' are not counted as zero.
#'
#' @return \code{TRUE} if at least one \code{0} is found among numeric
#'   elements, \code{FALSE} otherwise. Never \code{NA}.
#'
#' @examples
#' has_zero(c(1, 0, 3))            # TRUE
#' has_zero(c(1, 2, 3))            # FALSE
#' has_zero(c(NA, 0, 1))           # TRUE - NA is not zero, but 0 is present
#' has_zero(c(NA, NA))             # FALSE
#'
#' df <- data.frame(x = c(1, 0), y = c(3, 4))
#' has_zero(df)                     # TRUE
#'
#' @author John Lennon L. Calorio
#' @export
has_zero <- function(x, warn = TRUE) {
  .check_warn(warn)
  v <- .coerce_to_atomic(x, "has_zero", warn)
  if (is.null(v)) return(FALSE)
  v_num <- suppressWarnings(as.numeric(v))
  any(v_num == 0, na.rm = TRUE)
}


# -----------------------------------------------------------------------------
# has_negative
# -----------------------------------------------------------------------------

#' @title Check if an Object Contains Negative Values
#'
#' @description
#' Tests whether an object contains at least one element with a value strictly
#' less than \code{0}. Accepts tabular objects, plain lists, or atomic vectors.
#' Returns a single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test.
#' @param warn Logical scalar. If \code{TRUE} (default), a warning is emitted
#'   when \code{x} is \code{NULL} or of an unsupported type.
#'
#' @details
#' Only numeric (or coercible-to-numeric) elements are tested. \code{NA}
#' values are excluded from the comparison. \code{-Inf} is treated as
#' negative.
#'
#' @return \code{TRUE} if at least one strictly negative value is found,
#'   \code{FALSE} otherwise. Never \code{NA}.
#'
#' @examples
#' has_negative(c(1, -2, 3))       # TRUE
#' has_negative(c(1, 2, 3))        # FALSE
#' has_negative(c(0, -0.001, NA))  # TRUE
#' has_negative(c(-Inf, 1, 2))     # TRUE
#'
#' df <- data.frame(x = c(1, -1), y = c(3, 4))
#' has_negative(df)                 # TRUE
#'
#' @author John Lennon L. Calorio
#' @export
has_negative <- function(x, warn = TRUE) {
  .check_warn(warn)
  v <- .coerce_to_atomic(x, "has_negative", warn)
  if (is.null(v)) return(FALSE)
  v_num <- suppressWarnings(as.numeric(v))
  any(v_num < 0, na.rm = TRUE)
}


# -----------------------------------------------------------------------------
# has_inf
# -----------------------------------------------------------------------------

#' @title Check if an Object Contains Infinite Values
#'
#' @description
#' Tests whether an object contains at least one \code{Inf} or \code{-Inf}
#' value. Accepts tabular objects, plain lists, or atomic vectors. Returns a
#' single logical value (\code{TRUE} or \code{FALSE}).
#'
#' @param x An R object to test.
#' @param warn Logical scalar. If \code{TRUE} (default), a warning is emitted
#'   when \code{x} is \code{NULL} or of an unsupported type.
#'
#' @details
#' Both positive (\code{Inf}) and negative (\code{-Inf}) infinite values
#' trigger a \code{TRUE} result. \code{NaN} and \code{NA} are not infinite
#' (consistent with \code{is.infinite()} in base R). For tabular inputs the
#' entire table is scanned; for vectors or list inputs all elements are scanned
#' after unlisting.
#'
#' @return \code{TRUE} if at least one \code{Inf} or \code{-Inf} is found,
#'   \code{FALSE} otherwise. Never \code{NA}.
#'
#' @examples
#' has_inf(c(1, Inf, 3))           # TRUE
#' has_inf(c(1, -Inf, 3))          # TRUE
#' has_inf(c(1, NaN, NA))          # FALSE - NaN/NA are not Inf
#' has_inf(c(1, 2, 3))             # FALSE
#'
#' df <- data.frame(x = c(1, Inf), y = c(3, 4))
#' has_inf(df)                      # TRUE
#'
#' mat <- matrix(c(1, 2, -Inf, 4), nrow = 2)
#' has_inf(mat)                     # TRUE
#'
#' @author John Lennon L. Calorio
#' @export
has_inf <- function(x, warn = TRUE) {
  .check_warn(warn)
  v <- .coerce_to_atomic(x, "has_inf", warn)
  if (is.null(v)) return(FALSE)
  v_num <- suppressWarnings(as.numeric(v))
  any(is.infinite(v_num))
}