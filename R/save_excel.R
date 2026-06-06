# =============================================================================
# save_excel.R
# Part of the pondeR package
# =============================================================================

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.sanitize_filename_xl <- function(x) {
  x <- gsub('[/\\\\:*?"<>|]', "", x)
  x <- gsub("\\.{2,}", ".", x)
  x <- trimws(x)
  if (!nzchar(x)) stop("Filename resolves to an empty string after removing illegal characters.")
  x
}

#' @keywords internal
.sanitize_sheet_name_xl <- function(x) {
  # Strip Excel-illegal characters and leading apostrophes
  # NOTE: colon must be escaped as a literal inside the character class
  x <- gsub("[/\\\\?*:\\[\\]\\x3A]", "", x)   # belt-and-suspenders colon removal
  x <- gsub(":", "", x, fixed = TRUE)           # fixed = TRUE catches any remaining colons
  x <- sub("^'+", "", x)
  x <- trimws(x)
  if (!nzchar(x)) x <- "Sheet"
  # Truncate to 25 chars (reserves room for " (999)" dedup suffix = 6 chars)
  .truncate_sheet <- function(s, limit = 25L) {
    if (nchar(s) <= limit) return(s)
    need_remove <- nchar(s) - limit
    # Remove vowels from the end first
    chars <- strsplit(s, "")[[1L]]
    vowel_idx <- rev(grep("[aeiouAEIOU]", chars))
    removed <- 0L
    for (idx in vowel_idx) {
      if (removed >= need_remove) break
      chars[idx] <- NA_character_
      removed <- removed + 1L
    }
    s2 <- paste(chars[!is.na(chars)], collapse = "")
    if (nchar(s2) > limit) {
      # Still too long — remove consonants from the end
      chars2 <- strsplit(s2, "")[[1L]]
      cons_idx <- rev(grep("[^aeiouAEIOU\\s]", chars2))
      still_need <- nchar(s2) - limit
      removed2 <- 0L
      for (idx in cons_idx) {
        if (removed2 >= still_need) break
        chars2[idx] <- NA_character_
        removed2 <- removed2 + 1L
      }
      s2 <- paste(chars2[!is.na(chars2)], collapse = "")
    }
    substr(s2, 1L, limit)
  }
  .truncate_sheet(x, 25L)
}

#' @keywords internal
.unique_names_xl <- function(nms, max_suffix = 999L) {
  counts <- table(nms)
  dups   <- names(counts)[counts > 1L]
  if (!length(dups)) return(nms)
  seen <- list()
  for (i in seq_along(nms)) {
    nm <- nms[i]
    if (nm %in% dups) {
      idx <- seen[[nm]] %||% 0L
      if (idx == 0L) {
        seen[[nm]] <- 1L
      } else {
        if (idx > max_suffix) warning(
          sprintf("More than %d sheets share the name '%s'; numbering may exceed suffix limit.", max_suffix, nm),
          call. = FALSE
        )
        nms[i]    <- paste0(nm, " (", idx, ")")
        seen[[nm]] <- idx + 1L
      }
    }
  }
  nms
}

#' @keywords internal
.resolve_path_xl <- function(folder, base_name, ext,
                              fd_date_stamp, fd_time_stamp,
                              fn_date_stamp, fn_time_stamp,
                              overwrite) {
  now <- Sys.time()
  if (!is.null(folder)) {
    if (fd_date_stamp || fd_time_stamp) {
      stamp <- c(
        if (fd_date_stamp) format(now, "%b%d%Y"),
        if (fd_time_stamp) format(now, "%H%M")
      )
      folder <- paste0(folder, "_", paste(stamp, collapse = "_"))
    }
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  }
  if (fn_date_stamp || fn_time_stamp) {
    stamp <- c(
      if (fn_date_stamp) format(now, "%b%d%Y"),
      if (fn_time_stamp) format(now, "%H%M")
    )
    base_name <- paste0(base_name, "_", paste(stamp, collapse = "_"))
  }
  make_path <- function(bn) {
    fname <- paste0(bn, ".", ext)
    if (!is.null(folder)) file.path(folder, fname) else fname
  }
  path <- make_path(base_name)
  if (!overwrite && file.exists(path)) {
    k <- 1L
    repeat {
      cand <- make_path(paste0(base_name, " (", k, ")"))
      if (!file.exists(cand)) { path <- cand; break }
      k <- k + 1L
    }
  }
  path
}

# Null-coalescing operator (defined once; skip if already in namespace)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# ── Tabular type checks ───────────────────────────────────────────────────────

#' @keywords internal
.is_tabular <- function(obj) {
  is.data.frame(obj) ||
    inherits(obj, c("tbl_df", "tbl", "data.table")) ||
    is.matrix(obj)
}

#' @keywords internal
.is_ponder_result <- function(obj) {
  # pondeR run_* function results carry a class starting with "run_" or any
  # registered pondeR class.  We extract the first data-frame-like element.
  any(grepl("^run_", class(obj))) ||
    any(class(obj) %in% c("run_auc", "run_assoc", "run_diff"))
}

#' @keywords internal
.extract_tables_from_ponder <- function(obj) {
  # Return a named list of data frames found inside a pondeR result object
  candidates <- Filter(.is_tabular, obj)
  if (length(candidates) == 0L) {
    # Recurse one level into list slots
    candidates <- Filter(.is_tabular, unlist(obj, recursive = FALSE))
  }
  candidates
}

# ── Main function ─────────────────────────────────────────────────────────────

#' @title Save Data Frames, Tibbles, or Matrices to an Excel File
#'
#' @description
#' Exports one or more data frames, tibbles, or matrices to a single
#' \code{.xlsx} file (or \code{.csv} for single-object exports).  When
#' \code{x} is a named list, each tabular element is written to its own
#' worksheet.  pondeR result objects (classes beginning with \code{run_}) are
#' automatically unpacked to extract their tabular components.
#' If \code{x} contains no direct tables but holds a list of tables one level
#' deep, that nested list is transparently unwrapped.
#'
#' @param x A data frame, tibble, matrix, a pondeR result object (e.g. output
#'   of \code{run_auc()}, \code{run_diff()}), or a named list whose elements
#'   are any combination of those types.  One level of nesting is automatically
#'   unpacked when no direct tabular elements are found at the top level.
#' @param filename A single character string for the output file name.  Used
#'   as-is when \code{x} is a single object.  Ignored when \code{x} is a list
#'   with multiple elements — each element's name is used as its sheet name and
#'   the list object's name becomes the file name.  Illegal characters are
#'   removed silently.  Defaults to \code{"output"} if not supplied.
#' @param folder A single character string giving the target directory.
#'   Supports nested folders (e.g., \code{"Main folder/sub folder"}).
#'   Created recursively if absent. Defaults to \code{NULL} (current working
#'   directory).
#' @param overwrite Logical.  When \code{TRUE} an existing file is silently
#'   overwritten.  When \code{FALSE} (default) a Windows-style counter suffix
#'   is appended, e.g. \code{"results (1).xlsx"}.
#' @param freeze_first_row Logical.  When \code{TRUE} (default) the first row
#'   of every worksheet is frozen.
#' @param freeze_first_col Logical.  When \code{TRUE} the first column of
#'   every worksheet is frozen.  Defaults to \code{FALSE}.
#' @param auto_col_width Logical.  When \code{TRUE} (default) column widths
#'   are auto-fitted.
#' @param fd_date_stamp Logical.  Appends \code{"MmmDDYYYY"} date to
#'   \code{folder}.  Default \code{FALSE}.
#' @param fd_time_stamp Logical.  Appends \code{"HHMM"} time to
#'   \code{folder}.  Default \code{FALSE}.
#' @param fn_date_stamp Logical.  Appends date to the file name.  Default
#'   \code{FALSE}.
#' @param fn_time_stamp Logical.  Appends time to the file name.  Default
#'   \code{FALSE}.
#'
#' @details
#' \strong{Sheet naming.}  Sheet names are taken from \code{names(x)}.
#' Illegal Excel characters are stripped, names are truncated to 25 characters
#' (vowels removed from the end first, then consonants if needed), the reserved
#' name \code{"History"} becomes \code{"History_"}, and duplicates get a
#' \code{" (N)"} suffix (up to \code{" (999)"}, hence the 25-character base
#' limit).
#'
#' \strong{Filename priority.}  The \code{filename} argument is prioritized
#' unless it is left as the default \code{"output"}.  If \code{filename} is
#' \code{"output"} and \code{x} is a list (or a pondeR result containing
#' multiple tables), the function automatically uses the variable name of
#' \code{x} as the filename for convenience.
#'
#' \strong{pondeR result objects.}  Objects whose class starts with
#' \code{run_} (e.g. \code{run_auc}, \code{run_diff}) are unpacked
#' automatically and their tabular slots are written as separate sheets.
#'
#' @return
#' Called for its side effect.  Returns the full normalised file path
#' invisibly as a character string.
#'
#' @examples
#' \dontrun{
#' # Single data frame — filename is respected
#' save_excel(mtcars, filename = "mtcars_export")
#'
#' # Named list — filename is ignored; sheet names come from list names
#' sheets <- list(Cars = mtcars, Iris = iris)
#' save_excel(sheets, filename = "ignored")
#'
#' # pondeR result object
#' res <- run_auc(data = mydata, response = "group", predictors = c("A","B"))
#' save_excel(res, filename = "auc_results")
#'
#' # With date + time stamp on the file name
#' save_excel(iris, filename = "iris", fn_date_stamp = TRUE, fn_time_stamp = TRUE)
#' }
#'
#' @author John Lennon L. Calorio
#' @export
save_excel <- function(
    x,
    filename         = "output",
    folder           = NULL,
    overwrite        = FALSE,
    freeze_first_row = TRUE,
    freeze_first_col = FALSE,
    auto_col_width   = TRUE,
    fd_date_stamp    = FALSE,
    fd_time_stamp    = FALSE,
    fn_date_stamp    = FALSE,
    fn_time_stamp    = FALSE
) {

  # ── 0. Dependency ─────────────────────────────────────────────────────────
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("'openxlsx' is required. Install with: install.packages(\"openxlsx\")",
         call. = FALSE)

  # ── 1. Scalar logical validation ──────────────────────────────────────────
  .chk <- function(v, nm) {
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("'%s' must be a single non-NA logical (TRUE or FALSE).", nm),
           call. = FALSE)
  }
  .chk(overwrite,        "overwrite")
  .chk(freeze_first_row, "freeze_first_row")
  .chk(freeze_first_col, "freeze_first_col")
  .chk(auto_col_width,   "auto_col_width")
  .chk(fd_date_stamp,    "fd_date_stamp")
  .chk(fd_time_stamp,    "fd_time_stamp")
  .chk(fn_date_stamp,    "fn_date_stamp")
  .chk(fn_time_stamp,    "fn_time_stamp")

  # ── 2. filename validation ─────────────────────────────────────────────────
  if (!is.character(filename) || length(filename) != 1L ||
      is.na(filename) || !nzchar(trimws(filename)))
    stop("'filename' must be a single non-empty character string, e.g. \"results\".",
         call. = FALSE)
  filename <- trimws(filename)

  ext_pat   <- "\\.([a-zA-Z0-9]+)$"
  ext_hit   <- regmatches(filename, regexpr(ext_pat, filename))
  if (length(ext_hit) == 1L) {
    supplied_ext <- tolower(sub("^\\.", "", ext_hit))
    base_name    <- sub(ext_pat, "", filename)
  } else {
    supplied_ext <- "xlsx"
    base_name    <- filename
  }
  if (!supplied_ext %in% c("xlsx", "csv"))
    stop(sprintf(
      "Unsupported extension '.%s'. Use '.xlsx' (recommended) or '.csv'.",
      supplied_ext), call. = FALSE)

  base_name <- .sanitize_filename_xl(base_name)

  # ── 3. folder validation ───────────────────────────────────────────────────
  if (!is.null(folder)) {
    if (!is.character(folder) || length(folder) != 1L ||
        is.na(folder) || !nzchar(trimws(folder)))
      stop("'folder' must be a single non-empty character string or NULL.",
           call. = FALSE)
    # Remove illegal chars but preserve path separators (/ \) and drive colon (:)
    folder <- trimws(gsub('[*?"<>|]', "", folder))
    if (!nzchar(folder)) stop("'folder' resolves to an empty string after removing illegal characters.", call. = FALSE)
  }

  # ── 4. Capture the unevaluated name of x for fallback naming ─────────────
  x_name_raw <- tryCatch(deparse(substitute(x)), error = function(e) "output")
  # Collapse multi-line deparsing; keep only safe characters
  x_name_raw <- paste(x_name_raw, collapse = "")
  x_name_safe <- gsub("[^[:alnum:]._]", "_", x_name_raw)
  x_name_safe <- gsub("_+", "_", x_name_safe)
  x_name_safe <- sub("_+$", "", x_name_safe)
  if (!nzchar(x_name_safe)) x_name_safe <- "output"

  # ── 5. Coerce x into a named list of data frames ──────────────────────────
  coerce_list <- function(lst) lapply(lst, function(el) {
    if (is.matrix(el)) as.data.frame(el) else el
  })

  if (.is_tabular(x)) {
    # --- Single tabular object: honour filename ---
    sheets <- stats::setNames(list(as.data.frame(x)), base_name)
    use_filename <- if (base_name == "output") x_name_safe else base_name

  } else if (.is_ponder_result(x)) {
    # --- pondeR result: unpack tables, use filename ---
    tbls <- .extract_tables_from_ponder(x)
    if (!length(tbls))
      stop("The pondeR result object contains no tabular data to export.",
           call. = FALSE)
    sheets <- coerce_list(tbls)
    use_filename <- if (base_name == "output") x_name_safe else base_name

  } else if (is.list(x)) {
    if (!length(x)) {
      warning(
        "'x' is an empty list. An empty Excel file will be exported.",
        call. = FALSE
      )
      sheets       <- stats::setNames(list(data.frame()), base_name)
      use_filename <- base_name
      n_sheets     <- 1L
      # Skip the rest of list processing and jump straight to export
      return(invisible(.export_xl_workbook(
        sheets, use_filename, supplied_ext, folder,
        fd_date_stamp, fd_time_stamp, fn_date_stamp, fn_time_stamp,
        overwrite, freeze_first_row, freeze_first_col, auto_col_width
      )))
    }

    # Check if any direct element is tabular / pondeR
    direct_tab  <- vapply(x, .is_tabular,       logical(1L))
    direct_pond <- vapply(x, .is_ponder_result,  logical(1L))
    any_direct  <- direct_tab | direct_pond

    if (!any(any_direct)) {
      # One level deep: look for lists that contain tables
      nested <- Filter(is.list, x)
      if (!length(nested))
        stop(paste0(
          "'x' contains no data frames, tibbles, matrices, or pondeR result objects ",
          "(checked one level deep). Check that you are passing the correct object."
        ), call. = FALSE)
      # Unpack first qualifying nested list
      x_inner <- nested[[1L]]
      inner_tab <- Filter(.is_tabular, x_inner)
      if (!length(inner_tab))
        stop("No tabular data found one level inside 'x'. Check your input structure.",
             call. = FALSE)
      x         <- inner_tab
      any_direct <- rep(TRUE, length(x))
    }

    # Keep only tabular / pondeR elements; warn about others
    bad <- !any_direct
    if (any(bad)) {
      warning(sprintf(
        "Skipping %d element(s) in 'x' that are neither tabular nor pondeR results: %s",
        sum(bad),
        paste(which(bad), collapse = ", ")
      ), call. = FALSE)
      x          <- x[!bad]
      any_direct <- any_direct[!bad]
    }

    # Expand pondeR results into their tables inline
    expanded <- list()
    for (i in seq_along(x)) {
      el <- x[[i]]
      nm <- names(x)[i] %||% paste0("Sheet", i)
      if (.is_ponder_result(el)) {
        tbls <- .extract_tables_from_ponder(el)
        if (!length(tbls)) {
          warning(sprintf("pondeR result at position %d has no tabular data; skipping.", i),
                  call. = FALSE)
          next
        }
        # Prefix table names with the element name
        pfx <- if (nzchar(nm)) paste0(nm, "_") else ""
        names(tbls) <- paste0(pfx, names(tbls) %||% seq_along(tbls))
        expanded <- c(expanded, tbls)
      } else {
        expanded[[nm]] <- as.data.frame(el)
      }
    }
    sheets <- coerce_list(expanded)

    # Prioritize provided filename over variable name fallback
    use_filename <- if (base_name == "output") x_name_safe else base_name

  } else {
    stop(paste0(
      "'x' must be a data frame, tibble, matrix, a pondeR result object, ",
      "or a list of those types."
    ), call. = FALSE)
  }

  n_sheets <- length(sheets)

  # ── 6. CSV path ───────────────────────────────────────────────────────────
  csv_mode <- (supplied_ext == "csv")
  if (csv_mode && n_sheets > 1L) {
    warning(sprintf(paste0(
      "%d sheets were prepared, but '.csv' supports only one. ",
      "Only the first sheet ('%s') will be saved. ",
      "Use '.xlsx' to export all sheets in one file."
    ), n_sheets, names(sheets)[1L]), call. = FALSE)
    sheets   <- sheets[1L]
    n_sheets <- 1L
  }

  # ── 7. Sheet-name sanitisation ────────────────────────────────────────────
  if (!csv_mode) {
    sht_names <- vapply(names(sheets), .sanitize_sheet_name_xl, character(1L))
    sht_names[sht_names == "History"] <- "History_"
    sht_names <- .unique_names_xl(sht_names)
    names(sheets) <- sht_names

    if (n_sheets > 255L)
      warning(sprintf(paste0(
        "The workbook contains %d sheets, exceeding Excel's recommended 255-sheet limit. ",
        "The file will be saved but may not open correctly in all Excel versions."
      ), n_sheets), call. = FALSE)
  }

  # ── 8. Resolve output path ────────────────────────────────────────────────
  out_path <- .resolve_path_xl(
    folder        = folder,
    base_name     = use_filename,
    ext           = supplied_ext,
    fd_date_stamp = fd_date_stamp,
    fd_time_stamp = fd_time_stamp,
    fn_date_stamp = fn_date_stamp,
    fn_time_stamp = fn_time_stamp,
    overwrite     = overwrite
  )

  # ── 9a. CSV export ────────────────────────────────────────────────────────
  if (csv_mode) {
    utils::write.csv(sheets[[1L]], file = out_path, row.names = FALSE)
    message(sprintf("File saved: %s", normalizePath(out_path, mustWork = FALSE)))
    return(invisible(normalizePath(out_path, mustWork = FALSE)))
  }

  # ── 9b. XLSX export ───────────────────────────────────────────────────────
  wb <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(
    textDecoration = "bold",
    border         = "Bottom",
    borderStyle    = "medium",
    fgFill         = "#D9E1F2",
    halign         = "CENTER",
    valign         = "CENTER",
    wrapText       = FALSE
  )

  for (sht in names(sheets)) {
    df <- sheets[[sht]]
    openxlsx::addWorksheet(wb, sheetName = sht)
    openxlsx::writeData(wb, sheet = sht, x = df,
                        startRow = 1L, startCol = 1L,
                        headerStyle = hdr_style, rowNames = FALSE)
    if (freeze_first_row || freeze_first_col)
      openxlsx::freezePane(wb, sheet = sht,
                           firstRow = freeze_first_row,
                           firstCol = freeze_first_col)
    if (auto_col_width && ncol(df) > 0L)
      openxlsx::setColWidths(wb, sheet = sht,
                             cols   = seq_len(ncol(df)),
                             widths = "auto")
  }

  tryCatch(
    openxlsx::saveWorkbook(wb, file = out_path, overwrite = TRUE),
    error = function(e) stop(sprintf(paste0(
      "Failed to save '%s'.\n",
      "  - Ensure the path is valid and the file is not open in another program.\n",
      "  Original error: %s"
    ), out_path, conditionMessage(e)), call. = FALSE)
  )

  message(sprintf("File saved: %s", normalizePath(out_path, mustWork = FALSE)))
  invisible(normalizePath(out_path, mustWork = FALSE))
}