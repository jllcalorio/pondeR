# =============================================================================
# save_plots.R
# Part of the pondeR package
# =============================================================================

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.sanitize_plot_filename <- function(x) {
  x <- gsub('[/\\\\:*?"<>|]', "", x)
  x <- trimws(x)
  if (!nzchar(x)) stop("Plot name resolves to an empty string after removing illegal characters.")
  x
}

#' @keywords internal
.unique_plot_names <- function(nms) {
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
        nms[i]    <- paste0(nm, " (", idx, ")")
        seen[[nm]] <- idx + 1L
      }
    }
  }
  nms
}

#' @keywords internal
.resolve_plot_path <- function(folder, base_name, ext,
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

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# ── Plot-type detection ───────────────────────────────────────────────────────

#' @keywords internal
.is_ggplot <- function(obj) inherits(obj, "ggplot")

#' @keywords internal
.is_plotly <- function(obj) inherits(obj, "plotly")

#' @keywords internal
.is_grob <- function(obj) inherits(obj, c("grob", "gTree", "gtable"))

#' @keywords internal
.is_recordedplot <- function(obj) inherits(obj, "recordedplot")

#' @keywords internal
.is_base_plot_fn <- function(obj) is.function(obj)

# pondeR plot-producing classes (expand as new functions are added)
.PONDER_PLOT_CLASSES <- c(
  "plot_diff", "run_auc", "plot_score", "plot_scree", "plot_volcano"
)

#' @keywords internal
.is_ponder_plot_result <- function(obj) {
  # Matches any class starting with "plot_" or in the known list
  any(grepl("^plot_", class(obj))) ||
    any(class(obj) %in% .PONDER_PLOT_CLASSES)
}

#' @keywords internal
.extract_plots_from_ponder <- function(obj) {
  # Recursively extract ggplot / grob objects from a pondeR result
  candidates <- Filter(function(el)
    .is_ggplot(el) || .is_grob(el) || .is_plotly(el) || .is_recordedplot(el),
    obj)
  if (!length(candidates) && is.list(obj)) {
    # One level deeper (e.g. obj$plots$...)
    inner <- unlist(obj, recursive = FALSE)
    candidates <- Filter(function(el)
      .is_ggplot(el) || .is_grob(el) || .is_plotly(el) || .is_recordedplot(el),
      inner)
  }
  candidates
}

#' @keywords internal
.is_any_plot <- function(obj) {
  .is_ggplot(obj)          ||
    .is_plotly(obj)        ||
    .is_grob(obj)          ||
    .is_recordedplot(obj)  ||
    .is_base_plot_fn(obj)  ||
    .is_ponder_plot_result(obj)
}

# ── Device helpers ────────────────────────────────────────────────────────────

#' @keywords internal
.open_device <- function(path, ext, width, height, dpi, bg) {
  switch(ext,
    png  = grDevices::png( path, width = width, height = height,
                           units = "in", res = dpi, bg = bg),
    jpg  = ,
    jpeg = grDevices::jpeg(path, width = width, height = height,
                           units = "in", res = dpi, bg = bg, quality = 95),
    tiff = ,
    tif  = grDevices::tiff(path, width = width, height = height,
                           units = "in", res = dpi, bg = bg,
                           compression = "lzw"),
    bmp  = grDevices::bmp( path, width = width, height = height,
                           units = "in", res = dpi, bg = bg),
    pdf  = grDevices::pdf( path, width = width, height = height, bg = bg),
    svg  = grDevices::svg( path, width = width, height = height, bg = bg),
    eps  = grDevices::postscript(path, width = width, height = height,
                                 horizontal = FALSE, onefile = FALSE,
                                 paper = "special"),
    stop(sprintf("Cannot open a graphics device for extension '.%s'.", ext),
         call. = FALSE)
  )
}

#' @keywords internal
.write_one_plot <- function(obj, path, ext, width, height, dpi, transparent) {
  bg <- if (transparent) "transparent" else "white"

  # ── pondeR result: extract first renderable plot ──────────────────────────
  if (.is_ponder_plot_result(obj)) {
    plots <- .extract_plots_from_ponder(obj)
    if (!length(plots))
      stop(sprintf(
        "The pondeR result object (class '%s') contains no extractable plot. ",
        paste(class(obj), collapse = "/")
      ), call. = FALSE)
    if (length(plots) > 1L)
      message(sprintf(
        "  Note: '%s' contains %d plots; exporting the first one ('%s'). ",
        paste(class(obj), collapse = "/"),
        length(plots), names(plots)[1L] %||% "unnamed"
      ))
    obj <- plots[[1L]]
  }

  # ── ggplot ────────────────────────────────────────────────────────────────
  if (.is_ggplot(obj)) {
    if (ext %in% c("pdf", "svg")) {
      ggplot2::ggsave(filename = path, plot = obj,
                      width = width, height = height,
                      bg = bg, units = "in")
    } else if (ext == "eps") {
      ggplot2::ggsave(filename = path, plot = obj,
                      width = width, height = height,
                      device = grDevices::postscript, units = "in")
    } else {
      ggplot2::ggsave(filename = path, plot = obj,
                      width = width, height = height,
                      dpi = dpi, bg = bg, units = "in")
    }
    return(invisible(NULL))
  }

  # ── plotly ────────────────────────────────────────────────────────────────
  if (.is_plotly(obj)) {
    if (ext == "html") {
      if (!requireNamespace("htmlwidgets", quietly = TRUE))
        stop("Package 'htmlwidgets' needed for HTML export. ",
             "Install with: install.packages(\"htmlwidgets\")", call. = FALSE)
      htmlwidgets::saveWidget(obj, file = path, selfcontained = TRUE)
    } else {
      if (!requireNamespace("plotly", quietly = TRUE))
        stop("Package 'plotly' needed for static plotly export. ",
             "Install with: install.packages(\"plotly\")", call. = FALSE)
      plotly::save_image(obj, file = path,
                         width  = width  * dpi,
                         height = height * dpi,
                         scale  = 1)
    }
    return(invisible(NULL))
  }

  # ── grob / gtable ─────────────────────────────────────────────────────────
  if (.is_grob(obj)) {
    if (!requireNamespace("grid", quietly = TRUE))
      stop("Package 'grid' needed for grob export.", call. = FALSE)
    .open_device(path, ext, width, height, dpi, bg)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.draw(obj)
    return(invisible(NULL))
  }

  # ── recordedplot ──────────────────────────────────────────────────────────
  if (.is_recordedplot(obj)) {
    .open_device(path, ext, width, height, dpi, bg)
    on.exit(grDevices::dev.off(), add = TRUE)
    grDevices::replayPlot(obj)
    return(invisible(NULL))
  }

  # ── zero-argument base-R drawing function ─────────────────────────────────
  if (.is_base_plot_fn(obj)) {
    .open_device(path, ext, width, height, dpi, bg)
    on.exit(grDevices::dev.off(), add = TRUE)
    obj()
    return(invisible(NULL))
  }

  stop(sprintf(
    paste0(
      "Unsupported plot type '%s'. Accepted types: ggplot, plotly, grob/gtable, ",
      "recordedplot, a zero-argument base-R drawing function, or a pondeR result ",
      "object (classes: plot_diff, run_auc, plot_score, plot_scree, plot_volcano, ",
      "or any class beginning with 'plot_')."
    ),
    paste(class(obj), collapse = "/")
  ), call. = FALSE)
}

# ── Main function ──────────────────────────────────────────────────────────────

#' @title Save Plots to Image or Vector Files
#'
#' @description
#' Exports a single plot or a named list of plots to individual files.
#' Supported plot types include \pkg{ggplot2} objects, \pkg{plotly} widgets,
#' \pkg{grid} grobs/gtables, \code{recordedplot} objects, zero-argument
#' base-R drawing functions, and pondeR result objects produced by functions
#' such as \code{plot_diff()}, \code{run_auc()}, \code{plot_score()},
#' \code{plot_scree()}, and \code{plot_volcano()}.
#'
#' When \code{x} is a single plot, \code{filename} is used as the output file
#' name.  When \code{x} is a list of plots, \code{filename} is ignored and
#' each element's list name is used instead.  If \code{x} contains no direct
#' plot objects at the top level but holds a list of plots one level deep
#' (e.g. \code{result$plots}), that nested list is transparently unpacked.
#'
#' @param x A single plot object or a named list of plot objects.  Accepted
#'   types: \code{ggplot}, \code{plotly}, \code{grob}/\code{gtable},
#'   \code{recordedplot}, a zero-argument base-R drawing function, or a pondeR
#'   result object.  One level of list nesting is automatically unpacked when
#'   no direct plot objects are found at the top level.
#' @param filename A single character string for the output file name.  Used
#'   when \code{x} is a single plot.  When \code{x} is a multi-plot list this
#'   argument is ignored and each element's name is used as the file name.
#'   The extension in \code{filename} sets the format unless \code{format} is
#'   also supplied.  Defaults to \code{"plot"} when omitted.
#' @param format Optional single character string overriding the file
#'   extension for all exported files, e.g. \code{"pdf"}, \code{"svg"},
#'   \code{"png"}.  Case-insensitive; accepts with or without a leading dot.
#'   When supplied, the extension in \code{filename} is ignored.
#' @param width Numeric.  Output width in inches.  Default \code{15}.
#' @param height Numeric.  Output height in inches.  Default \code{12}.
#' @param dpi Numeric.  Resolution in dots per inch for raster formats.
#'   Default \code{300}.  Ignored for vector formats.
#' @param folder A single character string for the target directory.
#'   Supports nested folders (e.g., \code{"Main folder/sub folder"}).
#'   Created recursively if absent. Defaults to \code{NULL} (current working
#'   directory).
#' @param transparent Logical.  Use a transparent background.  Defaults to
#'   \code{FALSE}.  A warning is issued and white is used as fallback for
#'   JPEG, which does not support transparency.
#' @param overwrite Logical.  Overwrite existing files when \code{TRUE}.
#'   When \code{FALSE} (default) a Windows-style counter suffix is appended,
#'   e.g. \code{"plot (1).png"}.
#' @param fd_date_stamp Logical.  Appends \code{"MmmDDYYYY"} date to
#'   \code{folder}.  Default \code{FALSE}.
#' @param fd_time_stamp Logical.  Appends \code{"HHMM"} time to
#'   \code{folder}.  Default \code{FALSE}.
#' @param fn_date_stamp Logical.  Appends date to each file name.  Default
#'   \code{FALSE}.
#' @param fn_time_stamp Logical.  Appends time to each file name.  Default
#'   \code{FALSE}.
#'
#' @details
#' \strong{Single vs. multi-plot behaviour.}  A single plot object uses
#' \code{filename} exactly (after sanitisation).  A list with multiple
#' elements ignores \code{filename}; each element's list name becomes its
#' file name.  Illegal characters are stripped, and duplicate names are
#' resolved with \code{(1)}, \code{(2)}, etc.
#'
#' \strong{Unnamed list elements.}  Elements without a name receive the
#' placeholder name \code{"plot_N"} where \code{N} is their position.
#'
#' \strong{pondeR result objects.}  Objects whose class starts with
#' \code{"plot_"} or matches a known pondeR class are unpacked automatically.
#' The first extractable plot is exported; a note is printed when multiple
#' plots are present.
#'
#' \strong{plotly HTML export.}  Requires \pkg{htmlwidgets}.  Static raster
#' export of plotly objects requires \pkg{plotly} with a working
#' \pkg{reticulate}/kaleido installation:
#' \code{reticulate::py_install("kaleido")}.
#'
#' @return
#' Called for its side effect.  Returns the full normalised file path(s) of
#' the saved plot(s) invisibly as a character vector.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Single ggplot — filename is respected
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' save_plots(p, filename = "AUROC All (Neg)")
#'
#' # List of ggplots — filename is ignored; list names are used
#' p2 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'         geom_point()
#' save_plots(list(scatter = p, iris_plot = p2))
#'
#' # Accessing a list slot with multiple plots (auto-unpacked)
#' res <- run_auc(data = mydata, response = "group", predictors = c("A","B"))
#' save_plots(res$plots, folder = "figures")
#'
#' # Override format for all plots
#' save_plots(list(scatter = p, iris_plot = p2), format = "pdf")
#'
#' # Base-R plot via a zero-argument function
#' save_plots(list(hist_mpg = function() hist(mtcars$mpg, main = "MPG")))
#' }
#'
#' @author John Lennon L. Calorio
#' @export
save_plots <- function(
    x,
    filename      = "plot",
    format        = NULL,
    width         = 15,
    height        = 12,
    dpi           = 300,
    folder        = NULL,
    transparent   = FALSE,
    overwrite     = FALSE,
    fd_date_stamp = FALSE,
    fd_time_stamp = FALSE,
    fn_date_stamp = FALSE,
    fn_time_stamp = FALSE
) {

  # ── 0. Dependency ─────────────────────────────────────────────────────────
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("'ggplot2' is required. Install with: install.packages(\"ggplot2\")",
         call. = FALSE)

  # ── 1. Scalar logical validation ──────────────────────────────────────────
  .chk <- function(v, nm) {
    if (!is.logical(v) || length(v) != 1L || is.na(v))
      stop(sprintf("'%s' must be a single non-NA logical (TRUE or FALSE).", nm),
           call. = FALSE)
  }
  .chk(transparent,   "transparent")
  .chk(overwrite,     "overwrite")
  .chk(fd_date_stamp, "fd_date_stamp")
  .chk(fd_time_stamp, "fd_time_stamp")
  .chk(fn_date_stamp, "fn_date_stamp")
  .chk(fn_time_stamp, "fn_time_stamp")

  # ── 2. Numeric validation ─────────────────────────────────────────────────
  .chk_pos <- function(v, nm) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0)
      stop(sprintf("'%s' must be a single positive numeric value.", nm), call. = FALSE)
  }
  .chk_pos(width,  "width")
  .chk_pos(height, "height")
  .chk_pos(dpi,    "dpi")

  # ── 3. Format table ───────────────────────────────────────────────────────
  raster_fmts <- c("png", "jpg", "jpeg", "tiff", "tif", "bmp")
  vector_fmts <- c("pdf", "svg", "eps")
  all_fmts    <- c(raster_fmts, vector_fmts, "html")

  # ── 4. filename validation ─────────────────────────────────────────────────
  if (!is.character(filename) || length(filename) != 1L ||
      is.na(filename) || !nzchar(trimws(filename)))
    stop("'filename' must be a single non-empty character string.", call. = FALSE)
  filename <- trimws(filename)

  ext_pat  <- "\\.([a-zA-Z0-9]+)$"
  ext_hit  <- regmatches(filename, regexpr(ext_pat, filename))
  if (length(ext_hit) == 1L) {
    fn_ext    <- tolower(sub("^\\.", "", ext_hit))
    base_name <- sub(ext_pat, "", filename)
  } else {
    fn_ext    <- "png"
    base_name <- filename
  }
  base_name <- .sanitize_plot_filename(base_name)

  # ── 5. Resolve final format ───────────────────────────────────────────────
  if (!is.null(format)) {
    if (!is.character(format) || length(format) != 1L || is.na(format))
      stop("'format' must be a single character string, e.g. \"pdf\".", call. = FALSE)
    out_ext <- tolower(trimws(sub("^\\.", "", format)))
  } else {
    out_ext <- fn_ext
  }
  if (!out_ext %in% all_fmts)
    stop(sprintf("Unsupported format '.%s'. Accepted: %s.",
                 out_ext, paste(all_fmts, collapse = ", ")), call. = FALSE)

  if (transparent && out_ext %in% c("jpg", "jpeg")) {
    warning("JPEG does not support transparency; background will be white.",
            call. = FALSE)
    transparent <- FALSE
  }

  # ── 6. folder validation ───────────────────────────────────────────────────
  if (!is.null(folder)) {
    if (!is.character(folder) || length(folder) != 1L ||
        is.na(folder) || !nzchar(trimws(folder)))
      stop("'folder' must be a single non-empty character string or NULL.", call. = FALSE)
    # Remove illegal chars but preserve path separators (/ \) and drive colon (:)
    folder <- trimws(gsub('[*?"<>|]', "", folder))
    if (!nzchar(folder)) stop("'folder' resolves to an empty string after removing illegal characters.", call. = FALSE)
  }

  # ── 7. Capture unevaluated name of x ─────────────────────────────────────
  x_name_raw  <- tryCatch(deparse(substitute(x)), error = function(e) "plot")
  x_name_raw  <- paste(x_name_raw, collapse = "")
  # Replace any non-alphanumeric runs with a single underscore; strip edges
  x_name_safe <- gsub("[^[:alnum:].]+", "_", x_name_raw)
  x_name_safe <- gsub("^_+|_+$", "", x_name_safe)
  if (!nzchar(x_name_safe)) x_name_safe <- "plot"

  # ── 8. Coerce x into a named list of plot objects ─────────────────────────

  # Case A: single recognised plot object
  if (.is_any_plot(x)) {
    plots        <- stats::setNames(list(x), base_name)
    use_names    <- base_name          # honour filename
    multi_mode   <- FALSE

  } else if (is.list(x) && !inherits(x, "ggplot")) {

    if (!length(x))
      stop("'x' is an empty list. Provide at least one plot object.", call. = FALSE)

    direct_plot <- vapply(x, .is_any_plot, logical(1L))

    if (!any(direct_plot)) {
      # One level deep: flatten ALL nested lists that contain plots (not just the first)
      nested_lists <- Filter(is.list, x)
      if (!length(nested_lists))
        stop(paste0(
          "'x' contains no recognised plot objects (checked one level deep). ",
          "Accepted types: ggplot, plotly, grob/gtable, recordedplot, ",
          "a zero-argument drawing function, or a pondeR plot result."
        ), call. = FALSE)
      # Collect plots from every qualifying nested list, preserving names
      inner_plots <- unlist(
        lapply(nested_lists, function(lst) Filter(.is_any_plot, lst)),
        recursive = FALSE
      )
      if (!length(inner_plots))
        stop(paste0(
          "No plot objects found one level inside 'x'. ",
          "Check that you are passing the correct sub-list."
        ), call. = FALSE)
      x            <- inner_plots
      direct_plot  <- rep(TRUE, length(x))
    }

    # Warn about and drop non-plot elements
    bad <- !direct_plot
    if (any(bad)) {
      warning(sprintf(
        "Skipping %d element(s) that are not recognised plot objects: position(s) %s.",
        sum(bad), paste(which(bad), collapse = ", ")
      ), call. = FALSE)
      x           <- x[!bad]
      direct_plot <- direct_plot[!bad]
    }

    # Fill missing names
    nms <- if (is.null(names(x))) character(length(x)) else names(x)
    blank      <- !nzchar(trimws(nms))
    nms[blank] <- paste0("plot_", which(blank))
    names(x)   <- nms

    plots      <- x
    multi_mode <- length(plots) > 1L
    # For multi-plot lists: ignore filename; use per-element names
    # For single-element lists: use filename
    use_names  <- if (multi_mode) names(plots) else base_name

  } else {
    stop(paste0(
      "'x' must be a recognised plot object or a list of plot objects. ",
      "Accepted types: ggplot, plotly, grob/gtable, recordedplot, ",
      "a zero-argument drawing function, or a pondeR plot result."
    ), call. = FALSE)
  }

  # ── 9. Sanitise and deduplicate plot names ────────────────────────────────
  if (multi_mode) {
    clean_names <- vapply(names(plots), .sanitize_plot_filename, character(1L))
    clean_names <- .unique_plot_names(clean_names)
  } else {
    clean_names <- .sanitize_plot_filename(use_names)
  }

  # ── 10. Export ────────────────────────────────────────────────────────────
  saved_paths <- character(length(plots))

  for (i in seq_along(plots)) {
    nm       <- clean_names[i]
    out_path <- .resolve_plot_path(
      folder        = folder,
      base_name     = nm,
      ext           = out_ext,
      fd_date_stamp = fd_date_stamp,
      fd_time_stamp = fd_time_stamp,
      fn_date_stamp = fn_date_stamp,
      fn_time_stamp = fn_time_stamp,
      overwrite     = overwrite
    )

    tryCatch(
      .write_one_plot(
        obj         = plots[[i]],
        path        = out_path,
        ext         = out_ext,
        width       = width,
        height      = height,
        dpi         = dpi,
        transparent = transparent
      ),
      error = function(e) stop(sprintf(paste0(
        "Failed to save plot '%s' to '%s'.\n",
        "  - Ensure the directory is writable.\n",
        "  - For plotly static export, install kaleido: ",
        "reticulate::py_install(\"kaleido\").\n",
        "  Original error: %s"
      ), nm, out_path, conditionMessage(e)), call. = FALSE)
    )

    saved_paths[i] <- normalizePath(out_path, mustWork = FALSE)
    message(sprintf("Plot saved [%d/%d]: %s", i, length(plots), saved_paths[i]))
  }

  invisible(saved_paths)
}