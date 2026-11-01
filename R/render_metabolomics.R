#' Render a standalone metabolomics analysis report
#'
#' Runs the full direct-injection metabolomics workflow (preprocessing, batch
#' correction views, PCA, OPLS-DA / PLS-DA, fold change, t-test / ANOVA, volcano
#' plot, ROC) on one dataset and writes a single self-contained HTML report
#' (all tables and figures embedded) plus an \code{.rds} of the results to a
#' \code{results} folder next to the input file. The workflow lives in the
#' package, so projects never hold their own copy.
#'
#' @section Data layout:
#' Samples are rows and features are columns in \code{features}; \code{metadata}
#' has one row per sample. Supply either
#' \itemize{
#'   \item two files (\code{features}, \code{metadata}), or one workbook with two
#'     sheets (\code{features_sheet}, \code{metadata_sheet}), or
#'   \item one table holding metadata and features together (give only
#'     \code{features}); all non-numeric columns plus \code{sample_id_col},
#'     \code{group_col}, \code{batch_col}, \code{injection_col} and
#'     \code{metadata_cols} are treated as metadata, the rest as features.
#' }
#' Supported formats: \code{.csv}, \code{.tsv}/\code{.txt}, \code{.xlsx}/\code{.xls}.
#' A data frame can also be given instead of a path (then set \code{output_dir}).
#' Features may carry a sample ID column, or be in the same row order as
#' \code{metadata}; a features table with samples as columns is transposed
#' automatically.
#'
#' @section Overriding function arguments:
#' Settings of the underlying functions are passed through \code{...} as one
#' named list per function, so the same argument name in different functions is
#' never ambiguous:
#' \preformatted{
#' render_metabolomics("serum.xlsx", features_sheet = "features",
#'   metadata_sheet = "metadata",
#'   run_foldchange = list(arrange = c("Control", "Case"), up = 2),
#'   plot_volcano   = list(pval = 0.01),
#'   run_DIpreprocess = list(missing_threshold = 0.3))
#' }
#' Valid names: \code{run_DIpreprocess}, \code{plot_beforeafter},
#' \code{plot_dist_beforeafter}, \code{run_pca}, \code{run_pls},
#' \code{plot_score} (PCA and PLS plots), \code{run_foldchange}, \code{run_diff},
#' \code{get_volcanodata}, \code{plot_volcano}, \code{plot_meanmap},
#' \code{run_auc}. Argument names are checked against the function, and a
#' supplied value replaces the workflow's default (including \code{NULL}).
#' Data inputs (\code{x}) cannot be overridden.
#'
#' @param features Path to the features file (or the combined file), or a data frame.
#' @param metadata Path to the metadata file, or a data frame. Omit when
#'   \code{metadata_sheet} is in the \code{features} workbook, or for a combined table.
#' @param features_sheet,metadata_sheet Sheet names (workbooks only).
#' @param metadata_cols Extra metadata column names for a combined table.
#' @param name Report name; default the features file name.
#' @param output_dir Folder that will contain \code{results/}. Default: the
#'   folder of the \code{features} file.
#' @param sample_id_col,group_col,batch_col,injection_col Metadata column names.
#' @param qc_types Group labels identifying QC samples.
#' @param normalize_method Normalization method (e.g. \code{"median"}, \code{"sum"},
#'   \code{"pqn_global"}); see \code{dimsprepr::run_DIpreprocess}.
#' @param norm_factor_col Numeric metadata column with per-sample normalization
#'   values (e.g. specific gravity for urine). When given, it overrides
#'   \code{normalize_method}.
#' @param normalize_qc_method How QC samples are normalized when
#'   \code{norm_factor_col} is used (QC rows usually have no value): \code{"median"},
#'   \code{"mean"} or \code{"none"}. Default \code{"median"}; ignored otherwise.
#' @param control_group Reference group for ROC (two-group designs).
#' @param ... Per-function argument lists, see above.
#' @param overwrite If \code{FALSE} (default) an existing report is never replaced:
#'   a timestamp (\code{_YYYYmmdd-HHMMSS}) is appended to the new file names.
#'   \code{TRUE} replaces the existing report and results.
#' @param timestamp If \code{TRUE}, always append a timestamp to the file names.
#' @param quiet Set \code{FALSE} to see Quarto's progress output.
#'
#' @return Invisibly, the path of the HTML report.
#' @export
render_metabolomics <- function(features, metadata = NULL,
                                features_sheet = NULL, metadata_sheet = NULL,
                                metadata_cols = NULL, name = NULL, output_dir = NULL,
                                sample_id_col = "Sample", group_col = "Group",
                                batch_col = "Batch", injection_col = "InjectionSequence",
                                qc_types = c("QC", "SQC", "EQC"),
                                normalize_method = "median", norm_factor_col = NULL,
                                normalize_qc_method = "median", control_group = NULL,
                                ..., overwrite = FALSE, timestamp = FALSE, quiet = TRUE) {
  for (p in c("quarto", "dimsprepr", "gt"))
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required.", call. = FALSE)
  src <- system.file("workflows", package = "pondeR")
  if (!nzchar(src)) stop("Bundled workflow not found; reinstall pondeR.", call. = FALSE)

  over <- .wf_check_overrides(list(...))

  if (is.null(output_dir)) {
    if (!is.character(features))
      stop("Give `output_dir` when `features` is a data frame.", call. = FALSE)
    output_dir <- dirname(normalizePath(features))
  }
  if (is.null(name))
    name <- if (is.character(features)) tools::file_path_sans_ext(basename(features)) else "dataset"
  name <- gsub("[^A-Za-z0-9._-]+", "_", name)

  d <- .wf_prepare(features, metadata, features_sheet, metadata_sheet, metadata_cols,
                   sample_id_col, group_col, batch_col, injection_col)

  # Preprocess here so input problems surface immediately, not inside Quarto
  prep_args <- list(
    x = d$features, metadata = d$metadata, sample_id_col = sample_id_col,
    group_col = group_col, qc_types = qc_types, batch_col = batch_col,
    injection_col = injection_col, normalize_method = normalize_method,
    norm_factor_col = norm_factor_col,
    normalize_qc_method = if (is.null(norm_factor_col)) "none" else normalize_qc_method,
    verbose = FALSE)
  di <- do.call(dimsprepr::run_DIpreprocess, .wf_merge(prep_args, over$run_DIpreprocess))
  if (!is.null(di$error))
    stop("Preprocessing failed: ", di$error, call. = FALSE)
  if (anyNA(di$data_nonpls))
    warning("Preprocessed data still contain missing values (e.g. drift correction ",
            "failed for some features); PCA/PLS steps will fail. Check the QC samples per batch.",
            call. = FALSE)

  res_dir <- file.path(output_dir, "results")
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  res_dir <- normalizePath(res_dir)
  name <- .wf_output_name(res_dir, name, overwrite, timestamp)

  # Render from a temp copy so nothing but the outputs is written to disk
  tmp <- tempfile("pondeR_workflow_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file.copy(file.path(src, "metabolomics_report.qmd"), tmp)
  spec_file <- file.path(tmp, "spec.rds")
  saveRDS(list(di = di, over = over, name = name, group_col = group_col,
               sample_id_col = sample_id_col, batch_col = batch_col,
               injection_col = injection_col, qc_types = qc_types,
               control_group = control_group,
               results_file = file.path(res_dir, paste0(name, "_results.rds"))),
          spec_file)

  report <- paste0(name, "_report.html")
  tryCatch(
    quarto::quarto_render(file.path(tmp, "metabolomics_report.qmd"),
                          execute_params = list(spec_file = spec_file),
                          execute_dir = tmp, output_file = report, quiet = quiet),
    error = function(e) stop("Report rendering failed; rerun with `quiet = FALSE` for details.\n",
                             conditionMessage(e), call. = FALSE))

  out <- file.path(res_dir, report)
  file.copy(file.path(tmp, report), out, overwrite = TRUE)
  message("Report: ", out)
  invisible(out)
}

# ---- internals -------------------------------------------------------------

.wf_functions <- function() {
  s3 <- function(g, cls) formals(getS3method(g, cls, envir = asNamespace("pondeR")))
  list(
    run_DIpreprocess      = formals(dimsprepr::run_DIpreprocess),
    plot_beforeafter      = formals(dimsprepr::plot_beforeafter),
    plot_dist_beforeafter = formals(dimsprepr::plot_dist_beforeafter),
    run_pca               = formals(run_pca),
    run_pls               = formals(run_pls),
    plot_score            = c(s3("plot_score", "run_pca"), s3("plot_score", "run_pls")),
    run_foldchange        = formals(run_foldchange),
    run_diff              = formals(run_diff),
    get_volcanodata       = formals(get_volcanodata),
    plot_volcano          = formals(plot_volcano),
    plot_meanmap          = formals(plot_meanmap),
    run_auc               = formals(run_auc))
}

# Validate `...`: named lists keyed by function, with argument names that exist
.wf_check_overrides <- function(over) {
  if (!length(over)) return(list())
  fns <- .wf_functions()
  if (is.null(names(over)) || any(names(over) == ""))
    stop("Extra arguments must be named lists, e.g. run_foldchange = list(up = 2). ",
         "Valid names: ", paste(names(fns), collapse = ", "), call. = FALSE)
  for (fn in names(over)) {
    if (!fn %in% names(fns))
      stop("Unknown function '", fn, "'. Valid names: ", paste(names(fns), collapse = ", "),
           call. = FALSE)
    o <- over[[fn]]
    if (!is.list(o) || (length(o) && (is.null(names(o)) || any(names(o) == ""))))
      stop("`", fn, "` must be a named list of arguments, e.g. ", fn, " = list(arg = value).",
           call. = FALSE)
    ok <- names(fns[[fn]])
    # Functions taking `...` (forwarded elsewhere) cannot be checked by name
    bad <- if ("..." %in% ok && !fn %in% c("get_volcanodata")) character(0) else setdiff(names(o), ok)
    bad <- union(bad, intersect(names(o), c("x", "res")))
    if (length(bad))
      stop("Argument(s) not allowed for ", fn, ": ", paste(bad, collapse = ", "),
           ". Allowed: ", paste(setdiff(ok, c("x", "res", "...")), collapse = ", "),
           call. = FALSE)
  }
  over
}

# Replace defaults with user values (NULL values are kept)
.wf_merge <- function(defaults, user) {
  if (length(user)) defaults[names(user)] <- user
  defaults
}

# Called from the report: run `fn` with workflow defaults, overridden by the user's
.wf_call <- function(fn, defaults, over) {
  f <- if (fn %in% c("run_DIpreprocess", "plot_beforeafter", "plot_dist_beforeafter"))
    getExportedValue("dimsprepr", fn) else getExportedValue("pondeR", fn)
  do.call(f, .wf_merge(defaults, over[[fn]]))
}

.wf_read <- function(x, sheet, what) {
  if (is.data.frame(x) || is.matrix(x)) return(as.data.frame(x, check.names = FALSE))
  if (!is.character(x) || length(x) != 1 || !file.exists(x))
    stop(what, " file not found: ", x, call. = FALSE)
  ext <- tolower(tools::file_ext(x))
  switch(ext,
    csv = utils::read.csv(x, check.names = FALSE),
    tsv = , txt = utils::read.delim(x, check.names = FALSE),
    xlsx = , xls = {
      if (!requireNamespace("openxlsx", quietly = TRUE))
        stop("Package 'openxlsx' is required for Excel files.", call. = FALSE)
      sheets <- openxlsx::getSheetNames(x)
      if (is.null(sheet)) {
        if (length(sheets) == 1) sheet <- sheets
        else stop("Specify the ", what, " sheet. Sheets in ", basename(x), ": ",
                  paste(sheets, collapse = ", "), call. = FALSE)
      }
      if (!sheet %in% sheets)
        stop("Sheet '", sheet, "' not found in ", basename(x), ". Sheets: ",
             paste(sheets, collapse = ", "), call. = FALSE)
      openxlsx::read.xlsx(x, sheet = sheet, check.names = FALSE)
    },
    stop("Unsupported ", what, " file type '.", ext, "' (use csv, tsv/txt or xlsx).", call. = FALSE))
}

# Load + validate; returns list(features, metadata) with matching rows
.wf_prepare <- function(features, metadata, features_sheet, metadata_sheet,
                        metadata_cols, sample_id_col, group_col, batch_col, injection_col) {
  # Workbook with metadata in a second sheet of the same file
  if (is.null(metadata) && !is.null(metadata_sheet)) metadata <- features
  ft <- .wf_read(features, features_sheet, "features")

  if (is.null(metadata)) {
    # Combined table
    meta_names <- intersect(names(ft), c(sample_id_col, group_col, batch_col, injection_col, metadata_cols))
    non_num <- names(ft)[!vapply(ft, is.numeric, logical(1))]
    meta_names <- union(meta_names, non_num)
    message("Combined table: using ", length(meta_names), " metadata column(s): ",
            paste(meta_names, collapse = ", "), " (add others via `metadata_cols`).")
    md <- ft[meta_names]
    ft <- ft[setdiff(names(ft), meta_names)]
  } else {
    md <- .wf_read(metadata, metadata_sheet, "metadata")
  }

  if (!sample_id_col %in% names(md))
    stop("Metadata has no sample ID column '", sample_id_col, "'. Columns: ",
         paste(names(md), collapse = ", "), call. = FALSE)
  ids <- as.character(md[[sample_id_col]])
  if (anyDuplicated(ids)) stop("Sample IDs in metadata are not unique.", call. = FALSE)
  for (cl in c(group_col, batch_col, injection_col))
    if (!cl %in% names(md)) stop("Metadata has no column '", cl, "'. Columns: ",
                                 paste(names(md), collapse = ", "), call. = FALSE)

  # Orient features: samples in rows, ordered like metadata
  if (sample_id_col %in% names(ft)) {
    rownames(ft) <- as.character(ft[[sample_id_col]]); ft[[sample_id_col]] <- NULL
  } else if (all(ids %in% names(ft))) {
    first <- ft[[1]]
    if (!is.numeric(first)) { r <- ft[-1]; rownames(r) <- as.character(first); ft <- r }
    ft <- as.data.frame(t(ft[ids]), check.names = FALSE)
    message("Features had samples as columns; transposed.")
  } else if (!is.numeric(ft[[1]]) && all(ids %in% as.character(ft[[1]]))) {
    rownames(ft) <- as.character(ft[[1]]); ft <- ft[-1]
  } else if (nrow(ft) == nrow(md)) {
    message("No sample IDs in features; assuming the same row order as metadata.")
    rownames(ft) <- ids
  } else {
    stop("Cannot match features to metadata: features have ", nrow(ft), " rows, metadata ",
         nrow(md), ". Include a '", sample_id_col, "' column in the features.", call. = FALSE)
  }
  if (!all(ids %in% rownames(ft)))
    stop(sum(!ids %in% rownames(ft)), " metadata sample(s) missing from features, e.g. ",
         paste(utils::head(ids[!ids %in% rownames(ft)], 3), collapse = ", "), call. = FALSE)
  ft <- ft[ids, , drop = FALSE]
  rownames(md) <- ids

  bad <- names(ft)[!vapply(ft, is.numeric, logical(1))]
  if (length(bad))
    stop(length(bad), " non-numeric feature column(s), e.g. ", paste(utils::head(bad, 3), collapse = ", "),
         ". For a combined table list extra metadata columns in `metadata_cols`.", call. = FALSE)
  list(features = ft, metadata = md)
}

# Report/results base name; stamped when asked, or when an earlier report exists
.wf_output_name <- function(dir, name, overwrite, timestamp) {
  taken <- file.exists(file.path(dir, c(paste0(name, "_report.html"), paste0(name, "_results.rds"))))
  if (isTRUE(timestamp) || (!isTRUE(overwrite) && any(taken)))
    name <- paste0(name, "_", format(Sys.time(), "%Y%m%d-%H%M%S"))
  name
}
