#' Render the metabolomics analysis workflow for a project
#'
#' Renders the Quarto workflow bundled with pondeR (\code{inst/workflows}) using
#' the settings in a project-specific config file. The workflow lives only in the
#' package, so updating pondeR updates every project; a project folder needs just
#' its config file and data.
#'
#' @param config Path to an R script defining \code{datasets} (a named list, one
#'   entry per dataset) and optionally \code{common} and \code{out_dir}. See
#'   \code{system.file("workflows/config_template.R", package = "pondeR")}.
#'   \code{NULL} renders the built-in demo data.
#' @param output_dir Folder for the report and \code{results/}; also the working
#'   directory, so relative paths in the config resolve from here. Default: the
#'   folder containing \code{config}.
#' @param output_file Name of the HTML report.
#' @param quiet Set \code{FALSE} to see Quarto's progress and error output.
#'
#' @return Invisibly, the path of the rendered report.
#' @examples
#' \dontrun{
#' render_metabolomics("P1/P1_Serum/config.R")
#' }
#' @export
render_metabolomics <- function(config = NULL, output_dir = NULL,
                                output_file = "metabolomics_report.html", quiet = TRUE) {
  if (!requireNamespace("quarto", quietly = TRUE))
    stop("Package 'quarto' is required.", call. = FALSE)
  src <- system.file("workflows", package = "pondeR")
  if (!nzchar(src)) stop("Bundled workflow not found; reinstall pondeR.", call. = FALSE)
  if (!is.null(config)) {
    if (!file.exists(config)) stop("Config file not found: ", config, call. = FALSE)
    config <- normalizePath(config)
  }
  if (is.null(output_dir)) output_dir <- if (is.null(config)) getwd() else dirname(config)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_dir <- normalizePath(output_dir)

  # Render from a temp copy so nothing is written into the package or project
  tmp <- tempfile("pondeR_workflow_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file.copy(list.files(src, pattern = "[.]qmd$", full.names = TRUE), tmp)

  quarto::quarto_render(
    file.path(tmp, "metabolomics_workflow.qmd"),
    execute_params = list(config_file = config,
                          child_file = file.path(tmp, "_metabolomics_dataset.qmd")),
    execute_dir = output_dir, output_file = output_file, quiet = quiet)

  out <- file.path(output_dir, output_file)
  file.copy(file.path(tmp, output_file), out, overwrite = TRUE)
  message("Report: ", out)
  invisible(out)
}
