#' Performs Descriptive Statistics on a Dataframe
#'
#' @description
#' Generates a comprehensive descriptive statistics table from a given
#' dataframe, optionally stratified by a categorical variable. Provides
#' customisable summaries for both continuous and categorical variables,
#' including control over rounding, missing data reporting, labelling,
#' formatting, and table aesthetics. The output is a publication-ready
#' \code{gtsummary} / \code{gt} table.
#'
#' When \code{add_inferential_pvalues = TRUE}, p-values are computed
#' \emph{outside} \code{gtsummary} using \code{\link{run_diff}} (for
#' continuous variables) and \code{\link{run_assoc}} (for categorical
#' variables), then injected into the table. Each p-value cell carries a
#' letter superscript identifying the exact test used; a legend is appended
#' as table footnotes. Variables for which a p-value could not be computed
#' are shown as \code{---†} with a corresponding footnote.
#'
#' @param x Dataframe. The input dataframe.
#' @param summarize_what Character vector. Column names to include in the
#'   analysis. When \code{NULL} (default), all columns except those named in
#'   \code{split_by} and \code{strata_by} are included.
#' @param split_by String. A column name in \code{x} (must be categorical with
#'   at least 2 unique levels) by which to stratify the table. Required when
#'   \code{add_inferential_pvalues = TRUE}. Default is \code{NULL}.
#' @param filter Optional character vector of levels in \code{split_by} to 
#'   exclude from the analysis. If \code{x} is a \code{run_DIpreprocess} 
#'   object, defaults to the identified QC types.
#' @param split_by_header String. Display header for the \code{split_by}
#'   variable. Defaults to the value of \code{split_by}.
#' @param strata_by String. A column name in \code{x} (must be categorical
#'   with at least 2 unique levels) by which to stratify \emph{before}
#'   \code{split_by}. Cannot be combined with
#'   \code{add_inferential_pvalues = TRUE}. Default is \code{NULL}.
#' @param rename_variables List. Formulas of the form
#'   \code{list("original" ~ "new", ...)} used to relabel variables in the
#'   table.
#' @param continuous_statistics String. Summary statistic(s) for continuous
#'   variables. One of: \code{"meanSD"}, \code{"meanSD2"},
#'   \code{"medianIQR"}, \code{"mean"}, \code{"sd"}, \code{"median"},
#'   \code{"p0"}, \code{"p25"}, \code{"p50"}, \code{"p75"}, \code{"p100"},
#'   \code{"IQR"}. Default is \code{"meanSD"}.
#' @param categorical_statistics String. Summary statistic(s) for categorical
#'   variables. One of: \code{"n_percent"}, \code{"n"}, \code{"percent"}.
#'   Default is \code{"n_percent"}.
#' @param force_continuous Character vector. Column names to treat as
#'   continuous regardless of their storage type. The same coercion is applied
#'   before passing data to \code{\link{run_diff}}.
#' @param force_categorical Character vector. Column names to treat as
#'   categorical regardless of their storage type. The same coercion is
#'   applied before passing data to \code{\link{run_assoc}}.
#' @param n_digits_continuous Numeric vector \code{c(d1, d2)} or
#'   \code{"dynamic"}. Number of decimal places for the first and subsequent
#'   numeric tokens in continuous summary cells (e.g., mean and SD). Default
#'   is \code{c(2, 2)}.
#' @param n_digits_categorical Numeric vector \code{c(d1, d2)} or
#'   \code{"dynamic"}. Decimal places for count and percentage in categorical
#'   summary cells. Default is \code{c(0, 2)}.
#' @param zero_as_exp Logical. When \code{TRUE} (default), continuous summary
#'   tokens that round to \code{0} at the requested precision but are
#'   genuinely non-zero are displayed in scientific notation (e.g.,
#'   \code{1.00e-3}).
#' @param display_missing String. How to report missing values:
#'   \code{"ifany"} (default), \code{"no"}, or \code{"always"}.
#' @param missing_text String. Label for missing-value rows. Default is
#'   \code{"No data/missing"}.
#' @param missing_stat String. Format for missing-value counts: one of
#'   \code{"n_percent"} (default), \code{"n"}, \code{"percent"}.
#' @param include_missing_in_splits Logical. If \code{TRUE},
#'   missing values in \code{split_by} / \code{strata_by} columns are kept as
#'   an explicit level. Defaults to \code{FALSE}.
#' @param sort_categorical_variables_by String. Level ordering for categorical
#'   variables: \code{"alphanumeric"} (default) or \code{"frequency"}.
#' @param calc_percent_by String. Percentage base for categorical variables:
#'   \code{"column"} (default), \code{"row"}, \code{"cell"}, or \code{"total"}.
#'   \code{"total"} calculates percentages based on the total number of samples, study participants, or simply number of rows in `x`.
#' @param calc_col_percent_using String. Denominator when
#'   \code{calc_percent_by = "column"}: \code{"n_valid_in_column"} (default)
#'   or \code{"n_in_column"}.
#' @param add_inferential_pvalues Logical. If \code{TRUE}, p-values are
#'   computed via \code{\link{run_diff}} (continuous variables) and
#'   \code{\link{run_assoc}} (categorical variables) and injected into the
#'   table. Requires \code{split_by} to be non-\code{NULL} and
#'   \code{strata_by} to be \code{NULL}. Default is \code{FALSE}.
#' @param test_type_continuous String. Test strategy forwarded to
#'   \code{\link{run_diff}} when \code{add_inferential_pvalues = TRUE}. One
#'   of \code{"auto"} (default), \code{"parametric"}, or
#'   \code{"nonparametric"}. See \code{\link{run_diff}} for details.
#' @param test_type_categorical String. Test strategy forwarded to
#'   \code{\link{run_assoc}} when \code{add_inferential_pvalues = TRUE}. One
#'   of \code{"auto"} (default), \code{"chisq"}, or \code{"fisher"}. See
#'   \code{\link{run_assoc}} for details.
#' @param paired Logical. If \code{TRUE}, paired tests are used when computing
#'   p-values. Forwarded to both \code{\link{run_diff}} and
#'   \code{\link{run_assoc}}. Note that the descriptive statistics shown in
#'   the table body are always unpaired summaries. Default is \code{FALSE}.
#' @param n_digits_pvalues Numeric or \code{NULL}. Decimal places for
#'   displayed p-values. Default is \code{3}. Set to \code{NULL} for
#'   automatic formatting (up to 3 d.p. when \eqn{p < 0.10}, otherwise 1
#'   d.p.).
#' @param bold_significant_pvalues Logical. If \code{TRUE} (default), p-values
#'   below \code{bold_significant_pvalues_at} are bolded.
#' @param bold_significant_pvalues_at Numeric. Significance threshold for
#'   bolding. Default is \code{0.05}.
#' @param header String. Header for the variable-label column. Default is
#'   \code{"Variable"}.
#' @param show_n_header Logical. When \code{TRUE} (default), appends the total
#'   sample size to the \code{header} (e.g., "Variable (N = 100)").
#' @param bold_labels Logical. If \code{TRUE} (default), variable labels are
#'   bolded.
#' @param italicize_levels Logical. If \code{TRUE} (default), categorical
#'   level labels are italicised.
#' @param clean_table Logical. If \code{TRUE} (default), cells containing a
#'   zero count (e.g., \code{0 (0.00\%)}) are replaced with empty strings.
#' @param table_name String or \code{NULL}. Table title. \code{"auto"}
#'   (default) generates a title automatically. \code{NULL} suppresses the
#'   title.
#'
#' @return A \code{gt_tbl} object (via \code{gtsummary::as_gt()}) representing
#'   the descriptive (and optionally inferential) statistics table.
#'
#' @seealso
#' \code{\link{run_diff}} for the continuous-variable test engine.
#' \code{\link{run_assoc}} for the categorical-variable test engine.
#'
#' @importFrom openxlsx saveWorkbook writeData addWorksheet createWorkbook
#' @import gtsummary
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' library(gtsummary)
#' library(dplyr)
#'
#' set.seed(123)
#' sample_df <- data.frame(
#'   Age       = sample(c(18:65, NA), 100, replace = TRUE,
#'                      prob = c(rep(1, 48), 5)),
#'   Gender    = sample(c("Male", "Female", NA), 100, replace = TRUE,
#'                      prob = c(0.48, 0.48, 0.04)),
#'   Education = sample(c("Elementary", "High School",
#'                         "Bachelors", "Masters", NA), 100,
#'                      replace = TRUE,
#'                      prob = c(0.40, 0.30, 0.15, 0.10, 0.05)),
#'   Income    = sample(c(seq(5000, 40000, by = 5000), NA), 100,
#'                      replace = TRUE, prob = c(rep(1, 8), 2)),
#'   Smoker    = sample(c(1, 0, NA), 100, replace = TRUE,
#'                      prob = c(0.25, 0.65, 0.10)),
#'   Region    = sample(c("North", "South", "East", "West", NA), 100,
#'                      replace = TRUE,
#'                      prob = c(0.24, 0.24, 0.24, 0.24, .04))
#' )
#'
#' # Basic descriptive statistics (no grouping)
#' sample_df |> run_summarytable()
#'
#' # Descriptive statistics split by Gender
#' sample_df |>
#'   run_summarytable(
#'     split_by        = "Gender",
#'     split_by_header = "Participant Gender",
#'     rename_variables = list(Age    ~ "Age (Years)",
#'                             Income ~ "Annual Income",
#'                             Smoker ~ "Smoking Status")
#'   )
#'
#' # Median and IQR for continuous variables
#' sample_df |>
#'   run_summarytable(
#'     continuous_statistics  = "medianIQR",
#'     categorical_statistics = "n",
#'     n_digits_continuous    = c(0, 0),
#'     n_digits_categorical   = c(0, 0)
#'   )
#'
#' # Inferential p-values using run_diff (continuous) and run_assoc
#' # (categorical), auto test selection
#' sample_df |>
#'   run_summarytable(
#'     split_by                = "Gender",
#'     add_inferential_pvalues = TRUE
#'   )
#'
#' # Force non-parametric tests for continuous variables
#' sample_df |>
#'   run_summarytable(
#'     split_by                = "Gender",
#'     add_inferential_pvalues = TRUE,
#'     test_type_continuous    = "nonparametric",
#'     test_type_categorical   = "fisher"
#'   )
#'
#' # Paired inferential p-values (descriptive body remains unpaired)
#' sample_df |>
#'   run_summarytable(
#'     split_by                = "Gender",
#'     add_inferential_pvalues = TRUE,
#'     paired                  = TRUE
#'   )
#' }
#'
#' @export
run_summarytable <- function(
    x,
    summarize_what                = NULL,
    split_by                      = NULL,
    filter                        = NULL,
    split_by_header               = NULL,
    strata_by                     = NULL,
    rename_variables              = NULL,
    continuous_statistics         = "meanSD",
    categorical_statistics        = "n_percent",
    force_continuous              = NULL,
    force_categorical             = NULL,
    n_digits_continuous           = c(2, 2),
    n_digits_categorical          = c(0, 2),
    zero_as_exp                   = TRUE,
    display_missing               = "ifany",
    missing_text                  = "No data/missing",
    missing_stat                  = "n_percent",
    include_missing_in_splits     = FALSE,
    sort_categorical_variables_by = "alphanumeric",
    calc_percent_by               = "column",
    calc_col_percent_using        = "n_valid_in_column",
    add_inferential_pvalues       = FALSE,
    test_type_continuous          = c("auto", "parametric", "nonparametric"),
    test_type_categorical         = c("auto", "chisq", "fisher"),
    paired                        = FALSE,
    n_digits_pvalues              = 3,
    bold_significant_pvalues      = TRUE,
    bold_significant_pvalues_at   = 0.05,
    header                        = "Variable",
    show_n_header                 = TRUE,
    bold_labels                   = TRUE,
    italicize_levels              = TRUE,
    clean_table                   = TRUE,
    table_name                    = "auto"
) {

  # ---------------------------------------------------------------------------
  # 0. Coerce input to data frame
  # ---------------------------------------------------------------------------
  if (inherits(x, "run_DIpreprocess")) {
    if (is.null(filter)) filter <- eval(x$parameters$qc_types)
    if (is.null(split_by)) split_by <- x$parameters$group_col
    
    target_meta <- if (!is.null(x$metadata_merged)) x$metadata_merged else x$metadata
    target_data <- if (!is.null(x$data_nonpls_merged)) x$data_nonpls_merged else x$data_nonpls
    x <- cbind(target_meta, target_data)
  }

  x <- as.data.frame(x)

  if (!is.null(split_by) && !is.null(filter)) {
    if (!split_by %in% names(x)) stop(paste0("Split variable '", split_by, "' not found."))
    keep_rows <- !as.character(x[[split_by]]) %in% filter
    x <- x[keep_rows, , drop = FALSE]
  }

  # Filter by NAs in split/strata if they are not to be included in the analysis.
  # This ensures the total N (header and denominators) is consistent throughout.
  if (!include_missing_in_splits) {
    if (!is.null(split_by))  x <- x[!is.na(x[[split_by]]), , drop = FALSE]
    if (!is.null(strata_by)) x <- x[!is.na(x[[strata_by]]), , drop = FALSE]
  }

  # Add N to header if requested
  if (isTRUE(show_n_header)) {
    header <- sprintf("%s (N = %d)", header, dim(x)[1])
  }

  # Match engine test-type arguments early so bad values error immediately
  test_type_continuous  <- match.arg(test_type_continuous)
  test_type_categorical <- match.arg(test_type_categorical)

  # ---------------------------------------------------------------------------
  # 1. Piping guard
  # ---------------------------------------------------------------------------
  if (is.data.frame(summarize_what)) {
    stop(paste(
      "Argument 'summarize_what' cannot be a data frame.",
      "This error often happens when piping data and also passing",
      "the data frame as an argument.",
      "Incorrect: my_data |> run_summarytable(my_data, ...)",
      "Correct:   my_data |> run_summarytable(...)",
      sep = "\n"
    ))
  }

  # ---------------------------------------------------------------------------
  # 2. gtsummary version check
  # ---------------------------------------------------------------------------
  if (add_inferential_pvalues) {
    if (utils::packageVersion("gtsummary") < "1.7.0") {
      warning(
        "gtsummary >= 1.7.0 is recommended. ",
        "Currently installed: ", utils::packageVersion("gtsummary"), ". ",
        "Please update with: install.packages('gtsummary')"
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Apply force_categorical to the data frame first (affects split_by check)
  # ---------------------------------------------------------------------------
  if (!is.null(force_categorical)) {
    for (col_name in force_categorical) {
      x[[col_name]] <- as.factor(x[[col_name]])
    }
  }

  # ---------------------------------------------------------------------------
  # 4. split_by validation
  # ---------------------------------------------------------------------------
  if (!is.null(split_by)) {
    if (!is.character(split_by) || length(split_by) != 1)
      stop("'split_by' must be a single string representing a column name.")
    if (!(split_by %in% names(x)))
      stop(paste0("'split_by' column '", split_by, "' not found in data."))

    col_data <- x[[split_by]]
    if (!is.factor(col_data) && !is.character(col_data))
      stop(paste0("Column '", split_by,
                  "' must be categorical (factor or character) for splitting."))
    if (length(unique(stats::na.omit(col_data))) < 2)
      stop(paste0("Column '", split_by,
                  "' must have at least 2 unique non-missing levels."))
    if (is.null(split_by_header))
      split_by_header <- split_by
  }

  # ---------------------------------------------------------------------------
  # 5. strata_by validation
  # ---------------------------------------------------------------------------
  if (!is.null(strata_by)) {
    if (!is.character(strata_by) || length(strata_by) != 1)
      stop("'strata_by' must be a single string representing a column name.")
    if (!(strata_by %in% names(x)))
      stop(paste0("'strata_by' column '", strata_by, "' not found in data."))

    col_data <- x[[strata_by]]
    if (!is.factor(col_data) && !is.character(col_data))
      stop(paste0("Column '", strata_by,
                  "' must be categorical (factor or character)."))
    if (length(unique(stats::na.omit(col_data))) < 2)
      stop(paste0("Column '", strata_by,
                  "' must have at least 2 unique non-missing levels."))
  }

  # ---------------------------------------------------------------------------
  # 6. Inferential p-value pre-conditions
  # ---------------------------------------------------------------------------
  if (add_inferential_pvalues) {
    if (is.null(split_by))
      stop("'add_inferential_pvalues = TRUE' requires 'split_by' to be non-NULL.")
    if (!is.null(strata_by))
      stop("'add_inferential_pvalues' cannot be TRUE when 'strata_by' is present.")
  }

  # ---------------------------------------------------------------------------
  # 7. Resolve summarize_what
  # ---------------------------------------------------------------------------
  if (is.null(summarize_what)) {
    summarize_what <- setdiff(names(x), c(split_by, strata_by))
  } else {
    if (!is.character(summarize_what))
      stop("'summarize_what' must be a character vector of column names or NULL.")
    missing_cols <- setdiff(summarize_what, names(x))
    if (length(missing_cols) > 0L)
      stop(paste0("Columns in 'summarize_what' not found in data: ",
                  paste(missing_cols, collapse = ", "), "."))
  }
  if (length(summarize_what) == 0L)
    stop("No columns remain to summarise.")

  # ---------------------------------------------------------------------------
  # 8. Remaining parameter validation
  # ---------------------------------------------------------------------------
  if (!is.logical(include_missing_in_splits))
    stop("'include_missing_in_splits' must be TRUE or FALSE.")

  valid_calc_col <- c("n_in_column", "n_valid_in_column")
  if (!calc_col_percent_using %in% valid_calc_col)
    stop(paste0("'calc_col_percent_using' must be one of: ",
                paste(valid_calc_col, collapse = ", "), "."))

  valid_cont_stats <- c("meanSD", "meanSD2", "medianIQR", "mean", "sd",
                        "median", "p0", "p25", "p50", "p75", "p100", "IQR")
  if (!all(continuous_statistics %in% valid_cont_stats))
    stop(paste0("'continuous_statistics' must be one of: ",
                paste(valid_cont_stats, collapse = ", "), "."))

  valid_cat_stats <- c("n_percent", "n", "percent")
  if (!categorical_statistics %in% valid_cat_stats)
    stop(paste0("'categorical_statistics' must be one of: ",
                paste(valid_cat_stats, collapse = ", "), "."))

  if (!is.null(force_continuous)) {
    if (!is.character(force_continuous) ||
        !all(force_continuous %in% names(x)))
      stop("'force_continuous' must be a character vector of valid column names.")
  }
  if (!is.null(force_categorical)) {
    if (!is.character(force_categorical) ||
        !all(force_categorical %in% names(x)))
      stop("'force_categorical' must be a character vector of valid column names.")
  }
  if (!is.null(force_continuous) && !is.null(force_categorical)) {
    overlap <- intersect(force_continuous, force_categorical)
    if (length(overlap) > 0)
      stop(paste0("Variables in both 'force_continuous' and 'force_categorical': ",
                  paste(overlap, collapse = ", "), "."))
  }

  # Validate n_digits_continuous
  if (!(
    (is.numeric(n_digits_continuous) && length(n_digits_continuous) == 2 &&
     all(n_digits_continuous >= 0)) ||
    (is.character(n_digits_continuous) && length(n_digits_continuous) == 1 &&
     n_digits_continuous == "dynamic")
  ))
    stop("'n_digits_continuous' must be c(d1, d2) or \"dynamic\".")

  # Validate n_digits_categorical
  if (!(
    (is.numeric(n_digits_categorical) && length(n_digits_categorical) == 2 &&
     all(n_digits_categorical >= 0)) ||
    (is.character(n_digits_categorical) && length(n_digits_categorical) == 1 &&
     n_digits_categorical == "dynamic")
  ))
    stop("'n_digits_categorical' must be c(d1, d2) or \"dynamic\".")

  if (!is.logical(zero_as_exp) || length(zero_as_exp) != 1)
    stop("'zero_as_exp' must be TRUE or FALSE.")

  valid_display <- c("ifany", "no", "always")
  if (!display_missing %in% valid_display)
    stop(paste0("'display_missing' must be one of: ",
                paste(valid_display, collapse = ", "), "."))

  valid_missing_stat <- c("n_percent", "n", "percent")
  if (!missing_stat %in% valid_missing_stat)
    stop(paste0("'missing_stat' must be one of: ",
                paste(valid_missing_stat, collapse = ", "), "."))

  valid_sort <- c("alphanumeric", "frequency")
  if (!sort_categorical_variables_by %in% valid_sort)
    stop(paste0("'sort_categorical_variables_by' must be one of: ",
                paste(valid_sort, collapse = ", "), "."))

  valid_pct <- c("column", "row", "cell", "total")
  if (!calc_percent_by %in% valid_pct)
    stop(paste0("'calc_percent_by' must be one of: ",
                paste(valid_pct, collapse = ", "), "."))

  if (!is.logical(add_inferential_pvalues))
    stop("'add_inferential_pvalues' must be TRUE or FALSE.")

  if (!is.logical(paired) || length(paired) != 1)
    stop("'paired' must be a single logical value (TRUE or FALSE).")

  if (!is.null(n_digits_pvalues) &&
      (!is.numeric(n_digits_pvalues) || n_digits_pvalues <= 0))
    stop("'n_digits_pvalues' must be a positive numeric value or NULL.")

  if (!is.logical(bold_significant_pvalues))
    stop("'bold_significant_pvalues' must be TRUE or FALSE.")
  if (!is.numeric(bold_significant_pvalues_at))
    stop("'bold_significant_pvalues_at' must be numeric.")
  if (!is.logical(show_n_header) || length(show_n_header) != 1)
    stop("'show_n_header' must be TRUE or FALSE.")
  if (!is.logical(bold_labels) || length(bold_labels) != 1)
    stop("'bold_labels' must be TRUE or FALSE.")
  if (!is.logical(italicize_levels) || length(italicize_levels) != 1)
    stop("'italicize_levels' must be TRUE or FALSE.")
  if (!is.logical(clean_table) || length(clean_table) != 1)
    stop("'clean_table' must be TRUE or FALSE.")
  if (!is.character(table_name) && !is.null(table_name))
    stop("'table_name' must be a string or NULL.")

  if (!is.null(rename_variables)) {
    if (!is.list(rename_variables) ||
        any(!sapply(rename_variables, inherits, "formula")))
      stop("'rename_variables' must be a list of formulas.")
    orig_vars <- sapply(rename_variables, function(f) as.character(f)[[2]])
    missing_rv <- orig_vars[!orig_vars %in% names(x)]
    if (length(missing_rv) > 0)
      stop(paste0("Variables in 'rename_variables' not in data: ",
                  paste(missing_rv, collapse = ", "), "."))
  }

  # ---------------------------------------------------------------------------
  # 9. Coerce force_continuous columns to numeric
  # ---------------------------------------------------------------------------
  if (!is.null(force_continuous)) {
    for (col_name in force_continuous) {
      if (!is.numeric(x[[col_name]])) {
        converted <- suppressWarnings(as.numeric(x[[col_name]]))
        if (any(is.na(converted) & !is.na(x[[col_name]])))
          warning(paste0("Column '", col_name,
                         "' contains non-numeric values coerced to NA."))
        x[[col_name]] <- converted
      }
    }
  }

  # Keep a clean copy of the data for inferential tests (run_diff/run_assoc)
  # so that missing values transformed into levels for % calculation aren't
  # treated as real categories in statistical tests.
  x_clean <- x

  # ---------------------------------------------------------------------------
  # 10. Missing-value handling for total/column denominators
  # ---------------------------------------------------------------------------
  use_total_denom <- (calc_percent_by == "total") || 
                     (calc_percent_by == "column" && calc_col_percent_using == "n_in_column")
  
  # Handle warning for inferential p-values with total column denominators
  if (use_total_denom && add_inferential_pvalues && calc_percent_by == "column") {
    warning("`calc_col_percent_using = 'n_in_column'` is ignored when `add_inferential_pvalues = TRUE`.")
    calc_col_percent_using <- "n_valid_in_column"
    use_total_denom <- (calc_percent_by == "total") # Recalculate flag
  }

  if (use_total_denom) {
    # Identify categorical columns to convert NAs to levels.
    # This ensures gtsummary includes them in the denominator (total N).
    is_cat_check <- function(v_name) {
      col <- x[[v_name]]
      if (v_name %in% force_continuous) return(FALSE)
      if (v_name %in% force_categorical) return(TRUE)
      if (is.character(col) || is.factor(col)) return(TRUE)
      # Heuristic for numeric columns gtsummary would treat as categorical
      if (is.numeric(col) && length(unique(stats::na.omit(col))) < 10) return(TRUE)
      FALSE
    }
    
    to_convert <- summarize_what[vapply(summarize_what, is_cat_check, logical(1))]
    
    # Also include grouping variables if missing values should be explicit columns
    if (include_missing_in_splits) {
      to_convert <- unique(c(to_convert, split_by, strata_by)) |> purrr::compact()
    }
    
    if (length(to_convert) > 0) {
      x <- x |>
        dplyr::mutate(dplyr::across(
          dplyr::all_of(to_convert),
          ~ forcats::fct_na_value_to_level(factor(.), level = missing_text)
        ))
    }
  } else if (include_missing_in_splits) {
    # Only grouping variables need conversion if they aren't using total denom
    to_convert <- unique(c(split_by, strata_by)) |> purrr::compact()
    if (length(to_convert) > 0) {
      x <- x |>
        dplyr::mutate(dplyr::across(
          dplyr::all_of(to_convert),
          ~ forcats::fct_na_value_to_level(factor(.), level = missing_text)
        ))
    }
  }

  # ---------------------------------------------------------------------------
  # 11. gtsummary stat strings
  # ---------------------------------------------------------------------------
  stat_exprs <- c(
    meanSD    = "{mean} \u00b1 {sd}",
    meanSD2   = "{mean} ({sd})",
    medianIQR = "{median} ({p25}, {p75})",
    mean      = "{mean}",
    sd        = "{sd}",
    median    = "{median}",
    p0        = "{p0}",
    p25       = "{p25}",
    p50       = "{p50}",
    p75       = "{p75}",
    p100      = "{p100}",
    IQR       = "{p25}, {p75}"
  )
  gts_continuous_stat  <- paste(stat_exprs[continuous_statistics],
                                collapse = ", ")
  gts_categorical_stat <- switch(categorical_statistics,
    n_percent = "{n} ({p}%)",
    n         = "{n}",
    percent   = "{p}%"
  )
  gts_missing_stat <- switch(missing_stat,
    n_percent = "{N_miss} ({p_miss}%)",
    n         = "{N_miss}",
    percent   = "{p_miss}%"
  )

  # ---------------------------------------------------------------------------
  # 12. type_list for gtsummary
  # ---------------------------------------------------------------------------
  type_list <- if (!is.null(force_continuous)) {
    lapply(force_continuous, function(col_name)
      stats::as.formula(paste0("`", col_name, "` ~ 'continuous'")))
  } else list()

  if (!is.null(force_categorical)) {
    type_list <- c(type_list, lapply(force_categorical, function(col_name)
      stats::as.formula(paste0("`", col_name, "` ~ 'categorical'"))))
  }

  # ---------------------------------------------------------------------------
  # 13. zero_as_exp helpers
  # ---------------------------------------------------------------------------
  .maybe_exp <- function(token, dig) {
    val <- suppressWarnings(as.numeric(token))
    if (is.na(val) || val == 0) return(token)
    rounded_val <- as.numeric(formatC(round(val, dig), digits = dig, format = "f"))
    if (rounded_val == 0) {
      sci <- formatC(val, digits = dig, format = "e")
      sci <- gsub("e\\+0*(\\d+)", "e\\1", sci)
      sci <- gsub("e-0*(\\d+)",  "e-\\1", sci)
      return(sci)
    }
    token
  }

  .apply_exp_to_cell <- function(cell_str, dig1, dig2) {
    if (is.na(cell_str) || !nzchar(cell_str)) return(cell_str)
    tokens <- strsplit(
      cell_str,
      "(?<=\\d)(?=[^0-9eE.\\-+])|(?<=[^0-9eE.\\-+])(?=[0-9\\-])",
      perl = TRUE
    )[[1]]
    numeric_idx <- which(!is.na(suppressWarnings(as.numeric(tokens))))
    for (i in seq_along(numeric_idx)) {
      dig <- if (i == 1L) dig1 else dig2
      tokens[numeric_idx[i]] <- .maybe_exp(tokens[numeric_idx[i]], dig)
    }
    paste(tokens, collapse = "")
  }

  # ---------------------------------------------------------------------------
  # 14. Continuous digits argument for gtsummary
  # ---------------------------------------------------------------------------
  n_stats_required <- switch(continuous_statistics,
    meanSD = 2, meanSD2 = 2, medianIQR = 3,
    mean = 1, sd = 1, median = 1,
    p0 = 1, p25 = 1, p50 = 1, p75 = 1, p100 = 1,
    IQR = 2, 1
  )

  if (identical(n_digits_continuous, "dynamic")) {
    continuous_digits_arg <- NULL
  } else {
    actual_digits <- rep(n_digits_continuous, length.out = n_stats_required)
    if (zero_as_exp) {
      make_exp_formatter <- function(dig) {
        function(x_val) {
          purrr::map_chr(x_val, function(v) {
            if (is.na(v))   return(NA_character_)
            if (v == 0)     return(formatC(0, format = "f", digits = dig))
            rv <- round(v, dig)
            if (rv == 0 && v != 0) {
              sci <- formatC(v, digits = dig, format = "e")
              sci <- sub("e\\+0*", "e",  sci)
              sci <- sub("e-0*",   "e-", sci)
              return(sci)
            }
            formatC(v, format = "f", digits = dig)
          })
        }
      }
      continuous_digits_arg <- lapply(actual_digits, make_exp_formatter)
      if (length(continuous_digits_arg) == 1)
        continuous_digits_arg <- continuous_digits_arg[[1]]
    } else {
      continuous_digits_arg <- actual_digits
    }
  }

  if (identical(n_digits_categorical, "dynamic"))
    n_digits_categorical <- NULL

  # ---------------------------------------------------------------------------
  # 15. Internal summary-table builder
  # ---------------------------------------------------------------------------
  create_summary <- function(data) {
    gtsummary::tbl_summary(
      data         = data,
      include      = summarize_what,
      by           = split_by,
      label        = rename_variables,
      type         = type_list,
      statistic    = list(
        gtsummary::all_continuous()  ~ gts_continuous_stat,
        gtsummary::all_categorical() ~ gts_categorical_stat
      ),
      digits       = list(
        gtsummary::all_categorical() ~ n_digits_categorical,
        gtsummary::all_continuous()  ~ continuous_digits_arg
      ),
      missing      = display_missing,
      missing_text = missing_text,
      missing_stat = gts_missing_stat,
      sort         = list(gtsummary::all_categorical() ~
                            sort_categorical_variables_by),
      percent      = ifelse(calc_percent_by == "total", "cell", calc_percent_by)
    ) |>
      gtsummary::modify_header(label = paste0("**", header, "**"))
  }

  # ---------------------------------------------------------------------------
  # 16. Build base gtsummary table
  # ---------------------------------------------------------------------------
  result_table <- if (!is.null(strata_by)) {
    x |>
      dplyr::filter(!is.na(.data[[strata_by]])) |>
      gtsummary::tbl_strata(
        strata    = strata_by,
        include   = summarize_what,
        .tbl_fun  = \(.x) create_summary(.x)
      )
  } else {
    create_summary(x)
  }

  # Spanning header for split_by
  if (is.null(strata_by) && !is.null(split_by) && !is.null(split_by_header)) {
    result_table <- result_table |>
      gtsummary::modify_spanning_header(
        gtsummary::all_stat_cols() ~
          paste0("**", split_by_header, "**")
      )
  }

  # ---------------------------------------------------------------------------
  # 18.5 Re-calculate missing row percentages
  # ---------------------------------------------------------------------------
  # Ensure the 'No data/missing' row respects calc_percent_by and 
  # calc_col_percent_using. This is primarily for continuous variables or 
  # categorical variables where NAs were not converted to levels.
  if (display_missing != "no" && is.null(strata_by)) {
    total_n <- dim(x)[1]
    
    # Identify groups and their sizes
    if (!is.null(split_by)) {
      # Respect gtsummary ordering (factor levels or frequency)
      grp_factor <- factor(x[[split_by]])
      if (sort_categorical_variables_by == "frequency") {
        grp_factor <- forcats::fct_infreq(grp_factor)
      }
      grp_lvls <- levels(grp_factor)
      grp_ns   <- as.numeric(table(grp_factor)[grp_lvls])
    } else {
      grp_lvls <- NULL
      grp_ns   <- total_n
    }

    # Pre-calculate missing counts per variable per group
    vars_with_missing <- unique(result_table$table_body$variable[result_table$table_body$row_type == "missing"])
    
    miss_counts_map <- list()
    for (v in vars_with_missing) {
      if (is.null(split_by)) {
        miss_counts_map[[v]] <- sum(is.na(x[[v]]))
      } else {
        counts <- tapply(is.na(x[[v]]), grp_factor, sum)
        # Handle potential NA groups or missing levels
        counts[is.na(counts)] <- 0
        miss_counts_map[[v]] <- counts[grp_lvls]
      }
    }

    result_table <- result_table |>
      gtsummary::modify_table_body(function(tb) {
        stat_cols <- names(tb)[grep("^stat_", names(tb))]
        
        for (col_nm in stat_cols) {
          # Map stat_X to group index
          idx <- if (is.null(split_by)) 1 else as.numeric(sub("stat_", "", col_nm))
          if (is.na(idx) || idx > length(grp_ns)) next
          
          denom_group <- grp_ns[idx]
          
          tb[[col_nm]] <- purrr::imap_chr(tb[[col_nm]], function(cell, i) {
            if (is.na(cell) || !nzchar(cell) || tb$row_type[i] != "missing") return(cell)
            
            var_nm <- tb$variable[i]
            n_miss <- if (is.null(split_by)) miss_counts_map[[var_nm]] else miss_counts_map[[var_nm]][idx]
            if (is.na(n_miss)) n_miss <- 0
            
            # Determine denominator based on parameters
            final_denom <- if (calc_percent_by == "total") {
              total_n
            } else if (calc_percent_by == "column") {
              if (calc_col_percent_using == "n_in_column") denom_group else denom_group - n_miss
            } else {
              denom_group
            }
            
            if (is.na(final_denom) || final_denom <= 0) return(cell)
            
            p_val <- (n_miss / final_denom) * 100
            dig   <- if (!is.null(n_digits_categorical)) n_digits_categorical[2] else 2
            p_str <- formatC(p_val, digits = dig, format = "f")
            
            # Reconstruct string based on existing content
            if (grepl("\\(", cell)) {
              return(sprintf("%d (%s%%)", n_miss, p_str))
            } else if (grepl("%", cell)) {
              return(sprintf("%s%%", p_str))
            } else {
              return(as.character(n_miss))
            }
          })
        }
        tb
      })
  }

  # bold_labels / italicize_levels
  if (bold_labels)      result_table <- result_table |> gtsummary::bold_labels()
  if (italicize_levels) result_table <- result_table |> gtsummary::italicize_levels()

  # ---------------------------------------------------------------------------
  # 17. clean_table
  # ---------------------------------------------------------------------------
  if (clean_table) {
    zero_pct_fmt <- paste0(
      "0",
      if (!is.null(n_digits_categorical) && n_digits_categorical[2] > 0)
        paste0(".", strrep("0", n_digits_categorical[2]))
      else "",
      "%"
    )
    zero_strings <- switch(categorical_statistics,
      n_percent = paste0("0 (", zero_pct_fmt, ")"),
      n         = "0",
      percent   = zero_pct_fmt
    )
    if (missing_stat == "n_percent")
      zero_strings <- c(zero_strings, "0 (0.0%)", "0 (0%)", "0 (NA%)",
                        "NA \u00b1 NA")
    zero_strings <- unique(zero_strings)

    result_table <- result_table |>
      gtsummary::modify_table_body(
        ~ .x |>
          dplyr::mutate(dplyr::across(
            dplyr::starts_with("stat_"),
            ~ dplyr::case_when(.x %in% zero_strings ~ "", TRUE ~ .x)
          ))
      )
  }

  # ---------------------------------------------------------------------------
  # 18. zero_as_exp post-processing
  # ---------------------------------------------------------------------------
  if (zero_as_exp && !is.null(n_digits_continuous)) {
    dig1 <- n_digits_continuous[1]
    dig2 <- n_digits_continuous[2]
    forced_cont_vars <- if (!is.null(force_continuous)) force_continuous else character(0)

    result_table <- result_table |>
      gtsummary::modify_table_body(function(tb) {
        if (!all(c("var_type", "variable", "row_type") %in% names(tb)))
          return(tb)
        tb |>
          dplyr::mutate(dplyr::across(
            dplyr::starts_with("stat_"),
            ~ purrr::imap_chr(.x, function(cell, i) {
              if (is.na(cell) || !nzchar(trimws(cell))) return(cell)
              is_cont <- isTRUE(tb$var_type[i] == "continuous") ||
                         isTRUE(tb$variable[i] %in% forced_cont_vars)
              if (is_cont && isTRUE(tb$row_type[i] == "label"))
                .apply_exp_to_cell(cell, dig1, dig2)
              else
                cell
            })
          ))
      })
  }

  # ---------------------------------------------------------------------------
  # 19. Inferential p-values via run_diff / run_assoc
  # ---------------------------------------------------------------------------
  if (add_inferential_pvalues) {

    # -- 19a. Determine variable type for each summarised variable ------------
    # Pull the table body to inspect var_type assigned by gtsummary
    tb_body <- result_table$table_body

    # Helper: is a variable treated as continuous?
    .is_continuous <- function(var_nm) {
      if (!is.null(force_continuous) && var_nm %in% force_continuous) return(TRUE)
      if (!is.null(force_categorical) && var_nm %in% force_categorical) return(FALSE)
      rows <- tb_body[tb_body$variable == var_nm, ]
      if (dim(rows)[1] == 0) return(FALSE)
      isTRUE(rows$var_type[1] == "continuous")
    }

    # Only label rows carry the p-value; collect unique variable names
    label_vars <- unique(tb_body$variable[tb_body$row_type == "label"])
    label_vars <- intersect(label_vars, summarize_what)

    # -- 19b. Format p-value string ------------------------------------------
    .fmt_p <- function(p) {
      if (is.null(n_digits_pvalues)) {
        if (is.na(p))  return(NA_character_)
        if (p < 0.001) return("<0.001")
        if (p < 0.10)  return(formatC(p, digits = 3, format = "f"))
        return(formatC(p, digits = 1, format = "f"))
      }
      if (is.na(p))  return(NA_character_)
      if (p < 0.001) return("<0.001")
      formatC(p, digits = n_digits_pvalues, format = "f")
    }

    # -- 19c. Run engines and collect results ---------------------------------
    # We build two lookup tables keyed by variable name:
    #   p_map      : variable -> formatted p-value string (or NA)
    #   test_map   : variable -> test name string (or NA)
    #   error_map  : variable -> TRUE if engine failed

    p_map     <- stats::setNames(rep(NA_character_, length(label_vars)),
                                 label_vars)
    test_map  <- stats::setNames(rep(NA_character_, length(label_vars)),
                                 label_vars)
    error_map <- stats::setNames(rep(FALSE,          length(label_vars)),
                                 label_vars)

    for (var_nm in label_vars) {

      if (.is_continuous(var_nm)) {
        # ---- run_diff -------------------------------------------------------
        engine_result <- tryCatch({
          run_diff(
            x                     = x_clean[, var_nm, drop = FALSE],
            metadata              = x_clean,
            outcome               = var_nm,
            group                 = split_by,
            filter                = NULL, # already filtered globally
            paired                = paired,
            test_type             = test_type_continuous,
            verbose               = FALSE
          )
        }, error = function(e) {
          message(sprintf(
            "[run_summarytable] run_diff failed for '%s': %s",
            var_nm, conditionMessage(e)
          ))
          NULL
        })

        if (!is.null(engine_result) &&
            !is.null(engine_result$test_result) &&
            !is.na(engine_result$test_result$p.value)) {
          p_map[var_nm]   <- .fmt_p(engine_result$test_result$p.value)
          test_map[var_nm] <- engine_result$test_used
        } else {
          error_map[var_nm] <- TRUE
        }

      } else {
        # ---- run_assoc -------------------------------------------------------
        # Ensure the variable is factor/character before calling run_assoc
        x_assoc <- x_clean
        if (!is.factor(x_assoc[[var_nm]]) &&
            !is.character(x_assoc[[var_nm]])) {
          x_assoc[[var_nm]] <- as.factor(x_assoc[[var_nm]])
        }

        engine_result <- tryCatch({
          run_assoc(
            x                 = x_assoc,
            var1              = var_nm,
            var2              = split_by,
            paired            = paired,
            test_type         = test_type_categorical,
            verbose           = FALSE
          )
        }, error = function(e) {
          message(sprintf(
            "[run_summarytable] run_assoc failed for '%s': %s",
            var_nm, conditionMessage(e)
          ))
          NULL
        })

        if (!is.null(engine_result) &&
            !is.null(engine_result$test_result) &&
            !is.na(engine_result$test_result$p.value)) {
          p_map[var_nm]    <- .fmt_p(engine_result$test_result$p.value)
          test_map[var_nm] <- engine_result$test_used
        } else {
          error_map[var_nm] <- TRUE
        }
      }
    }

    # -- 19d. Build superscript legend ----------------------------------------
    # Assign a letter (a, b, c, …) to each unique test name, in order of
    # first appearance across label_vars.
    distinct_tests  <- unique(test_map[!is.na(test_map)])
    test_letters    <- stats::setNames(
      letters[seq_along(distinct_tests)],
      distinct_tests
    )
    # Map variable -> letter (NA for failed)
    var_letter <- stats::setNames(
      vapply(label_vars, function(v) {
        tn <- test_map[v]
        if (is.na(tn)) NA_character_ else test_letters[[tn]]
      }, character(1)),
      label_vars
    )

    any_failed <- any(error_map)

    # -- 19e. Inject p.value column into table body ---------------------------
    result_table <- result_table |>
      gtsummary::modify_table_body(function(tb) {
        # Add p.value column if absent (Option B — we never called add_p())
        if (!"p.value" %in% names(tb))
          tb$p.value <- NA_character_

        for (i in seq_len(dim(tb)[1])) {
          var_nm   <- tb$variable[i]
          row_type <- tb$row_type[i]
          if (row_type != "label") next
          if (!var_nm %in% label_vars) next

          if (error_map[var_nm]) {
            # em-dash + dagger superscript
            tb$p.value[i] <- "\u2014\u2020"
          } else {
            pv  <- p_map[var_nm]
            ltr <- var_letter[var_nm]
            # Bold wrapper (we apply bold later via modify_table_body;
            # store raw value + superscript here)
            tb$p.value[i] <- if (!is.na(ltr))
              paste0(pv, " (", ltr, ")")
            else
              pv
          }
        }
        tb
      })

    # -- 19f. Modify p.value column header ------------------------------------
    result_table <- result_table |>
      gtsummary::modify_header(p.value = "**p-value**")

    # -- 19g. Convert to gt --------------------------------------------------
    result_table <- result_table |> gtsummary::as_gt()

    # Reformat p-value cells: "0.043 (a)" -> "0.043<sup>a</sup>"
    result_table <- result_table |>
      gt::text_transform(
        locations = gt::cells_body(columns = "p.value"),
        fn = function(x_cell) {
          gsub(
            "([^\\s]+) \\(([a-z])\\)",
            "\\1<sup>\\2</sup>",
            x_cell,
            perl = TRUE
          )
        }
      )

    # Bold significant p-values if requested
    if (bold_significant_pvalues) {
      num_p_map <- stats::setNames(
        vapply(label_vars, function(v) {
          pv_str <- p_map[v]
          if (is.na(pv_str)) return(NA_real_)
          if (pv_str == "<0.001") return(0.0005)
          suppressWarnings(as.numeric(pv_str))
        }, numeric(1)),
        label_vars
      )

      tb_body2 <- result_table[["_data"]]
      if (!is.null(tb_body2)) {
        sig_rows <- which(vapply(seq_len(dim(tb_body2)[1]), function(i) {
          var_nm <- tb_body2$variable[i]
          if (is.null(var_nm) || !var_nm %in% label_vars) return(FALSE)
          pv <- num_p_map[var_nm]
          if (is.na(pv)) return(FALSE)
          pv < bold_significant_pvalues_at
        }, logical(1)))

        if (length(sig_rows) > 0) {
          result_table <- result_table |>
            gt::tab_style(
              style     = gt::cell_text(weight = "bold"),
              locations = gt::cells_body(
                columns = "p.value",
                rows    = sig_rows
              )
            )
        }
      }
    }

    # -- 19h. Add footnotes ---------------------------------------------------
    # Unicode superscript letters a-z for use in footnote labels.
    # This avoids gt adding its own numeric counter alongside our letters.
    superscript_map <- c(
      a = "\u1d43", b = "\u1d47", c = "\u1d9c", d = "\u1d48",
      e = "\u1d49", f = "\u1da0", g = "\u1d4d", h = "\u02b0",
      i = "\u2071",  j = "\u02b2", k = "\u1d4f", l = "\u02e1",
      m = "\u1d50", n = "\u207f", o = "\u1d52", p = "\u1d56",
      r = "\u02b3", s = "\u02e2", t = "\u1d57", u = "\u1d58",
      v = "\u1d5b", w = "\u02b7", x = "\u02e3", y = "\u02b8",
      z = "\u1dbb"
    )

    # Remove gt's automatic footnote numbering from the p.value column label
    # by suppressing source notes added via tab_footnote on column labels.
    # We instead append a single clean footnote block as a source note so
    # gt does not inject superscript numbers into the column header.
    footnote_lines <- character(0)
    for (test_nm in names(test_letters)) {
      ltr     <- test_letters[[test_nm]]
      sup_chr <- superscript_map[ltr]
      if (is.na(sup_chr)) sup_chr <- paste0("(", ltr, ")")
      footnote_lines <- c(footnote_lines, paste0(sup_chr, " ", test_nm))
    }
    if (any_failed) {
      footnote_lines <- c(footnote_lines,
        "\u2020 P-value could not be computed for this variable.")
    }

    # Append all footnote lines as a single source note so gt renders them
    # cleanly without adding numeric superscripts to the column header.
    result_table <- result_table |>
      gt::tab_source_note(
        source_note = gt::html(paste(footnote_lines, collapse = "<br>"))
      )

  } else {
    # No inferential p-values: just convert to gt for title addition below
    result_table <- result_table |> gtsummary::as_gt()
  }

  # ---------------------------------------------------------------------------
  # 20. Table title
  # ---------------------------------------------------------------------------
  if (!requireNamespace("gt", quietly = TRUE))
    stop("Package 'gt' is required. Install with: install.packages('gt')")

  if (!is.null(table_name)) {
    title_text <- if (identical(table_name, "auto")) {
      prefix <- if (add_inferential_pvalues)
        "Inferential Statistics" else "Descriptive Statistics"
      if (!is.null(strata_by) && !is.null(split_by_header)) {
        paste0(prefix, " Stratified by ", strata_by,
               " and ", split_by_header, " Columns")
      } else if (is.null(strata_by) && !is.null(split_by_header)) {
        paste0(prefix, " Stratified by ", split_by_header, " Column")
      } else {
        prefix
      }
    } else {
      table_name
    }

    result_table <- result_table |>
      gt::tab_header(title = gt::md(paste0("**", title_text, "**")))
  }

  # ---------------------------------------------------------------------------
  # 21. Return
  # ---------------------------------------------------------------------------

  result_table
}