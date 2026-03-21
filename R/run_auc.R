#' @title Comprehensive Area Under the ROC Curve (AUROC) Analysis
#'
#' @description
#' Performs a comprehensive Receiver Operating Characteristic (ROC) analysis
#' for one or more continuous predictors against a binary response variable.
#' Wraps \code{pROC::roc()} to compute AUC and confidence intervals, and
#' optionally produces publication-ready ROC plots via \pkg{ggplot2}.
#'
#' @param x A \code{data.frame}, \code{tibble}, or named \code{matrix} whose
#'   columns are the candidate predictor variables to be evaluated.
#' @param y A single character string naming the column in \code{x} that
#'   contains the binary grouping variable (cases vs. controls). The column
#'   may be categorical (\code{character} / \code{factor}) with exactly two
#'   levels after optional removal via \code{remove}, or numeric containing
#'   only \code{0} (control) and \code{1} (case) after optional removal.
#' @param remove Optional. A vector of values to remove from \code{y} before
#'   analysis, used to reduce a multi-level variable down to exactly two
#'   groups. May be a character vector (for categorical \code{y}), a numeric
#'   vector (for numeric \code{y}), or a mixed list. Ignored when \code{y}
#'   already contains exactly two categories or only \code{0}/\code{1} values.
#'   Default \code{NULL}.
#' @param metadata An optional object of the same class and row count as
#'   \code{x} carrying additional sample annotations. Currently reserved for
#'   future use and is ignored by the analysis. Default \code{NULL}.
#' @param set_control A single string specifying which level of a categorical
#'   \code{y} represents the control group. Ignored when \code{y} is numeric.
#'   Default \code{NULL} (the first level alphabetically or by factor order).
#' @param include_cols A character vector of column names, or an integer
#'   vector of column positions, selecting the predictor columns of \code{x}
#'   to be analysed. Default \code{NULL} includes all non-\code{y} columns.
#' @param include_rows A character vector of row names, or an integer vector
#'   of row positions, selecting the rows of \code{x} to be analysed. Row
#'   names are only considered when \code{rownames(x)} are not the default
#'   integer sequence. Default \code{NULL} includes all rows.
#' @param direction A single string controlling the direction of comparison
#'   passed to \code{pROC::roc()}. Choices: \code{"auto"} (default),
#'   \code{">"} (controls > cases), or \code{"<"} (controls < cases). See
#'   \code{?pROC::roc} for details.
#' @param compute_auc Logical. Whether to compute the AUC. Default \code{TRUE}.
#' @param compute_ci Logical. Whether to compute bootstrap or DeLong
#'   confidence intervals for the AUC. Default \code{TRUE}.
#' @param ci_level A numeric scalar in \eqn{(0, 1)} specifying the confidence
#'   level. Default \code{0.95}.
#' @param plot Controls ROC plot generation. Accepted values:
#'   \describe{
#'     \item{\code{"none"}}{No plots are produced.}
#'     \item{\code{"all"}}{All selected predictors are overlaid on one plot.}
#'     \item{\code{"separate"}}{One plot per predictor, returned as a named
#'       list.}
#'     \item{Integer \eqn{\ge 1}}{Plot the top \emph{n} predictors (by AUC)
#'       on one combined chart.}
#'     \item{Numeric in \eqn{[0, 1)}}{Plot only predictors with
#'       \eqn{AUC \ge} this threshold.}
#'     \item{\code{c("all", n)} or \code{c("separate", n)}}{Combine a string
#'       mode with an integer or numeric filter as described above.}
#'   }
#'   Default \code{"none"}.
#' @param theme A single string specifying the \pkg{ggplot2} theme for plots.
#'   Supported values: \code{"nature"} (default), \code{"minimal"},
#'   \code{"bw"}, \code{"classic"}, \code{"plos"}, \code{"gray"},
#'   \code{"light"}, \code{"dark"}. Colour-blind-friendly palettes are applied
#'   automatically.
#' @param ncol Integer. Number of legend columns when \code{plot = "all"} or
#'   \code{plot = c("all", ...)}. Default \code{NULL} (auto).
#' @param nrow Integer. Number of legend rows when \code{plot = "all"} or
#'   \code{plot = c("all", ...)}. Default \code{NULL} (auto).
#' @param plot_title A character vector of plot titles. \code{NULL} (default)
#'   uses \code{"ROC Curve"} for combined plots and
#'   \code{"ROC Curve for <predictor>"} for separate plots. When
#'   \code{plot = "separate"}, the vector must match the number of predictors
#'   being plotted.
#' @param plot_subtitle A single string used as the plot subtitle. Default
#'   \code{NULL} (no subtitle).
#' @param plot_legend A single string for the legend title. Default
#'   \code{"AUC (95\% CI)"}.
#' @param xlab A single string for the x-axis label. Default \code{NULL}
#'   renders \code{"1 - Specificity (False Positive Rate)"}.
#' @param ylab A single string for the y-axis label. Default \code{NULL}
#'   renders \code{"Sensitivity (True Positive Rate)"}.
#' @param linewidth A positive numeric value controlling the ROC curve line
#'   width. Default \code{1}.
#' @param global_font_size A positive numeric value setting the base font size
#'   for all plot text elements. Default \code{20}. Derived sizes:
#'   title \eqn{= 1.20 \times}, subtitle \eqn{= 0.90 \times},
#'   axis labels / ticks \eqn{= 1.00 \times}.
#' @param title_font_size Override for the title font size. Default \code{NULL}
#'   (derived from \code{global_font_size}).
#' @param subtitle_font_size Override for the subtitle font size. Default
#'   \code{NULL} (derived from \code{global_font_size}).
#' @param xlab_font_size Override for the x-axis label and tick font size.
#'   Default \code{NULL} (derived from \code{global_font_size}).
#' @param ylab_font_size Override for the y-axis label and tick font size.
#'   Default \code{NULL} (derived from \code{global_font_size}).
#' @param inplot_font_size Font size of the in-plot legend text (e.g., AUC
#'   labels in combined plots). Default \code{NULL} (derived from
#'   \code{global_font_size}).
#' @param ... Additional arguments forwarded to \code{pROC::roc()},
#'   \code{pROC::auc()}, or \code{pROC::ci()}.
#'
#' @details
#' \strong{Workflow overview}\cr
#' \enumerate{
#'   \item \code{x} is coerced to a \code{data.frame}; columns and rows are
#'     subset according to \code{include_cols} / \code{include_rows}.
#'   \item \code{y} values listed in \code{remove} are dropped, leaving
#'     exactly two groups.
#'   \item When \code{y} is categorical, the two remaining levels are mapped to
#'     \code{0} (control, as specified by \code{set_control}) and \code{1}
#'     (case). If \code{y} already encodes \code{0}/\code{1} as a character or
#'     factor, it is silently converted to numeric.
#'   \item \code{pROC::roc()} is called for each predictor column; AUC and,
#'     optionally, DeLong confidence intervals (\code{pROC::ci()}) are
#'     extracted.
#'   \item Plots are assembled with \pkg{ggplot2} using a colour-blind-friendly
#'     palette (Okabe–Ito). AUC labels are rounded to two decimal places;
#'     values that round to \code{1.00} are progressively displayed to more
#'     decimal places (up to five), after which scientific notation with two
#'     significant figures is used.
#' }
#'
#' \strong{Reference line}\cr
#' A grey dashed diagonal line representing an AUC of 0.5 (random classifier)
#' is always drawn.
#'
#' \strong{Legend placement}\cr
#' For combined plots (\code{plot = "all"}) the legend is positioned at the
#' bottom-right inside the plot area by default.
#'
#' \strong{pROC dependency}\cr
#' This function requires \pkg{pROC} (\eqn{\ge} 1.18) and \pkg{ggplot2}
#' (\eqn{\ge} 3.4). Both are listed as \code{Imports} in the package
#' \code{DESCRIPTION} file.
#'
#' @return A named list of class \code{"run_auc"} with the following elements:
#' \describe{
#'   \item{\code{roc_list}}{Named list of \code{pROC::roc} objects, one per
#'     analysed predictor.}
#'   \item{\code{auc_table}}{A \code{data.frame} with columns
#'     \code{predictor}, \code{auc}, \code{ci_lower}, \code{ci_upper},
#'     \code{ci_level}, \code{direction}, and \code{n_cases} /
#'     \code{n_controls}.}
#'   \item{\code{plots}}{A named list of \code{ggplot} objects (or
#'     \code{NULL} when \code{plot = "none"}). The element \code{"combined"}
#'     is present for combined-mode plots; individual predictor names are used
#'     as keys for separate-mode plots.}
#'   \item{\code{params}}{A list echoing the validated call parameters for
#'     reproducibility.}
#' }
#'
#' @examples
#' ## -----------------------------------------------------------------------
#' ## Example 1 — numeric binary response (0/1), all predictors
#' ## -----------------------------------------------------------------------
#' set.seed(42)
#' n  <- 120
#' df <- data.frame(
#'   group    = sample(0:1, n, replace = TRUE),
#'   marker_A = rnorm(n, mean = ifelse(sample(0:1, n, replace = TRUE) == 1, 1, 0)),
#'   marker_B = rnorm(n),
#'   marker_C = runif(n)
#' )
#'
#' res <- run_auc(x = df, y = "group", plot = "none")
#' print(res$auc_table)
#'
#' ## -----------------------------------------------------------------------
#' ## Example 2 — categorical response, remove one group, combined plot
#' ## -----------------------------------------------------------------------
#' df2 <- data.frame(
#'   status   = sample(c("control", "case", "unknown"), n, replace = TRUE),
#'   score_1  = rnorm(n),
#'   score_2  = rnorm(n)
#' )
#'
#' res2 <- run_auc(
#'   x           = df2,
#'   y           = "status",
#'   remove      = "unknown",
#'   set_control = "control",
#'   plot        = "none",
#'   compute_ci  = FALSE
#' )
#' print(res2$auc_table)
#'
#' ## -----------------------------------------------------------------------
#' ## Example 3 — column names with spaces, select top-2 by AUC
#' ## -----------------------------------------------------------------------
#' df3 <- data.frame(
#'   `my group`  = sample(0:1, n, replace = TRUE),
#'   `feature 1` = rnorm(n, mean = 0.5),
#'   `feature 2` = rnorm(n),
#'   check.names = FALSE
#' )
#'
#' res3 <- run_auc(
#'   x    = df3,
#'   y    = "my group",
#'   plot = "none"
#' )
#' print(res3$auc_table)
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom pROC roc auc ci
#' @importFrom ggplot2 ggplot aes geom_line geom_abline annotate labs theme_minimal
#'   theme_bw theme_classic theme_gray theme_light theme_dark theme element_text
#'   element_rect element_blank scale_colour_manual guides guide_legend
#'   scale_x_continuous scale_y_continuous coord_equal .pt
#' @importFrom stats setNames
#' @importFrom utils head
#'
#' @export
run_auc <- function(
    x,
    y,
    remove             = NULL,
    metadata           = NULL,
    set_control        = NULL,
    include_cols       = NULL,
    include_rows       = NULL,
    direction          = "auto",
    compute_auc        = TRUE,
    compute_ci         = TRUE,
    ci_level           = 0.95,
    plot               = "none",
    theme              = "nature",
    ncol               = NULL,
    nrow               = NULL,
    plot_title         = NULL,
    plot_subtitle      = NULL,
    plot_legend        = "AUC (95% CI)",
    xlab               = NULL,
    ylab               = NULL,
    linewidth          = 1,
    global_font_size   = 20,
    title_font_size    = NULL,
    subtitle_font_size = NULL,
    xlab_font_size     = NULL,
    ylab_font_size     = NULL,
    inplot_font_size   = NULL,
    ...
) {

  # ---------------------------------------------------------------------------
  # 0.  Namespace helpers (avoid hard importing every symbol)
  # ---------------------------------------------------------------------------
  .roc <- pROC::roc
  .auc <- pROC::auc
  .ci  <- pROC::ci

  # ---------------------------------------------------------------------------
  # 1.  Validate x
  # ---------------------------------------------------------------------------
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop(
      "`x` must be a data.frame, tibble, or named matrix. ",
      "Received an object of class: ", paste(class(x), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop(
        "`x` is a matrix but has no column names. ",
        "Please set column names with `colnames(x) <- ...` before calling run_auc().",
        call. = FALSE
      )
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  if (nrow(x) == 0L) stop("`x` has zero rows.", call. = FALSE)
  if (ncol(x) == 0L) stop("`x` has zero columns.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2.  Validate y
  # ---------------------------------------------------------------------------
  if (missing(y) || is.null(y)) {
    stop("`y` is required. Provide the column name in `x` that contains the grouping variable.",
         call. = FALSE)
  }
  if (!is.character(y) || length(y) != 1L) {
    stop("`y` must be a single character string (a column name of `x`).", call. = FALSE)
  }
  if (!y %in% colnames(x)) {
    stop(
      "`y` = \"", y, "\" is not a column name in `x`. ",
      "Available columns: ", paste(colnames(x), collapse = ", "), ".",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 3.  Validate direction
  # ---------------------------------------------------------------------------
  direction <- match.arg(direction, choices = c("auto", ">", "<"))

  # ---------------------------------------------------------------------------
  # 4.  Validate scalar logical / numeric arguments
  # ---------------------------------------------------------------------------
  .check_bool <- function(val, name) {
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop("`", name, "` must be a single non-NA logical value (TRUE or FALSE).", call. = FALSE)
    }
  }
  .check_bool(compute_auc, "compute_auc")
  .check_bool(compute_ci,  "compute_ci")

  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      is.na(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single numeric value strictly between 0 and 1 (e.g., 0.95).",
         call. = FALSE)
  }

  .check_pos_num <- function(val, name, allow_null = TRUE) {
    if (allow_null && is.null(val)) return(invisible(NULL))
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop("`", name, "` must be a single positive numeric value.", call. = FALSE)
    }
  }
  .check_pos_num(linewidth,          "linewidth",          allow_null = FALSE)
  .check_pos_num(global_font_size,   "global_font_size",   allow_null = FALSE)
  .check_pos_num(title_font_size,    "title_font_size")
  .check_pos_num(subtitle_font_size, "subtitle_font_size")
  .check_pos_num(xlab_font_size,     "xlab_font_size")
  .check_pos_num(ylab_font_size,     "ylab_font_size")
  .check_pos_num(inplot_font_size,   "inplot_font_size")

  # ---------------------------------------------------------------------------
  # 5.  Resolve font sizes
  #    Proportions relative to global_font_size:
  #      title      -> 1.20 x
  #      subtitle   -> 0.90 x
  #      xlab/ylab  -> 1.00 x
  #      inplot     -> 0.85 x
  #    Hard defaults (used when global_font_size is also NULL):
  #      title=14, subtitle=11, xlab/ylab=12, inplot=11
  # ---------------------------------------------------------------------------
  gfs   <- global_font_size
  fs_title    <- if (!is.null(title_font_size))    title_font_size    else gfs * 1.20
  fs_subtitle <- if (!is.null(subtitle_font_size)) subtitle_font_size else gfs * 0.90
  fs_xlab     <- if (!is.null(xlab_font_size))     xlab_font_size     else gfs * 1.00
  fs_ylab     <- if (!is.null(ylab_font_size))     ylab_font_size     else gfs * 1.00
  fs_inplot   <- if (!is.null(inplot_font_size))   inplot_font_size   else gfs * 0.85

  # ---------------------------------------------------------------------------
  # 6.  Subset rows
  # ---------------------------------------------------------------------------
  rn <- rownames(x)
  has_named_rows <- !is.null(rn) && !identical(rn, as.character(seq_len(nrow(x))))

  if (!is.null(include_rows)) {
    if (is.character(include_rows)) {
      if (!has_named_rows) {
        stop(
          "`include_rows` was supplied as a character vector but `x` has no named rows. ",
          "Use integer positions instead.", call. = FALSE
        )
      }
      bad <- setdiff(include_rows, rn)
      if (length(bad) > 0L) {
        stop(
          "The following values in `include_rows` are not row names of `x`: ",
          paste(bad, collapse = ", "), ".", call. = FALSE
        )
      }
      x <- x[include_rows, , drop = FALSE]
    } else if (is.numeric(include_rows)) {
      include_rows <- as.integer(include_rows)
      bad <- include_rows[include_rows < 1L | include_rows > nrow(x)]
      if (length(bad) > 0L) {
        stop(
          "The following `include_rows` positions are out of range [1, ", nrow(x), "]: ",
          paste(bad, collapse = ", "), ".", call. = FALSE
        )
      }
      x <- x[include_rows, , drop = FALSE]
    } else {
      stop("`include_rows` must be a character or integer vector.", call. = FALSE)
    }
    if (nrow(x) == 0L) stop("`include_rows` selection resulted in zero rows.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 7.  Subset columns — keep y always, only subset predictors
  # ---------------------------------------------------------------------------
  all_cols <- colnames(x)

  if (!is.null(include_cols)) {
    if (is.character(include_cols)) {
      bad <- setdiff(include_cols, all_cols)
      if (length(bad) > 0L) {
        stop(
          "The following `include_cols` names are not columns in `x`: ",
          paste(bad, collapse = ", "), ".", call. = FALSE
        )
      }
      pred_cols <- setdiff(include_cols, y)
    } else if (is.numeric(include_cols)) {
      include_cols <- as.integer(include_cols)
      bad <- include_cols[include_cols < 1L | include_cols > length(all_cols)]
      if (length(bad) > 0L) {
        stop(
          "The following `include_cols` positions are out of range [1, ",
          length(all_cols), "]: ", paste(bad, collapse = ", "), ".", call. = FALSE
        )
      }
      pred_cols <- setdiff(all_cols[include_cols], y)
    } else {
      stop("`include_cols` must be a character or integer vector.", call. = FALSE)
    }
  } else {
    pred_cols <- setdiff(all_cols, y)
  }

  if (length(pred_cols) == 0L) {
    stop(
      "No predictor columns remain after excluding `y` (\"", y, "\") and applying `include_cols`. ",
      "Ensure at least one predictor column is selected.", call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 8.  Extract and process the response column
  # ---------------------------------------------------------------------------
  response_raw <- x[[y]]

  ## -- 8a. Apply `remove` ----------------------------------------------------
  if (!is.null(remove)) {
    if (!is.vector(remove) && !is.list(remove)) {
      stop("`remove` must be a vector or list of values to exclude from `y`.", call. = FALSE)
    }
    remove_vals <- unlist(remove)
    keep_idx    <- !(response_raw %in% remove_vals)
    if (sum(keep_idx) == 0L) {
      stop(
        "After applying `remove`, no rows remain. ",
        "Check that the values in `remove` are a proper subset of those in `y`.",
        call. = FALSE
      )
    }
    x            <- x[keep_idx, , drop = FALSE]
    response_raw <- x[[y]]
  }

  ## -- 8b. Determine unique levels / values -----------------------------------
  is_numeric_y <- is.numeric(response_raw)

  if (is_numeric_y) {
    uniq_vals <- sort(unique(response_raw[!is.na(response_raw)]))
    ## Detect 0/1 coding
    if (!all(uniq_vals %in% c(0, 1))) {
      stop(
        "When `y` is numeric it must contain only 0 (control) and 1 (case) after removal. ",
        "Found values: ", paste(uniq_vals, collapse = ", "), ". ",
        "Use `remove` to exclude unwanted levels, or re-code your variable.", call. = FALSE
      )
    }
    if (length(uniq_vals) < 2L) {
      stop(
        "After applying `remove`, `y` (numeric) contains only one distinct value (",
        uniq_vals, "). Two values (0 and 1) are required.", call. = FALSE
      )
    }
    response <- response_raw
    ctrl_label <- "0"
    case_label <- "1"

  } else {
    ## Categorical path
    response_char <- as.character(response_raw)
    uniq_levels   <- sort(unique(response_char[!is.na(response_char)]))

    ## Silent conversion if the two levels are literally "0" and "1"
    if (setequal(uniq_levels, c("0", "1"))) {
      warning(
        "`y` appears to encode 0/1 as a character/factor. Converting to numeric. ",
        "0 = control, 1 = case.", call. = FALSE
      )
      response   <- as.numeric(response_char)
      ctrl_label <- "0"
      case_label <- "1"
      is_numeric_y <- TRUE
    } else {
      if (length(uniq_levels) != 2L) {
        stop(
          "After applying `remove`, `y` has ", length(uniq_levels),
          " level(s): ", paste(uniq_levels, collapse = ", "), ". ",
          "Exactly 2 levels are required. Use `remove` to eliminate unwanted categories.",
          call. = FALSE
        )
      }
      ## Assign control / case
      if (!is.null(set_control)) {
        if (!is.character(set_control) || length(set_control) != 1L) {
          stop("`set_control` must be a single character string.", call. = FALSE)
        }
        if (!set_control %in% uniq_levels) {
          stop(
            "`set_control` = \"", set_control,
            "\" is not one of the two remaining levels in `y`: ",
            paste(uniq_levels, collapse = ", "), ".", call. = FALSE
          )
        }
        ctrl_label <- set_control
        case_label <- setdiff(uniq_levels, set_control)
      } else {
        ctrl_label <- uniq_levels[1L]
        case_label <- uniq_levels[2L]
        message(
          "No `set_control` specified. Treating \"", ctrl_label,
          "\" as control and \"", case_label, "\" as case."
        )
      }
      response <- ifelse(response_char == case_label, 1L, 0L)
    }
  }

  ## -- 8c. Final NA check on response ----------------------------------------
  na_resp <- sum(is.na(response))
  if (na_resp > 0L) {
    warning(
      na_resp, " NA value(s) found in `y` after processing. ",
      "Rows with NA in `y` will be excluded from each ROC analysis.", call. = FALSE
    )
  }

  n_cases    <- sum(response == 1L, na.rm = TRUE)
  n_controls <- sum(response == 0L, na.rm = TRUE)

  if (n_cases == 0L)    stop("No cases (value = 1) remain in `y`. Check `remove` and `y`.", call. = FALSE)
  if (n_controls == 0L) stop("No controls (value = 0) remain in `y`. Check `remove` and `y`.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 9.  Validate `plot` argument
  # ---------------------------------------------------------------------------
  .parse_plot_arg <- function(plot_arg) {
    ## Returns a list: list(mode = "none"|"all"|"separate",
    ##                      filter_type = "none"|"top_n"|"auc_thresh",
    ##                      filter_val  = NA_real_)
    valid_modes <- c("none", "all", "separate")

    if (is.character(plot_arg) && length(plot_arg) == 1L) {
      plot_arg <- match.arg(plot_arg, choices = valid_modes)
      return(list(mode = plot_arg, filter_type = "none", filter_val = NA_real_))
    }

    if (is.numeric(plot_arg) && length(plot_arg) == 1L) {
      if (plot_arg >= 1) {
        n <- as.integer(plot_arg)
        if (n < 1L) stop("`plot` integer value must be >= 1.", call. = FALSE)
        return(list(mode = "all", filter_type = "top_n", filter_val = n))
      } else if (plot_arg >= 0 && plot_arg < 1) {
        return(list(mode = "all", filter_type = "auc_thresh", filter_val = plot_arg))
      } else {
        stop("`plot` numeric must be >= 0.", call. = FALSE)
      }
    }

    ## vector case: c("all"|"separate", numeric)
    if ((is.character(plot_arg) || is.list(plot_arg)) || length(plot_arg) == 2L) {
      pv <- unlist(plot_arg)
      str_part <- pv[!suppressWarnings(is.na(as.numeric(pv)))]
      num_part <- suppressWarnings(as.numeric(pv[!pv %in% valid_modes]))
      str_part <- pv[pv %in% valid_modes]

      if (length(str_part) != 1L || length(num_part) != 1L) {
        stop(
          "`plot` as a vector must be of the form c(\"all\", n) or c(\"separate\", n) ",
          "where n is an integer >= 1 or a numeric in [0, 1).", call. = FALSE
        )
      }
      mode <- match.arg(str_part, choices = c("all", "separate"))
      nval <- num_part
      if (nval >= 1) {
        return(list(mode = mode, filter_type = "top_n",     filter_val = as.integer(nval)))
      } else if (nval >= 0 && nval < 1) {
        return(list(mode = mode, filter_type = "auc_thresh", filter_val = nval))
      } else {
        stop("Numeric element in `plot` must be >= 0.", call. = FALSE)
      }
    }

    stop(
      "Unrecognised `plot` argument. See ?run_auc for valid options.",
      call. = FALSE
    )
  }

  plot_spec <- tryCatch(
    .parse_plot_arg(plot),
    error = function(e) stop(conditionMessage(e), call. = FALSE)
  )

  # ---------------------------------------------------------------------------
  # 10.  Validate theme
  # ---------------------------------------------------------------------------
  valid_themes <- c("nature", "minimal", "bw", "classic", "plos", "gray", "light", "dark")
  if (!is.character(theme) || length(theme) != 1L || !theme %in% valid_themes) {
    stop(
      "`theme` must be one of: ", paste(valid_themes, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 11.  Validate ncol / nrow
  # ---------------------------------------------------------------------------
  .check_int_null <- function(val, name) {
    if (is.null(val)) return(invisible(NULL))
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val < 1L) {
      stop("`", name, "` must be a single positive integer or NULL.", call. = FALSE)
    }
  }
  .check_int_null(ncol, "ncol")
  .check_int_null(nrow, "nrow")

  # ---------------------------------------------------------------------------
  # 12.  Core ROC computation loop
  # ---------------------------------------------------------------------------
  extra_args <- list(...)

  ## Helper: build a safe pROC::roc call, suppressing its internal messages.
  ## We filter extra_args against the formals of the exported pROC::roc
  ## generic — NOT pROC::roc.default, which is not exported and cannot be
  ## introspected reliably across pROC versions.
  .roc_formal_names <- tryCatch(
    names(formals(pROC::roc)),
    error = function(e) character(0L)
  )

  .safe_roc <- function(resp, pred, dir) {
    roc_args <- c(
      list(response  = resp,
           predictor = pred,
           direction = dir,
           levels    = c(0, 1),
           quiet     = TRUE),
      extra_args[names(extra_args) %in% .roc_formal_names]
    )
    suppressMessages(do.call(pROC::roc, roc_args))
  }

  roc_list   <- vector("list", length(pred_cols))
  names(roc_list) <- pred_cols

  auc_vec    <- numeric(length(pred_cols))
  ci_lower   <- rep(NA_real_, length(pred_cols))
  ci_upper   <- rep(NA_real_, length(pred_cols))
  dir_used   <- character(length(pred_cols))

  for (i in seq_along(pred_cols)) {
    col_nm  <- pred_cols[i]
    pred_i  <- x[[col_nm]]

    if (!is.numeric(pred_i)) {
      warning(
        "Column \"", col_nm, "\" is not numeric and will be skipped.",
        call. = FALSE
      )
      auc_vec[i]  <- NA_real_
      dir_used[i] <- NA_character_
      next
    }

    valid_idx <- !is.na(response) & !is.na(pred_i)
    if (sum(valid_idx) < 2L) {
      warning(
        "Column \"", col_nm,
        "\" has fewer than 2 complete observations and will be skipped.",
        call. = FALSE
      )
      auc_vec[i]  <- NA_real_
      dir_used[i] <- NA_character_
      next
    }

    roc_i <- tryCatch(
      .safe_roc(response[valid_idx], pred_i[valid_idx], direction),
      error = function(e) {
        warning("pROC::roc() failed for \"", col_nm, "\": ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )

    if (is.null(roc_i)) {
      auc_vec[i]  <- NA_real_
      dir_used[i] <- NA_character_
      next
    }

    roc_list[[i]] <- roc_i
    dir_used[i]   <- as.character(roc_i$direction)

    if (compute_auc) {
      auc_i       <- as.numeric(.auc(roc_i))
      auc_vec[i]  <- auc_i
    }

    if (compute_ci && compute_auc) {
      ci_args <- c(
        list(roc_i, conf.level = ci_level),
        extra_args[names(extra_args) %in% names(formals(pROC::ci.auc))]
      )
      ci_i <- tryCatch(
        suppressMessages(do.call(.ci, ci_args)),
        error = function(e) {
          warning("pROC::ci() failed for \"", col_nm, "\": ", conditionMessage(e),
                  call. = FALSE)
          NULL
        }
      )
      if (!is.null(ci_i)) {
        ci_lower[i] <- as.numeric(ci_i)[1L]
        ci_upper[i] <- as.numeric(ci_i)[3L]
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 13.  Build AUC summary table
  # ---------------------------------------------------------------------------
  auc_table <- data.frame(
    predictor  = pred_cols,
    auc        = auc_vec,
    ci_lower   = ci_lower,
    ci_upper   = ci_upper,
    ci_level   = ci_level,
    direction  = dir_used,
    n_cases    = n_cases,
    n_controls = n_controls,
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # 14.  AUC label formatter (smart decimal / sci-notation)
  # ---------------------------------------------------------------------------
  .fmt_auc_val <- function(val) {
    if (is.na(val)) return(NA_character_)
    for (dp in 2L:5L) {
      fmt_val <- formatC(val, digits = dp, format = "f")
      ## Does it round to exactly 1.000...?
      if (as.numeric(fmt_val) < 1) return(fmt_val)
    }
    ## Fall back to scientific notation
    formatC(val, digits = 2L, format = "e")
  }

  .fmt_auc_label <- function(auc_v, lo_v, hi_v, ci_lv) {
    pct <- round(ci_lv * 100)
    paste0(
      .fmt_auc_val(auc_v),
      " (", pct, "% CI: ",
      .fmt_auc_val(lo_v), "\u2013", .fmt_auc_val(hi_v), ")"
    )
  }

  # ---------------------------------------------------------------------------
  # 15.  Colour-blind-friendly palette (Okabe-Ito, 8 colours)
  # ---------------------------------------------------------------------------
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )

  .get_palette <- function(n) {
    if (n <= length(okabe_ito)) return(okabe_ito[seq_len(n)])
    ## Recycle with a warning
    warning(
      "More predictors than palette colours (", length(okabe_ito),
      "). Colours will be recycled.", call. = FALSE
    )
    rep_len(okabe_ito, n)
  }

  # ---------------------------------------------------------------------------
  # 16.  ggplot2 theme builder
  # ---------------------------------------------------------------------------
  .build_theme <- function(theme_name) {
    base <- switch(
      theme_name,
      nature  = , plos = ggplot2::theme_classic(base_size = gfs),
      minimal = ggplot2::theme_minimal(base_size = gfs),
      bw      = ggplot2::theme_bw(base_size = gfs),
      classic = ggplot2::theme_classic(base_size = gfs),
      gray    = ggplot2::theme_gray(base_size = gfs),
      light   = ggplot2::theme_light(base_size = gfs),
      dark    = ggplot2::theme_dark(base_size = gfs)
    )
    base +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(size = fs_title,    face = "bold",
                                              hjust = 0),
        plot.subtitle = ggplot2::element_text(size = fs_subtitle, hjust = 0,
                                              colour = "grey40"),
        axis.title.x  = ggplot2::element_text(size = fs_xlab),
        axis.title.y  = ggplot2::element_text(size = fs_ylab),
        axis.text.x   = ggplot2::element_text(size = fs_xlab),
        axis.text.y   = ggplot2::element_text(size = fs_ylab),
        legend.text   = ggplot2::element_text(size = fs_inplot),
        legend.title  = ggplot2::element_text(size = fs_inplot, face = "bold")
      )
  }

  # ---------------------------------------------------------------------------
  # 17.  Determine which predictors to plot
  # ---------------------------------------------------------------------------
  plots_out <- NULL

  if (plot_spec$mode != "none") {

    ## Complete-case AUC table for plotting decisions
    valid_auc <- auc_table[!is.na(auc_table$auc), ]

    ## Apply filters
    plot_pred_cols <- if (plot_spec$filter_type == "top_n") {
      ord  <- order(valid_auc$auc, decreasing = TRUE)
      top  <- head(ord, plot_spec$filter_val)
      valid_auc$predictor[top]
    } else if (plot_spec$filter_type == "auc_thresh") {
      valid_auc$predictor[valid_auc$auc >= plot_spec$filter_val]
    } else {
      valid_auc$predictor
    }

    if (length(plot_pred_cols) == 0L) {
      warning("No predictors satisfy the plot filter criteria. No plots will be produced.",
              call. = FALSE)
      plot_spec$mode <- "none"
    }
  }

  # ---------------------------------------------------------------------------
  # 18.  Build ROC data frames for ggplot
  # ---------------------------------------------------------------------------
  if (plot_spec$mode != "none") {
    x_lbl <- if (!is.null(xlab)) xlab else "1 - Specificity (False Positive Rate)"
    y_lbl <- if (!is.null(ylab)) ylab else "Sensitivity (True Positive Rate)"

    ## Gather per-predictor curve data (unsorted here; sorting happens inside
    ## .make_roc_plot so both "all" and "separate" modes benefit).
    curve_data_list <- lapply(plot_pred_cols, function(col_nm) {
      roc_obj <- roc_list[[col_nm]]
      if (is.null(roc_obj)) return(NULL)
      data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities,
        stringsAsFactors = FALSE
      )
    })
    names(curve_data_list) <- plot_pred_cols
    curve_data_list        <- Filter(Negate(is.null), curve_data_list)

    ## Validate plot_title length for separate mode
    if (plot_spec$mode == "separate" && !is.null(plot_title)) {
      if (length(plot_title) != length(plot_pred_cols)) {
        stop(
          "When `plot = \"separate\"`, `plot_title` must have the same length as the number ",
          "of predictors being plotted (", length(plot_pred_cols), "). ",
          "Supplied: ", length(plot_title), " title(s).", call. = FALSE
        )
      }
    }

    pal <- .get_palette(length(curve_data_list))

    ## Helper: build AUC label for legend
    .auc_lbl <- function(col_nm) {
      row <- auc_table[auc_table$predictor == col_nm, ]
      if (compute_ci && !is.na(row$ci_lower)) {
        .fmt_auc_label(row$auc, row$ci_lower, row$ci_upper, ci_level)
      } else if (compute_auc) {
        .fmt_auc_val(row$auc)
      } else {
        col_nm
      }
    }

    ## Helper: single ROC ggplot (combined or single-predictor)
    .make_roc_plot <- function(df_list, colours, title_str, subtitle_str,
                               legend_title_str, in_legend = TRUE) {
      ## Build name -> AUC label mapping (predictor name is always unique;
      ## two predictors can share the same rounded label string, which would
      ## cause duplicated factor levels if used directly as the colour key).
      #lbl_map <- sapply(names(df_list), .auc_lbl)   # named: predictor -> label
      lbl_map <- sapply(names(df_list), function(nm) {
        paste0(nm, .auc_lbl(nm))
      })

      ## Stack curve data; sort each curve by FPR then TPR so geom_line
      ## draws a smooth, non-jagged path (mirrors the perform_AUROC approach).
      all_curves <- do.call(rbind, lapply(names(df_list), function(nm) {
        df <- df_list[[nm]]
        df <- df[order(df$fpr, df$tpr), ]
        df$predictor <- nm          # keep original predictor name as colour key
        df
      }))

      ## Predictor column as an ordered factor so legend order matches AUC rank
      all_curves$predictor <- factor(all_curves$predictor,
                                     levels = names(df_list))

      p <- ggplot2::ggplot(
        all_curves,
        ggplot2::aes(x = fpr, y = tpr,
                     colour = predictor,
                     group  = predictor)
      ) +
        ## Reference line
        ggplot2::geom_abline(
          slope = 1, intercept = 0,
          linetype = "dashed", colour = "grey60", linewidth = linewidth * 0.7
        ) +
        ## ROC curves
        ggplot2::geom_line(linewidth = linewidth) +
        ## Map predictor names -> colours, but show AUC labels in the legend
        ggplot2::scale_colour_manual(
          name   = legend_title_str,
          values = setNames(colours, names(df_list)),
          labels = lbl_map           # display AUC label, key on predictor name
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title    = title_str,
          subtitle = subtitle_str,
          x        = x_lbl,
          y        = y_lbl
        ) +
        .build_theme(theme)

      if (in_legend) {
        ## Bottom-right legend inside plot
        p <- p + ggplot2::theme(
          legend.position        = c(0.98, 0.02),
          legend.justification   = c(1, 0),
          legend.background      = ggplot2::element_rect(fill = "white",
                                                          colour = "grey80",
                                                          linewidth = 0.4),
          legend.key             = ggplot2::element_blank(),
          legend.margin          = ggplot2::margin(4, 6, 4, 6)
        )
        if (!is.null(ncol) || !is.null(nrow)) {
          guide_args <- list(ncol = ncol, nrow = nrow)
          guide_args <- Filter(Negate(is.null), guide_args)
          p <- p + ggplot2::guides(
            colour = do.call(ggplot2::guide_legend, guide_args)
          )
        }
      } else {
        p <- p + ggplot2::theme(legend.position = "none")
      }

      p
    }

    # -------------------------------------------------------------------------
    # 19.  Produce plots according to mode
    # -------------------------------------------------------------------------
    plots_out <- list()
    legend_title_str <- if (!is.null(plot_legend)) plot_legend else "AUC (95% CI)"

    if (plot_spec$mode == "all") {
      title_str <- if (!is.null(plot_title)) plot_title[1L] else "ROC Curve"
      p_all <- .make_roc_plot(
        df_list          = curve_data_list,
        colours          = pal,
        title_str        = title_str,
        subtitle_str     = plot_subtitle,
        legend_title_str = legend_title_str,
        in_legend        = TRUE
      )
      plots_out[["combined"]] <- p_all

    } else if (plot_spec$mode == "separate") {
      for (j in seq_along(curve_data_list)) {
        col_nm   <- names(curve_data_list)[j]
        t_str    <- if (!is.null(plot_title)) {
          plot_title[j]
        } else {
          paste0("ROC Curve for ", col_nm)
        }
        p_sep <- .make_roc_plot(
          df_list          = curve_data_list[j],
          colours          = pal[j],
          title_str        = t_str,
          subtitle_str     = plot_subtitle,
          legend_title_str = legend_title_str,
          in_legend        = FALSE
        )

        ## Build the in-plot AUC annotation label for this predictor
        row_i     <- auc_table[auc_table$predictor == col_nm, ]
        auc_label <- if (compute_ci && !is.na(row_i$ci_lower)) {
          .fmt_auc_label(row_i$auc, row_i$ci_lower, row_i$ci_upper, ci_level)
        } else if (compute_auc && !is.na(row_i$auc)) {
          .fmt_auc_val(row_i$auc)
        } else {
          NULL
        }

        if (!is.null(auc_label)) {
          ## legend_title_str (e.g. "AUC (95% CI)") as a bold header,
          ## then the value on the next line — placed at bottom-right
          ## using annotation_custom + a ggplot2 text grob so that the
          ## position is in data coordinates (0–1 on both axes).
          ann_label <- paste0(legend_title_str, "\n", auc_label)
          p_sep <- p_sep +
            ggplot2::annotate(
              geom     = "label",
              x        = 0.98,
              y        = 0.04,
              label    = ann_label,
              hjust    = 1,
              vjust    = 0,
              size     = fs_inplot / ggplot2::.pt,   # convert pt -> mm for ggplot2
              fill     = "white",
              colour   = "grey20",
              label.size = 0.3                        # border thickness
            )
        }

        plots_out[[col_nm]] <- p_sep
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 20.  Assemble return value
  # ---------------------------------------------------------------------------
  result <- structure(
    list(
      roc_list  = roc_list,
      auc_table = auc_table,
      plots     = plots_out,
      params    = list(
        y           = y,
        set_control = if (!is_numeric_y) ctrl_label else NULL,
        remove      = remove,
        include_cols = pred_cols,
        direction   = direction,
        compute_auc = compute_auc,
        compute_ci  = compute_ci,
        ci_level    = ci_level,
        n_cases     = n_cases,
        n_controls  = n_controls,
        call        = match.call()
      )
    ),
    class = "run_auc"
  )

  result
}


# =============================================================================
#  S3 print method
# =============================================================================

#' @title Print Method for \code{run_auc} Objects
#'
#' @description Compact summary of a \code{run_auc} result.
#'
#' @param x An object of class \code{"run_auc"}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @author John Lennon L. Calorio
#' @export
print.run_auc <- function(x, ...) {
  cat("── run_auc results ─────────────────────────────────────────\n")
  cat(sprintf(
    "  Response : %s  |  Cases: %d  |  Controls: %d\n",
    x$params$y, x$params$n_cases, x$params$n_controls
  ))
  if (!is.null(x$params$set_control)) {
    cat(sprintf("  Control  : \"%s\"\n", x$params$set_control))
  }
  cat(sprintf(
    "  CI level : %.0f%%   |  Direction: %s\n",
    x$params$ci_level * 100, x$params$direction
  ))
  cat("\n  AUC Table:\n")
  tbl <- x$auc_table
  tbl$auc      <- round(tbl$auc,      4)
  tbl$ci_lower <- round(tbl$ci_lower, 4)
  tbl$ci_upper <- round(tbl$ci_upper, 4)
  print(tbl[, c("predictor", "auc", "ci_lower", "ci_upper", "direction")],
        row.names = FALSE)
  if (!is.null(x$plots) && length(x$plots) > 0L) {
    cat(sprintf("\n  Plots stored: %s\n", paste(names(x$plots), collapse = ", ")))
  }
  invisible(x)
}