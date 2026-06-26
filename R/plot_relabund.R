#' Plot Relative Abundance as Stacked Bar Chart
#'
#' @title Plot Relative Abundance as Stacked Bar Chart
#'
#' @description
#' Generates a stacked bar chart of relative abundances from a feature matrix,
#' optionally restricting the display to the top or bottom \code{limit}
#' features by \emph{total} abundance and collapsing the remainder into an
#' \code{"Others"} category. Samples may be shown individually or aggregated
#' by a grouping variable by summing raw counts per group before normalising.
#' Both horizontal and vertical orientations are supported, along with full
#' control over colours, fonts, legend layout, and axis labels. Column names
#' containing spaces or special characters are fully supported. When plotting
#' individual samples, bracket annotations can be drawn beneath (or beside)
#' the bars to visually group samples by a metadata variable.
#'
#' @param x A \code{data.frame}, \code{tibble}, or \code{matrix} with named
#'   columns representing features (e.g., taxa, metabolites). Rows are samples;
#'   columns are features. Column names may contain spaces and special
#'   characters. \strong{When \code{group_aggregate = FALSE} and
#'   \code{ignore_normalization = FALSE}}, \code{x} is expected to already
#'   contain row-wise relative abundances (each row summing to 1). A warning
#'   is issued for any row whose sum is not 0 or 1, indicating that row-wise
#'   relative abundance normalisation may not have been applied.
#' @param metadata A \code{data.frame} containing sample-level metadata. Must
#'   have the same number of rows as \code{x}; row order is assumed to
#'   correspond to \code{x}.
#' @param limit Integer or \code{NULL}. Number of features to display, selected
#'   by \emph{total} abundance across all samples. When \code{sort = "high"}
#'   (or \code{NULL}/\code{"alphanumeric"} with a limit), the top \code{limit}
#'   highest-total features are shown and the rest collapsed to \code{"Others"}.
#'   When \code{sort = "low"}, the bottom \code{limit} lowest-total features
#'   are shown instead. If \code{NULL}, all features are plotted. Default is
#'   \code{10}.
#' @param sort Character or \code{NULL}. Controls both which features are
#'   selected (when \code{limit} is set) and the order features appear in the
#'   legend/stack. \code{NULL} (default) preserves the original column order of
#'   \code{x}. Options: \code{"alphanumeric"} (lexicographic),
#'   \code{"high"} (highest total abundance first — selects top \code{limit}),
#'   \code{"low"} (lowest total abundance first — selects bottom \code{limit}).
#'   \strong{Note:} \code{"high"} and \code{"low"} rank by column totals and
#'   are meaningful only when \code{x} contains raw counts or intensities
#'   (row sums not already 0 or 1). If normalised data is detected, a warning
#'   is issued and \code{sort} is ignored. When using ranked sorting, supply
#'   raw count data and set \code{normalize_to_relabund = TRUE}.
#' @param others_position Character or \code{NULL}. Where to place the
#'   \code{"Others"} category in the legend/stack order. \code{NULL} (default)
#'   treats \code{"Others"} according to \code{sort}. \code{"first"} forces it
#'   to the top; \code{"last"} forces it to the bottom.
#' @param show_x_as Character or \code{NULL}. Column name in \code{metadata}
#'   to use as the x-axis variable (sample IDs or group labels). If \code{NULL},
#'   the function searches for \code{"UID"}, \code{"ID"}, or \code{"SampleID"}
#'   (case-insensitive) in \code{metadata}. An error is raised if none are
#'   found.
#' @param group_aggregate Logical. If \code{TRUE}, \code{show_x_as} is treated
#'   as a grouping variable and raw feature values in \code{x} are summed per
#'   group before normalising to relative abundances. \strong{Assumes \code{x}
#'   contains raw (unnormalised) data.} If \code{FALSE} (default), each row in
#'   \code{x} is treated as an individual sample. Cannot be \code{TRUE}
#'   simultaneously with \code{group_brackets}.
#' @param normalize_to_relabund Logical. Relevant when \code{group_aggregate =
#'   TRUE} or \code{sort} is \code{"high"}/\code{"low"}. If \code{TRUE}
#'   (default), row-wise relative abundance normalisation is applied to \code{x}
#'   via \code{run_normalize(method = "row_rel_abundance")} before aggregation
#'   or ranking. Set to \code{FALSE} if \code{x} is already appropriately
#'   normalised. Ignored when \code{ignore_normalization = TRUE}.
#' @param ignore_normalization Logical. If \code{TRUE}, all normalisation checks
#'   and pre-normalisation steps are bypassed and \code{x} is plotted as
#'   supplied. Default is \code{FALSE}. \strong{Note:} Setting this to
#'   \code{TRUE} means the plotted values are not relative abundances; use this
#'   when plotting already-transformed data such as CLR-transformed values where
#'   the data originated from normalised data and should be visualised as-is.
#' @param order Character vector, \code{"alphanumeric"}, or \code{NULL}.
#'   Left-to-right order of samples/groups on the x-axis. A character vector
#'   gives an explicit sequence. \code{"alphanumeric"} sorts lexicographically.
#'   \code{NULL} (default) preserves natural order of appearance in
#'   \code{metadata}.
#' @param group_brackets A formula or \code{NULL}. Draws bracket annotations
#'   beneath (or beside, when \code{abundance_in_y = FALSE}) the bars to
#'   visually group individual samples by a metadata variable. Specify as a
#'   one-sided formula of the form \code{GroupCol ~ c("Level1", "Level2")},
#'   where \code{GroupCol} is a column in \code{metadata} and the right-hand
#'   side gives the display order of the groups. \strong{Cannot be used when
#'   \code{group_aggregate = TRUE}.} Default is \code{NULL}.
#' @param angle Numeric or \code{NULL}. Rotation angle (0--360 degrees) for
#'   sample/group tick labels. \code{NULL} (default): \code{0} when
#'   \code{group_aggregate = TRUE}, \code{45} otherwise.
#' @param abundance_in_y Logical. If \code{TRUE} (default), relative abundance
#'   is on the y-axis and samples/groups are on the x-axis. If \code{FALSE},
#'   axes are flipped: relative abundance on the x-axis, samples/groups on the
#'   y-axis. The legend order is also reversed when \code{FALSE} so the
#'   top-of-stack feature appears first (top) in the legend.
#' @param legend_position Character or \code{NULL}. One of \code{"top"},
#'   \code{"bottom"}, \code{"left"}, \code{"right"}. Auto-detected when
#'   \code{NULL}: \code{"top"} for \code{abundance_in_y = TRUE}; \code{"right"}
#'   for \code{abundance_in_y = FALSE}.
#' @param legend_nrow Integer or \code{NULL}. Rows in legend key grid.
#' @param legend_ncol Integer or \code{NULL}. Columns in legend key grid.
#' @param italicize_legend_items Logical. Render legend labels in italic.
#'   Requires \code{ggtext}; falls back to \code{element_text(face = "italic")}
#'   with a warning. Default \code{FALSE}.
#' @param exclude_unknown Logical. If \code{TRUE}, any feature whose name is
#'   exactly \code{"unknown"} (case-insensitive) is omitted from the plot
#'   entirely. Default is \code{FALSE}.
#' @param exclude_others Logical. If \code{TRUE}, the \code{"Others"} category
#'   is omitted from the plot entirely. Default \code{FALSE}.
#' @param unknown_before_others Logical. If \code{TRUE} (default), any feature
#'   whose name is exactly \code{"unknown"} (case-insensitive, e.g.,
#'   \code{"Unknown"}, \code{"UNKNOWN"}, \code{"unknown"}) is placed immediately
#'   before \code{"Others"} in the legend/stack, or last when
#'   \code{exclude_others = TRUE}. Partial matches are not considered; only
#'   column names that are \emph{entirely} the word "unknown" are affected.
#' @param theme Character. ggplot2 theme. One of \code{"nature"},
#'   \code{"minimal"}, \code{"classic"}, \code{"bw"}, \code{"light"},
#'   \code{"dark"}. Default \code{"nature"}.
#' @param colors Character vector or \code{NULL}. Custom fill colours. Length
#'   must equal the number of plotted fill levels. Named vectors matched by
#'   name. \code{NULL} uses the Okabe-Ito palette (extended with HCL colours).
#' @param plot_title Character. Main title. Supports \code{ggtext} markdown:
#'   \code{**bold**}, \code{*italic*}. Default \code{""}.
#' @param plot_subtitle Character. Subtitle. Supports same markdown. Default
#'   \code{""}.
#' @param legend_title Character. Legend title. Default \code{""}.
#' @param xlab Character or \code{NULL}. X-axis label. Defaults to
#'   \code{show_x_as}.
#' @param ylab Character. Y-axis label. Default \code{"Relative abundance"}.
#' @param global_font_size Numeric or \code{NULL}. Scales all text elements
#'   proportionally relative to an internal 11 pt reference. Default \code{20}.
#' @param font_family Character. Font family for all text. Default \code{"sans"}.
#' @param title_size Numeric or \code{NULL}. Title font size (pt).
#' @param subtitle_size Numeric or \code{NULL}. Subtitle font size (pt).
#' @param legend_size Numeric or \code{NULL}. Legend title font size (pt).
#' @param xlab_size Numeric or \code{NULL}. X-axis label font size (pt).
#' @param ylab_size Numeric or \code{NULL}. Y-axis label font size (pt).
#' @param ... Additional arguments passed to \code{ggplot2::geom_bar()}.
#'
#' @details
#' \strong{Special characters in column names.} Column names containing spaces,
#' slashes, parentheses, or other special characters are fully supported. The
#' function avoids \code{stats::reshape()} (which can mangle such names) and
#' instead builds the long-format data frame manually.
#'
#' \strong{Feature selection with \code{limit}.} Features are ranked by their
#' \emph{column total} across all samples (or aggregated groups). When
#' \code{sort = "high"} or \code{NULL} or \code{"alphanumeric"}, the top
#' \code{limit} features by total are retained; all others are collapsed into
#' \code{"Others"}. When \code{sort = "low"}, the bottom \code{limit} features
#' by total are retained instead.
#'
#' \strong{Normalisation flow.} When \code{ignore_normalization = FALSE} and
#' \code{normalize_to_relabund = TRUE}, row-relative abundance is applied via
#' \code{run_normalize(method = "row_rel_abundance")} before any aggregation or
#' ranking. Row-wise proportions are then computed so each bar sums to 100\%.
#' When \code{ignore_normalization = TRUE}, the data are plotted as supplied
#' and the y-axis label should be interpreted accordingly.
#'
#' \strong{Group aggregation.} When \code{group_aggregate = TRUE}, all rows
#' belonging to the same level of \code{show_x_as} are summed (raw counts),
#' and the result is then normalised to row-wise relative abundances. This is
#' incompatible with \code{group_brackets}.
#'
#' \strong{Group brackets.} When \code{group_brackets} is supplied (and
#' \code{group_aggregate = FALSE}), bracket annotations are drawn below each
#' set of sample bars that belong to the same group, with the group label
#' centred beneath the bracket. The formula right-hand side controls both
#' which groups appear and their left-to-right order. Any group present in
#' \code{metadata} but absent from the right-hand side is appended at the end.
#'
#' \strong{Legend reversal when flipped.} When \code{abundance_in_y = FALSE},
#' \code{coord_flip()} causes the stack to render left-to-right. The legend
#' fill levels are reversed so the first (top-of-stack) feature also appears
#' first (top) in the legend, maintaining visual correspondence.
#'
#' @return A \code{ggplot} object, also drawn to the active graphics device.
#'
#' @examples
#' \dontrun{
#' library(pondeR)
#' set.seed(42)
#' n    <- 20
#' taxa <- c("Firmicutes", "Bacteroidetes", "Proteobacteria (class A)",
#'           "Actino/bacteroides", "Verrucomicrobia", "Fusobacteria sp.",
#'           "Tenericutes", "Spirochaetes", "Unknown", "Cyanobacteria",
#'           "Chloroflexi", "Acidobacteria")
#'
#' counts <- as.data.frame(
#'   matrix(sample(0L:5000L, n * length(taxa), replace = TRUE),
#'          nrow = n, dimnames = list(NULL, taxa))
#' )
#' meta <- data.frame(
#'   SampleID = paste0("S", seq_len(n)),
#'   Group    = rep(c("Control", "Treatment"), each = n / 2L)
#' )
#'
#' # --- Pipeline: col_rel_abundance -> CLR -> plot as-is ---
#' norm_out <- run_normalize(counts, meta, method = "col_rel_abundance",
#'                           verbose = FALSE)
#' clr_out  <- run_transform(norm_out$data, method = "clr")
#'
#' plot_relabund(
#'   x                    = clr_out$data,
#'   metadata             = meta,
#'   limit                = 10L,
#'   ignore_normalization = TRUE,
#'   plot_title           = "Gut Microbiome (CLR-transformed)"
#' )
#'
#' # --- Group-aggregated from raw counts (group_aggregate = TRUE) ---
#' plot_relabund(
#'   x               = counts,
#'   metadata        = meta,
#'   show_x_as       = "Group",
#'   group_aggregate = TRUE,
#'   order           = c("Treatment", "Control"),
#'   limit           = 8L,
#'   sort            = "high",
#'   plot_title      = "**Relative Abundance** by Group",
#'   legend_title    = "Phylum"
#' )
#'
#' # --- Individual samples with group brackets (group_aggregate = FALSE) ---
#' plot_relabund(
#'   x                = counts,
#'   metadata         = meta,
#'   show_x_as        = "SampleID",
#'   group_aggregate  = FALSE,
#'   group_brackets   = Group ~ c("Control", "Treatment"),
#'   limit            = 8L,
#'   sort             = "high",
#'   plot_title       = "Relative Abundance by Sample",
#'   legend_title     = "Phylum"
#' )
#'
#' # --- Horizontal bars, individual samples, group brackets ---
#' plot_relabund(
#'   x                     = counts,
#'   metadata              = meta,
#'   show_x_as             = "SampleID",
#'   group_aggregate       = FALSE,
#'   group_brackets        = Group ~ c("Control", "Treatment"),
#'   limit                 = 6L,
#'   exclude_others        = FALSE,
#'   unknown_before_others = TRUE,
#'   abundance_in_y        = FALSE
#' )
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual scale_x_discrete
#' @importFrom ggplot2 scale_y_continuous coord_flip labs theme_bw theme_minimal theme_classic
#' @importFrom ggplot2 theme_light theme_dark theme element_text element_blank element_line
#' @importFrom ggplot2 element_rect guide_legend expansion waiver unit margin
#' @importFrom stats median setNames
#'
#' @export
plot_relabund <- function(
    x,
    metadata,
    limit                  = 10L,
    sort                   = NULL,
    others_position        = NULL,
    show_x_as              = NULL,
    group_aggregate        = FALSE,
    normalize_to_relabund  = TRUE,
    ignore_normalization   = FALSE,
    order                  = NULL,
    group_brackets         = NULL,
    angle                  = NULL,
    abundance_in_y         = TRUE,
    legend_position        = NULL,
    legend_nrow            = NULL,
    legend_ncol            = NULL,
    italicize_legend_items = FALSE,
    exclude_unknown        = FALSE,
    exclude_others         = FALSE,
    unknown_before_others  = TRUE,
    theme                  = "nature",
    colors                 = NULL,
    plot_title             = "",
    plot_subtitle          = "",
    legend_title           = "",
    xlab                   = NULL,
    ylab                   = "Relative abundance",
    global_font_size       = 20,
    font_family            = "sans",
    title_size             = NULL,
    subtitle_size          = NULL,
    legend_size            = NULL,
    xlab_size              = NULL,
    ylab_size              = NULL,
    ...
) {

  # ---------------------------------------------------------------------------
  # 0. Package availability
  # ---------------------------------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Required package 'ggplot2' is not installed. ",
      "Install it with: install.packages('ggplot2')"
    )
  }

  # ---------------------------------------------------------------------------
  # 1. Helpers
  # ---------------------------------------------------------------------------
  .is_tabular <- function(z) is.data.frame(z) || is.matrix(z)
  .check_pos_num <- function(val, name, allow_null = FALSE) {
    if (allow_null && is.null(val)) return(invisible(NULL))
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop("'", name, "' must be a single positive number.")
  }
  .check_pos_int <- function(val, name, allow_null = FALSE) {
    if (allow_null && is.null(val)) return(invisible(NULL))
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val < 1L || val != floor(val))
      stop("'", name, "' must be a single positive integer or NULL.")
  }

  # ---------------------------------------------------------------------------
  # 2. Validate x
  # ---------------------------------------------------------------------------
  if (!.is_tabular(x)) stop("'x' must be a data.frame, tibble, or matrix.")
  if (nrow(x) == 0L || ncol(x) == 0L)
    stop("'x' must have at least one row and one column.")

  x_cols <- if (is.matrix(x)) colnames(x) else names(x)
  if (is.null(x_cols) || any(is.na(x_cols)) || any(x_cols == ""))
    stop("All columns of 'x' must be named (no empty or NA column names).")

  x <- as.data.frame(x, check.names = FALSE)

  non_numeric <- vapply(x, function(col) !is.numeric(col), logical(1L))
  if (any(non_numeric))
    stop(
      "All columns of 'x' must be numeric. Non-numeric column(s): ",
      paste(x_cols[non_numeric], collapse = ", "), "."
    )

  # ---------------------------------------------------------------------------
  # 3. Validate metadata
  # ---------------------------------------------------------------------------
  if (!is.data.frame(metadata)) stop("'metadata' must be a data.frame.")
  if (nrow(metadata) != nrow(x))
    stop(
      "'metadata' must have the same number of rows as 'x'. ",
      "Got: metadata = ", nrow(metadata), " rows, x = ", nrow(x), " rows."
    )

  # ---------------------------------------------------------------------------
  # 4. Resolve show_x_as
  # ---------------------------------------------------------------------------
  if (is.null(show_x_as)) {
    meta_lower  <- tolower(names(metadata))
    id_patterns <- c("uid", "id", "sampleid")
    matched_idx <- NA_integer_
    for (pat in id_patterns) {
      hit <- which(meta_lower == pat)
      if (length(hit) > 0L) { matched_idx <- hit[1L]; break }
    }
    if (is.na(matched_idx))
      stop(
        "'show_x_as' is NULL and no column matching 'UID', 'ID', or ",
        "'SampleID' (case-insensitive) was found in 'metadata'. ",
        "Please supply a column name via 'show_x_as'."
      )
    show_x_as <- names(metadata)[matched_idx]
    message("Using '", show_x_as, "' from 'metadata' for the x-axis.")
  } else {
    if (!is.character(show_x_as) || length(show_x_as) != 1L)
      stop("'show_x_as' must be a single character string or NULL.")
    if (!show_x_as %in% names(metadata))
      stop(
        "Column '", show_x_as, "' not found in 'metadata'. ",
        "Available columns: ", paste(names(metadata), collapse = ", "), "."
      )
  }

  # ---------------------------------------------------------------------------
  # 5. Validate remaining scalar arguments
  # ---------------------------------------------------------------------------

  # --- group_brackets: single authoritative validation block ---
  if (!is.null(group_brackets)) {
    if (isTRUE(group_aggregate))
      stop("'group_brackets' can only be used when 'group_aggregate' is FALSE (plotting individual samples).")

    if (inherits(group_brackets, "formula")) {
      gb_formula <- group_brackets
    } else if (is.list(group_brackets) &&
               length(group_brackets) > 0L &&
               inherits(group_brackets[[1L]], "formula")) {
      gb_formula <- group_brackets[[1L]]
    } else {
      stop("'group_brackets' must be a formula (e.g., Group ~ c('Control', 'Case')) or a list containing one.")
    }

    gb_lhs <- all.vars(gb_formula[[2L]])
    if (length(gb_lhs) != 1L)
      stop("The left-hand side of 'group_brackets' must specify exactly one column name.")
    gb_col <- gb_lhs

    if (!gb_col %in% names(metadata))
      stop("Column '", gb_col, "' specified in 'group_brackets' not found in 'metadata'.")

    gb_rhs <- eval(gb_formula[[3L]], envir = parent.frame())
    if (!is.character(gb_rhs) && !is.factor(gb_rhs))
      stop("The right-hand side of 'group_brackets' must evaluate to a character or factor vector.")
    gb_rhs <- as.character(gb_rhs)
  }

  if (!is.null(limit)) {
    .check_pos_int(limit, "limit")
    limit <- as.integer(limit)
    if (limit > ncol(x)) {
      warning(
        "'limit' (", limit, ") exceeds the number of columns in 'x' (",
        ncol(x), "). All features will be plotted without 'Others'."
      )
      limit <- ncol(x)
    }
  }

  valid_sorts <- c("alphanumeric", "high", "low")
  if (!is.null(sort)) {
    if (!is.character(sort) || length(sort) != 1L || !sort %in% valid_sorts)
      stop("'sort' must be one of: ", paste(valid_sorts, collapse = ", "),
           ", or NULL. Got: '", sort, "'.")
  }

  if (!is.null(others_position)) {
    if (!is.character(others_position) || length(others_position) != 1L ||
        !others_position %in% c("first", "last"))
      stop("'others_position' must be 'first', 'last', or NULL.")
  }

  for (lname in c("group_aggregate", "normalize_to_relabund",
                  "ignore_normalization", "abundance_in_y",
                  "italicize_legend_items", "exclude_others",
                  "unknown_before_others")) {
    val <- get(lname)
    if (!is.logical(val) || length(val) != 1L)
      stop("'", lname, "' must be a single logical value (TRUE or FALSE).")
  }

  if (!is.null(order) && !identical(order, "alphanumeric")) {
    if (!is.character(order))
      stop("'order' must be a character vector, 'alphanumeric', or NULL.")
  }

  # angle default: 0 when aggregating by group, 45 for individual samples
  if (is.null(angle)) {
    angle <- if (isTRUE(group_aggregate)) 0L else 45L
  } else {
    if (!is.numeric(angle) || length(angle) != 1L || is.na(angle))
      stop("'angle' must be a single numeric value (0-360) or NULL.")
    angle <- angle %% 360
  }

  if (!is.null(legend_position)) {
    if (!is.character(legend_position) || length(legend_position) != 1L ||
        !legend_position %in% c("top", "bottom", "left", "right"))
      stop("'legend_position' must be 'top', 'bottom', 'left', 'right', or NULL.")
  }

  .check_pos_int(legend_nrow, "legend_nrow", allow_null = TRUE)
  .check_pos_int(legend_ncol, "legend_ncol", allow_null = TRUE)

  valid_themes <- c("nature", "minimal", "classic", "bw", "light", "dark")
  theme_choice <- tolower(theme)
  if (!theme_choice %in% valid_themes)
    stop("'theme' must be one of: ", paste(valid_themes, collapse = ", "),
         ". Got: '", theme, "'.")

  if (!is.character(font_family) || length(font_family) != 1L)
    stop("'font_family' must be a single character string.")

  for (sarg in list(
    list(v = plot_title,    n = "plot_title"),
    list(v = plot_subtitle, n = "plot_subtitle"),
    list(v = legend_title,  n = "legend_title"),
    list(v = ylab,          n = "ylab")
  )) {
    if (!is.character(sarg$v) || length(sarg$v) != 1L)
      stop("'", sarg$n, "' must be a single character string.")
  }
  if (!is.null(xlab) && (!is.character(xlab) || length(xlab) != 1L))
    stop("'xlab' must be a single character string or NULL.")

  .check_pos_num(global_font_size, "global_font_size", allow_null = TRUE)
  .check_pos_num(title_size,       "title_size",       allow_null = TRUE)
  .check_pos_num(subtitle_size,    "subtitle_size",    allow_null = TRUE)
  .check_pos_num(legend_size,      "legend_size",      allow_null = TRUE)
  .check_pos_num(xlab_size,        "xlab_size",        allow_null = TRUE)
  .check_pos_num(ylab_size,        "ylab_size",        allow_null = TRUE)

  # ---------------------------------------------------------------------------
  # 6. Normalisation checks (skipped when ignore_normalization = TRUE)
  # ---------------------------------------------------------------------------
  row_sums_check  <- rowSums(x, na.rm = TRUE)
  already_normed  <- all(row_sums_check == 0 |
                           abs(row_sums_check - 1) < .Machine$double.eps^0.5)

  if (!ignore_normalization) {
    if (!isTRUE(group_aggregate)) {
      not_norm_rows <- which(!(row_sums_check == 0 |
                                 abs(row_sums_check - 1) < .Machine$double.eps^0.5))
      if (length(not_norm_rows) > 0L)
        warning(
          "The following row(s) in 'x' do not sum to 0 or 1, suggesting ",
          "row-wise relative abundance normalisation has not been applied: ",
          paste(not_norm_rows, collapse = ", "), ". ",
          "Consider running run_normalize(x, metadata, ",
          "method = 'row_rel_abundance') first."
        )
    }

    if (!is.null(sort) && sort %in% c("high", "low") && already_normed) {
      warning(
        "'sort = \"", sort, "\"' ranks by column totals, which is only ",
        "meaningful for raw count data. Row sums of 'x' are already ",
        "normalised (0 or 1), so ranked sorting is not meaningful. ",
        "'sort' will be ignored. Supply raw counts with ",
        "'normalize_to_relabund = TRUE' to use ranked sorting."
      )
      sort <- NULL
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Optional pre-normalisation (row_rel_abundance)
  # ---------------------------------------------------------------------------
  needs_prenorm <- !ignore_normalization &&
    normalize_to_relabund &&
    !already_normed &&
    (isTRUE(group_aggregate) ||
       (!is.null(sort) && sort %in% c("high", "low")))

  if (needs_prenorm) {
    dimsprepr_ns   <- tryCatch(asNamespace("dimsprepr"), error = function(e) NULL)
    run_norm_fn <- if (!is.null(dimsprepr_ns))
      get0("run_normalize", envir = dimsprepr_ns, mode = "function")
    else NULL

    if (is.null(run_norm_fn)) {
      warning(
        "'normalize_to_relabund = TRUE' requires 'run_normalize' from the ",
        "'dimsprepr' package, which could not be found. Skipping pre-normalisation."
      )
    } else {
      norm_result <- run_norm_fn(
        x        = x,
        metadata = metadata,
        method   = "row_rel_abundance",
        verbose  = FALSE
      )
      x <- as.data.frame(norm_result$data, check.names = FALSE)
      message(
        "Row-wise relative abundance normalisation applied via ",
        "run_normalize() before aggregation/ranking."
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 8. Group aggregation
  # ---------------------------------------------------------------------------
  group_var <- as.character(metadata[[show_x_as]])

  if (isTRUE(group_aggregate)) {
    grp_levels <- unique(group_var)
    agg_mat    <- do.call(rbind, lapply(grp_levels, function(g) {
      rows <- which(group_var == g)
      vapply(x[rows, , drop = FALSE], base::sum, numeric(1L), na.rm = TRUE)
    }))
    x           <- as.data.frame(agg_mat, check.names = FALSE)
    colnames(x) <- x_cols
    group_var   <- grp_levels
  }

  n_samples <- nrow(x)

  # ---------------------------------------------------------------------------
  # 9. Feature selection by total abundance & handling sort = "low"
  # ---------------------------------------------------------------------------
  col_totals <- colSums(x, na.rm = TRUE)
  n_features <- length(col_totals)
  has_others <- FALSE

  if (!is.null(limit) && limit < n_features) {
    if (!is.null(sort) && sort == "low") {
      selected_idx      <- base::order(col_totals, decreasing = FALSE)[seq_len(limit)]
      selected_features <- x_cols[selected_idx]

      row_sums_full                      <- rowSums(x, na.rm = TRUE)
      row_sums_full[row_sums_full == 0]  <- NA
      rel_full <- x / row_sums_full

      rel_x         <- rel_full[, selected_features, drop = FALSE]
      row_sums_sub  <- rowSums(rel_x, na.rm = TRUE)
      row_sums_sub[row_sums_sub == 0] <- NA
      rel_x         <- rel_x / row_sums_sub

      message("With sort = 'low', selected low-abundance features have been re-normalised to 100% among themselves to maximize visibility.")
      has_others <- FALSE
    } else {
      selected_idx      <- base::order(col_totals, decreasing = TRUE)[seq_len(limit)]
      selected_features <- x_cols[selected_idx]
      other_cols        <- setdiff(x_cols, selected_features)

      row_sums_full                     <- rowSums(x, na.rm = TRUE)
      row_sums_full[row_sums_full == 0] <- NA
      rel_full <- x / row_sums_full

      rel_x      <- rel_full[, selected_features, drop = FALSE]
      others_vec <- rowSums(rel_full[, other_cols, drop = FALSE], na.rm = TRUE)

      if (!exclude_others) {
        rel_x[["Others"]] <- others_vec
        has_others <- TRUE
      }
    }
  } else {
    row_sums_full                     <- rowSums(x, na.rm = TRUE)
    row_sums_full[row_sums_full == 0] <- NA
    rel_x <- x / row_sums_full
  }

  # --- Apply exclude_unknown ---
  if (exclude_unknown) {
    unknown_cols <- colnames(rel_x)[tolower(colnames(rel_x)) == "unknown"]
    if (length(unknown_cols) > 0L)
      rel_x <- rel_x[, !colnames(rel_x) %in% unknown_cols, drop = FALSE]
  }

  all_na_rows <- which(rowSums(!is.na(rel_x)) == 0L)
  if (length(all_na_rows) > 0L)
    warning(
      "Row(s) ", paste(all_na_rows, collapse = ", "),
      " are entirely NA after relative abundance normalisation ",
      "(possible zero-sum rows)."
    )

  feat_names  <- colnames(rel_x)

  # ---------------------------------------------------------------------------
  # 10. Fill level ordering
  # ---------------------------------------------------------------------------
  non_others <- feat_names[feat_names != "Others"]

  fill_levels_base <- if (is.null(sort)) {
    x_cols[x_cols %in% non_others]
  } else if (sort == "alphanumeric") {
    base::sort(non_others)
  } else if (sort == "high") {
    names(base::sort(col_totals[non_others], decreasing = TRUE))
  } else {  # "low"
    names(base::sort(col_totals[non_others], decreasing = FALSE))
  }

  if (unknown_before_others) {
    unknown_hits     <- fill_levels_base[tolower(fill_levels_base) == "unknown"]
    fill_levels_base <- fill_levels_base[tolower(fill_levels_base) != "unknown"]
  } else {
    unknown_hits <- character(0L)
  }

  fill_levels <- if (has_others) {
    c(fill_levels_base, unknown_hits, "Others")
  } else {
    c(fill_levels_base, unknown_hits)
  }

  if (has_others && !is.null(others_position)) {
    fill_levels <- fill_levels[fill_levels != "Others"]
    fill_levels <- if (others_position == "first") {
      c("Others", fill_levels)
    } else {
      c(fill_levels, "Others")
    }
  }

  fill_levels <- fill_levels[fill_levels %in% colnames(rel_x)]

  # ---------------------------------------------------------------------------
  # 11. X-axis ordering
  # ---------------------------------------------------------------------------
  if (!is.null(group_brackets)) {
    current_grp_values <- unique(as.character(metadata[[gb_col]]))
    expected_order     <- gb_rhs[gb_rhs %in% current_grp_values]
    missing_from_rhs   <- setdiff(current_grp_values, expected_order)

    if (length(missing_from_rhs) > 0L)
      expected_order <- c(expected_order, missing_from_rhs)

    sample_groups_natural <- as.character(metadata[[gb_col]])
    runs_natural          <- rle(sample_groups_natural)

    if (identical(runs_natural$values, expected_order)) {
      if (!is.null(order) && !identical(order, "alphanumeric")) {
        missing_lev <- setdiff(order, group_var)
        if (length(missing_lev) > 0L)
          warning("Value(s) in 'order' not found in '", show_x_as, "': ",
                  paste(missing_lev, collapse = ", "), ". Ignored.")
        order     <- order[order %in% group_var]
        extra_lev <- setdiff(group_var, order)
        if (length(extra_lev) > 0L) {
          warning("Sample/group label(s) not in 'order' will be appended: ",
                  paste(extra_lev, collapse = ", "), ".")
          order <- c(order, extra_lev)
        }
        x_levels <- order
      } else if (identical(order, "alphanumeric")) {
        x_levels <- base::sort(unique(group_var))
      } else {
        x_levels <- unique(group_var)
      }
    } else {
      # Reorder samples so groups are contiguous per gb_rhs sequence
      sample_order_idx <- base::order(                    # explicit base::order — avoids collision with `order` param
        match(as.character(metadata[[gb_col]]), expected_order)
      )
      x_levels <- unique(as.character(metadata[[show_x_as]])[sample_order_idx])
    }

  } else {
    if (!is.null(order) && !identical(order, "alphanumeric")) {
      missing_lev <- setdiff(order, group_var)
      if (length(missing_lev) > 0L)
        warning("Value(s) in 'order' not found in '", show_x_as, "': ",
                paste(missing_lev, collapse = ", "), ". Ignored.")
      order     <- order[order %in% group_var]
      extra_lev <- setdiff(group_var, order)
      if (length(extra_lev) > 0L) {
        warning("Sample/group label(s) not in 'order' will be appended: ",
                paste(extra_lev, collapse = ", "), ".")
        order <- c(order, extra_lev)
      }
      x_levels <- order
    } else if (identical(order, "alphanumeric")) {
      x_levels <- base::sort(unique(group_var))
    } else {
      x_levels <- unique(group_var)
    }
  }

  # ---------------------------------------------------------------------------
  # 12. Build long-format data frame — manual pivot (supports special chars)
  # ---------------------------------------------------------------------------
  n_feats  <- length(feat_names)

  long_sample  <- rep(group_var,       times = n_feats)
  long_feature <- rep(feat_names, each = n_samples)
  long_abund   <- unlist(
    lapply(feat_names, function(fn) rel_x[[fn]]),
    use.names = FALSE
  )

  long_df <- data.frame(
    .sample    = factor(long_sample,  levels = x_levels),
    .feature   = factor(long_feature, levels = fill_levels),
    .abundance = long_abund,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  long_df <- long_df[!is.na(long_df[[".abundance"]]), , drop = FALSE]

  # ---------------------------------------------------------------------------
  # 13. Colour palette
  # ---------------------------------------------------------------------------
  okabe_base <- c(
    "#0072B2", "#D55E00", "#009E73", "#56B4E9",
    "#E69F00", "#CC79A7", "#F0E442", "#000000",
    "#44AA99", "#332288", "#DDCC77", "#882255",
    "#88CCEE", "#AA4499", "#117733", "#999933",
    "#661100", "#6699CC", "#888888", "#AA4466"
  )

  n_fill_levels <- length(fill_levels)

  if (!is.null(colors)) {
    if (!is.character(colors))
      stop("'colors' must be a character vector or NULL.")
    if (length(colors) != n_fill_levels)
      stop("'colors' must have length ", n_fill_levels,
           " (number of plotted fill levels). Got length ", length(colors), ".")
    fill_colors <- if (is.null(names(colors))) {
      stats::setNames(colors, fill_levels)
    } else {
      colors
    }
  } else {
    palette_vec <- if (n_fill_levels <= length(okabe_base)) {
      okabe_base[seq_len(n_fill_levels)]
    } else {
      c(okabe_base,
        grDevices::hcl.colors(n_fill_levels - length(okabe_base),
                              palette = "Dark 3"))
    }
    fill_colors <- stats::setNames(palette_vec, fill_levels)
    if ("Others" %in% fill_levels) fill_colors[["Others"]] <- "#BBBBBB"
  }

  # ---------------------------------------------------------------------------
  # 14. Font size resolution
  # ---------------------------------------------------------------------------
  .resolve_size <- function(explicit, scale_factor, base_default) {
    if (!is.null(explicit))     return(explicit)
    if (!is.null(scale_factor)) return(base_default * scale_factor)
    base_default
  }

  ref_base   <- 11
  gfs_factor <- if (!is.null(global_font_size)) global_font_size / ref_base else NULL

  r_title_size    <- .resolve_size(title_size,    gfs_factor, ref_base + 2)
  r_subtitle_size <- .resolve_size(subtitle_size, gfs_factor, ref_base + 1)
  r_legend_size   <- .resolve_size(legend_size,   gfs_factor, ref_base)
  r_xlab_size     <- .resolve_size(xlab_size,     gfs_factor, ref_base)
  r_ylab_size     <- .resolve_size(ylab_size,     gfs_factor, ref_base)
  r_axis_text     <- .resolve_size(NULL,          gfs_factor, ref_base - 1)
  r_legend_text   <- .resolve_size(NULL,          gfs_factor, ref_base - 1)
  r_strip_text    <- .resolve_size(NULL,          gfs_factor, ref_base)

  # ---------------------------------------------------------------------------
  # 15. Legend position and fill level order for legend
  # ---------------------------------------------------------------------------
  if (is.null(legend_position))
    legend_position <- if (abundance_in_y) "top" else "right"

  fill_levels_legend <- if (!abundance_in_y) rev(fill_levels) else fill_levels
  fill_colors_legend <- fill_colors[fill_levels_legend]

  # ---------------------------------------------------------------------------
  # 16. Axis labels — swap when abundance_in_y = FALSE
  # ---------------------------------------------------------------------------
  x_axis_label <- if (is.null(xlab)) show_x_as else xlab

  if (abundance_in_y) {
    labs_x <- x_axis_label
    labs_y <- ylab
  } else {
    labs_x <- ylab
    labs_y <- x_axis_label
  }

  # ---------------------------------------------------------------------------
  # 17. Base ggplot2 theme
  # ---------------------------------------------------------------------------
  base_thm <- switch(
    theme_choice,
    "nature"  = ggplot2::theme_bw(base_size = ref_base, base_family = font_family) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(color = "grey70", fill = NA,
                                                 linewidth = 0.6),
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
        strip.text       = ggplot2::element_text(face = "bold", family = font_family,
                                                 size = r_strip_text)
      ),
    "minimal" = ggplot2::theme_minimal(base_size = ref_base, base_family = font_family),
    "classic" = ggplot2::theme_classic(base_size = ref_base, base_family = font_family),
    "bw"      = ggplot2::theme_bw(base_size      = ref_base, base_family = font_family),
    "light"   = ggplot2::theme_light(base_size   = ref_base, base_family = font_family),
    "dark"    = ggplot2::theme_dark(base_size    = ref_base, base_family = font_family)
  )

  # ---------------------------------------------------------------------------
  # 18. Legend text element
  # ---------------------------------------------------------------------------
  if (italicize_legend_items) {
    if (requireNamespace("ggtext", quietly = TRUE)) {
      legend_labels_fmt <- stats::setNames(
        paste0("*", fill_levels_legend, "*"), fill_levels_legend
      )
      legend_text_el <- ggtext::element_markdown(size = r_legend_text,
                                                 family = font_family)
    } else {
      warning(
        "'italicize_legend_items = TRUE' requires 'ggtext'. ",
        "Falling back to element_text(face = 'italic'). ",
        "Install with: install.packages('ggtext')"
      )
      legend_labels_fmt <- NULL
      legend_text_el    <- ggplot2::element_text(face = "italic",
                                                 size = r_legend_text,
                                                 family = font_family)
    }
  } else {
    legend_labels_fmt <- NULL
    legend_text_el    <- ggplot2::element_text(size = r_legend_text,
                                               family = font_family)
  }

  # ---------------------------------------------------------------------------
  # 19. Axis text and title margin elements
  # ---------------------------------------------------------------------------
  hjust_val    <- if (angle == 0) 0.5 else 1
  angle_rad    <- angle * pi / 180
  extra_margin <- ceiling(r_axis_text * abs(sin(angle_rad)) * 0.8)

  text_x_el <- ggplot2::element_text(
    size   = r_axis_text,
    family = font_family,
    angle  = if (abundance_in_y) angle else 0,
    hjust  = if (abundance_in_y) hjust_val else 0.5,
    vjust  = 0.5
  )
  text_y_el <- ggplot2::element_text(
    size   = r_axis_text,
    family = font_family,
    angle  = if (!abundance_in_y) angle else 0,
    hjust  = if (!abundance_in_y) hjust_val else 1,
    vjust  = 0.5
  )

  axis_title_x_el <- ggplot2::element_text(
    size   = r_xlab_size,
    family = font_family,
    margin = if (abundance_in_y && angle != 0)
      ggplot2::margin(t = extra_margin, unit = "pt")
    else
      ggplot2::margin()
  )
  axis_title_y_el <- ggplot2::element_text(
    size   = r_ylab_size,
    family = font_family,
    margin = if (!abundance_in_y && angle != 0)
      ggplot2::margin(r = extra_margin, unit = "pt")
    else
      ggplot2::margin()
  )

  # ---------------------------------------------------------------------------
  # 20. Markdown title / subtitle
  # ---------------------------------------------------------------------------
  .uses_markdown <- function(txt) grepl("\\*|_|<[a-z]", txt, perl = TRUE)
  has_ggtext     <- requireNamespace("ggtext", quietly = TRUE)

  title_el <- if (.uses_markdown(plot_title) && has_ggtext)
    ggtext::element_markdown(size = r_title_size, hjust = 0.5, family = font_family)
  else
    ggplot2::element_text(size = r_title_size, hjust = 0.5, family = font_family,
                          face = "bold")

  subtitle_el <- if (.uses_markdown(plot_subtitle) && has_ggtext)
    ggtext::element_markdown(size = r_subtitle_size, hjust = 0.5, family = font_family)
  else
    ggplot2::element_text(size = r_subtitle_size, hjust = 0.5, family = font_family)

  # ---------------------------------------------------------------------------
  # 21. Assemble ggplot
  # ---------------------------------------------------------------------------
  legend_guide <- ggplot2::guide_legend(
    title  = legend_title,
    nrow   = legend_nrow,
    ncol   = legend_ncol,
    byrow  = TRUE
  )

  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x    = .data[[".sample"]],
      y    = .data[[".abundance"]],
      fill = .data[[".feature"]]
    )
  ) +
    ggplot2::geom_bar(
      stat     = "identity",
      position = "stack",
      width    = 0.8,
      ...
    ) +
    ggplot2::scale_fill_manual(
      values = fill_colors_legend,
      limits = fill_levels_legend,
      labels = if (!is.null(legend_labels_fmt)) legend_labels_fmt else ggplot2::waiver(),
      guide  = legend_guide,
      drop   = FALSE
    ) +
    ggplot2::scale_x_discrete(limits = x_levels) +
    ggplot2::scale_y_continuous(
      limits = if (!is.null(group_brackets)) c(-0.15, 1.02) else NULL,
      expand = if (!is.null(group_brackets)) c(0, 0) else ggplot2::expansion(mult = c(0, 0.02)),
      labels = function(v) {
        if (!is.null(group_brackets)) {
          ifelse(v >= 0, paste0(round(v * 100, 0), "%"), "")
        } else {
          paste0(round(v * 100, 0), "%")
        }
      }
    ) +
    ggplot2::labs(
      title    = plot_title,
      subtitle = plot_subtitle,
      x        = labs_x,
      y        = labs_y,
      fill     = legend_title
    ) +
    base_thm +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = title_el,
      plot.subtitle    = subtitle_el,
      axis.title.x     = axis_title_x_el,
      axis.title.y     = axis_title_y_el,
      axis.text        = ggplot2::element_text(size = r_axis_text,
                                               family = font_family),
      axis.text.x      = text_x_el,
      axis.text.y      = text_y_el,
      legend.position  = legend_position,
      legend.title     = ggplot2::element_text(size = r_legend_size,
                                               family = font_family,
                                               face = "bold"),
      legend.text      = legend_text_el
    )

  if (!abundance_in_y) p <- p + ggplot2::coord_flip()

  # ---------------------------------------------------------------------------
  # 21.5. Render group brackets (if supplied)
  # ---------------------------------------------------------------------------
  if (!is.null(group_brackets)) {
    sample_to_group <- data.frame(
      sample = x_levels,
      group  = as.character(
        metadata[[gb_col]][match(x_levels, as.character(metadata[[show_x_as]]))]
      ),
      stringsAsFactors = FALSE
    )

    runs      <- rle(sample_to_group$group)
    end_pos   <- cumsum(runs$lengths)
    start_pos <- end_pos - runs$lengths + 1

    bracket_df <- data.frame(
      group = runs$values,
      start = start_pos,
      end   = end_pos,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(bracket_df))) {
      b_start <- bracket_df$start[i] - 0.45
      b_end   <- bracket_df$end[i]   + 0.45
      b_mid   <- (b_start + b_end) / 2
      b_label <- bracket_df$group[i]

      if (abundance_in_y) {
        p <- p +
          ggplot2::annotate("segment", x = b_start, xend = b_end,   y = -0.06, yend = -0.06, color = "black", linewidth = 0.6) +
          ggplot2::annotate("segment", x = b_start, xend = b_start, y = -0.06, yend = -0.03, color = "black", linewidth = 0.6) +
          ggplot2::annotate("segment", x = b_end,   xend = b_end,   y = -0.06, yend = -0.03, color = "black", linewidth = 0.6) +
          ggplot2::annotate("text",    x = b_mid,   y = -0.11,      label = b_label,
                            fontface = "bold", size = r_axis_text * 0.353,
                            family = font_family, vjust = 0.5)
      } else {
        p <- p +
          ggplot2::annotate("segment", x = b_start, xend = b_end,   y = -0.06, yend = -0.06, color = "black", linewidth = 0.6) +
          ggplot2::annotate("segment", x = b_start, xend = b_start, y = -0.06, yend = -0.03, color = "black", linewidth = 0.6) +
          ggplot2::annotate("segment", x = b_end,   xend = b_end,   y = -0.06, yend = -0.03, color = "black", linewidth = 0.6) +
          ggplot2::annotate("text",    x = b_mid,   y = -0.11,      label = b_label,
                            fontface = "bold", size = r_axis_text * 0.353,
                            family = font_family, angle = 90, vjust = 0.5)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 22. Draw and return
  # ---------------------------------------------------------------------------
  print(p)
  message("Relative abundance plot created successfully.")
  return(invisible(p))
}