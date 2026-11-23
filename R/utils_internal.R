# Internal helpers shared by the plot_* and run_* functions. Not exported.

# ---- Plot text --------------------------------------------------------------

# Text sizes in points. Everything follows `base_size` unless given explicitly,
# so one number rescales a whole figure. The default of 11 suits the 7 x 5 in
# figures of a default Quarto or R Markdown document.
.text_sizes <- function(base_size = 11, title = NULL, subtitle = NULL,
                        xlab = NULL, ylab = NULL, axis_text = NULL,
                        legend = NULL, caption = NULL, strip = NULL,
                        label = NULL) {
  args <- list(base_size = base_size, title_size = title,
               subtitle_size = subtitle, xlab_size = xlab, ylab_size = ylab,
               axis_text_size = axis_text, legend_size = legend,
               caption_size = caption, strip_size = strip, label_size = label)
  for (nm in names(args)) {
    v <- args[[nm]]
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0))
      stop("`", nm, "` must be a single positive number.", call. = FALSE)
  }
  list(
    base      = base_size,
    title     = title     %||% (base_size * 1.2),
    subtitle  = subtitle  %||% base_size,
    xlab      = xlab      %||% base_size,
    ylab      = ylab      %||% base_size,
    axis_text = axis_text %||% (base_size * 0.9),
    legend    = legend    %||% (base_size * 0.9),
    caption   = caption   %||% (base_size * 0.8),
    strip     = strip     %||% (base_size * 0.9),
    label     = label     %||% (base_size * 0.8)
  )
}

# Text geoms take their size in millimetres, not points.
.pt_to_mm <- function(pt) pt / ggplot2::.pt

# The theme used by every plot_* function; `sizes` comes from .text_sizes().
.ponder_theme <- function(theme = "nature", sizes = .text_sizes(),
                          font_family = "") {
  valid <- c("nature", "classic", "plos", "minimal", "bw", "light", "dark", "gray")
  if (!is.character(theme) || length(theme) != 1L || !tolower(theme) %in% valid)
    stop("`theme` must be one of: ", paste(valid, collapse = ", "), ".",
         call. = FALSE)
  theme <- tolower(theme)
  base_fn <- switch(theme,
    nature = , bw = ggplot2::theme_bw,
    classic = , plos = ggplot2::theme_classic,
    minimal = ggplot2::theme_minimal,
    light = ggplot2::theme_light,
    dark = ggplot2::theme_dark,
    gray = ggplot2::theme_gray
  )
  out <- base_fn(base_size = sizes$base, base_family = font_family) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = sizes$title, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = sizes$subtitle, colour = "grey30"),
      plot.caption  = ggplot2::element_text(size = sizes$caption, colour = "grey40"),
      axis.title.x  = ggplot2::element_text(size = sizes$xlab),
      axis.title.y  = ggplot2::element_text(size = sizes$ylab),
      axis.text     = ggplot2::element_text(size = sizes$axis_text),
      legend.title  = ggplot2::element_text(size = sizes$legend, face = "bold"),
      legend.text   = ggplot2::element_text(size = sizes$legend),
      strip.text    = ggplot2::element_text(size = sizes$strip, face = "bold")
    )
  if (theme == "nature")
    out <- out + ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70")
    )
  out
}

# ---- Report tables ----------------------------------------------------------

# Report-ready headers for test-result tables; drops columns that only repeat
# another column as text.
.report_names <- function(d) {
  hdr <- c(group = "Group", mean = "Mean", sd = "SD", se = "SE", median = "Median",
           q25 = "Q1", q75 = "Q3", min = "Min", max = "Max", check = "Check",
           result = "Result", p_value = "p-value", decision = "Decision",
           group1 = "Group 1", group2 = "Group 2", n1 = "n (group 1)",
           n2 = "n (group 2)", mean1 = "Mean (group 1)", mean2 = "Mean (group 2)",
           statistic = "Statistic", estimate = "Estimate", conf.low = "95% CI lower",
           conf.high = "95% CI upper", p = "p-value", p.adj = "Adjusted p-value",
           p.adj.signif = "Significance", posthoc_test = "Test",
           significant = "Significant", interpretation = "Interpretation",
           effect = "Effect", metric = "Metric", ci_low = "95% CI lower",
           ci_high = "95% CI upper", magnitude = "Magnitude")
  d <- as.data.frame(d)
  d$p.format <- d$p.adj.format <- d$.y. <- d$method <- d$alternative <- NULL
  hit <- names(d) %in% names(hdr)
  names(d)[hit] <- hdr[names(d)[hit]]
  d
}

# ---- dimsprepr input --------------------------------------------------------

.is_di <- function(x) inherits(x, c("run_DIpreprocess", "run_DIpreprocess_multi"))

# Unpacks a dimsprepr::run_DIpreprocess() result. Its metadata always uses the
# standard column names, whatever the original group column was called.
# stage: "scaled" (final non-PLS), "pls" (final PLS), or an earlier step.
.di_unpack <- function(x, stage = "scaled") {
  if (inherits(x, "run_DIpreprocess_multi"))
    stop("`x` holds several preprocessing combinations. Pick one first, e.g. ",
         "x[[\"", names(x)[1L], "\"]]. Available: ",
         paste(names(x), collapse = ", "), ".", call. = FALSE)
  if (!is.null(x$error))
    stop("`x` is a failed preprocessing run: ", x$error, call. = FALSE)
  # Only the final scaled data have merged-replicate versions
  use_merged <- !is.null(x$metadata_merged) && stage %in% c("scaled", "pls")
  data <- switch(stage,
    scaled      = if (use_merged) x$data_nonpls_merged else x$data_nonpls,
    pls         = if (use_merged) x$data_pls_merged else x$data_pls,
    normalized  = x$data_normalized,
    corrected   = x$data_corrected,
    transformed = x$data_transformed
  )
  if (is.null(data))
    stop("Data stage '", stage, "' is not available in `x`.", call. = FALSE)
  list(
    data     = data,
    metadata = if (use_merged) x$metadata_merged else x$metadata,
    group    = "Group",
    qc       = x$parameters$qc_types %||% c("QC", "SQC", "EQC")
  )
}
