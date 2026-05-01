#' @title Principal Coordinate Analysis (PCoA) with Optional PERMANOVA and Post-hoc
#'
#' @description
#' Performs Principal Coordinate Analysis (PCoA), also known as metric
#' multidimensional scaling (MDS), as a sequential wrapper around
#' \code{\link[vegan]{vegdist}} and \code{\link[ape]{pcoa}}. Optionally runs a
#' Permutational Multivariate Analysis of Variance (PERMANOVA) via
#' \code{\link[vegan]{adonis2}} on the same dissimilarity matrix. Can also perform
#' pairwise post-hoc testing using \code{OmicFlow::pairwise_adonis}. Both
#' analyses share a single dissimilarity matrix so the method is applied exactly
#' once.
#'
#' @param x A numeric \code{data.frame}, \code{tibble}, or \code{matrix} of
#'   observations (rows) by variables (columns). Column names may contain
#'   special characters. Alternatively, a precomputed dissimilarity object of
#'   class \code{"dist"}; when \code{x} is a \code{"dist"} object the
#'   \code{method} argument is ignored and the matrix is used directly for both
#'   PCoA and PERMANOVA.
#' @param method A single character string specifying the dissimilarity index
#'   passed to \code{\link[vegan]{vegdist}}. Common choices include
#'   \code{"bray"} (Bray–Curtis), \code{"euclidean"}, \code{"aitchison"}, and
#'   \code{"robust.aitchison"}; see \code{?vegan::vegdist} for the full list.
#'   Ignored when \code{x} is already a \code{"dist"} object. The same
#'   dissimilarity matrix is reused for PERMANOVA — the method is never applied
#'   twice. Default: \code{"bray"}.
#' @param correction Correction method for negative eigenvalues passed to
#'   \code{\link[ape]{pcoa}}. One of \code{"lingoes"}, \code{"cailliez"}, or
#'   \code{"none"}. Default: \code{"none"}.
#' @param metadata A \code{data.frame} with one row per observation, in the
#'   same row order as \code{x}. Required when \code{rhs} is specified
#'   (PERMANOVA); ignored otherwise. The grouping variable(s) named in
#'   \code{rhs} must be columns of \code{metadata}.
#' @param rhs A single character string naming a column of \code{metadata} to
#'   use as the right-hand side of the PERMANOVA formula
#'   (\code{dist_matrix ~ rhs}). When \code{NULL} (default), PERMANOVA is
#'   skipped entirely.
#' @param perms A positive integer. Number of permutations for PERMANOVA.
#'   Passed as the \code{permutations} argument to \code{\link[vegan]{adonis2}}.
#'   Default: \code{9999}.
#' @param groups A single character string naming a column of \code{metadata} to
#'   use for pairwise post-hoc testing via \code{OmicFlow::pairwise_adonis}.
#'   When \code{NULL} (default), post-hoc analysis is skipped.
#' @param p_adjust A single character string specifying the p-value adjustment
#'   method for post-hoc testing. Passed to \code{\link[stats]{p.adjust}}.
#'   Default: \code{"BH"} (Benjamini-Hochberg).
#' @param ... Additional arguments forwarded to \code{\link[vegan]{vegdist}}
#'   (e.g., \code{binary}, \code{na.rm}) when \code{x} is not a \code{"dist"}
#'   object, to \code{\link[ape]{pcoa}} (e.g., \code{rn}), and to
#'   \code{\link[vegan]{adonis2}} (e.g., \code{by}, \code{parallel},
#'   \code{strata}) when PERMANOVA is requested, and to
#'   \code{\link[OmicFlow]{pairwise_adonis}}. Arguments are routed to each
#'   function by matching against their respective \code{formals()}; unrecognised
#'   names trigger a warning.
#'
#' @details
#' ## Algorithm
#' When \code{x} is a numeric matrix or data frame the function proceeds as:
#' \enumerate{
#'   \item \strong{Dissimilarity}: \code{vegan::vegdist(x, method = method, ...)}
#'     produces a \code{"dist"} object.
#'   \item \strong{PCoA}: \code{ape::pcoa(dist_matrix, correction = correction, ...)}
#'     decomposes the matrix into principal coordinates.
#'   \item \strong{PERMANOVA} (when \code{rhs} is not \code{NULL}):
#'     \code{vegan::adonis2(dist_matrix ~ rhs, data = metadata, ...)} tests
#'     group differences using the already-computed dissimilarity matrix.
#'     Because the \code{"dist"} object is passed directly, \code{method} is
#'     not re-applied.
#' }
#' When \code{x} is already a \code{"dist"} object, stage 1 is skipped.
#'
#' ## PERMANOVA
#' PERMANOVA (Anderson, 2001; McArdle & Anderson, 2001) partitions multivariate
#' variance using a distance matrix and tests significance by permuting rows.
#' It is robust to non-normality but sensitive to heterogeneous within-group
#' dispersions (Warton et al., 2012); use alongside a test of homogeneity of
#' dispersions (\code{vegan::betadisper}) when groups differ in spread.
#'
#' Key output fields from \code{adonis2}:
#' \itemize{
#'   \item \strong{R²}: proportion of total variance explained by the grouping
#'     factor.
#'   \item \strong{F}: pseudo-F statistic.
#'   \item \strong{Pr(>F)}: permutation p-value.
#' }
#'
#' ## Negative eigenvalues
#' Non-Euclidean dissimilarity coefficients may produce negative eigenvalues.
#' Apply \code{correction = "lingoes"} or \code{"cailliez"} to obtain a fully
#' Euclidean representation (Lingoes, 1971; Cailliez, 1983).
#'
#' ## Variance explained
#' Computed as \eqn{\lambda_i / \sum \lambda^+}{lambda_i / sum(positive lambdas)},
#' using only positive eigenvalues, consistent with \pkg{ape} convention.
#'
#' ## Compatibility with \code{plot_score}
#' The returned object carries class \code{"run_pcoa"} and a \code{scores}
#' element with columns \code{PC1}, \code{PC2}, \ldots so that
#' \code{pondeR::plot_score} can dispatch correctly. When PERMANOVA results are
#' present, \code{plot_score.run_pcoa} formats and injects them automatically
#' as the plot subtitle (R² and p-value).
#'
#' @return A named list of class \code{c("run_pcoa", "list")} with the following
#'   elements:
#' \describe{
#'   \item{\code{scores}}{A \code{data.frame} of principal coordinate scores,
#'     columns \code{PC1}, \code{PC2}, \ldots Row names match those of \code{x}.}
#'   \item{\code{variance_explained}}{Named numeric vector of percentage variance
#'     explained per axis (\code{PC1}, \code{PC2}, \ldots).}
#'   \item{\code{eigenvalues}}{Raw eigenvalue vector from \code{ape::pcoa}.}
#'   \item{\code{dist_matrix}}{The \code{"dist"} object used for both PCoA and
#'     PERMANOVA.}
#'   \item{\code{method}}{Dissimilarity method used, or \code{"user-supplied dist"}.}
#'   \item{\code{correction}}{Correction method applied.}
#'   \item{\code{n_obs}}{Integer. Number of observations.}
#'   \item{\code{n_axes}}{Integer. Number of positive principal coordinate axes.}
#'   \item{\code{pcoa_raw}}{Complete \code{ape::pcoa} output for downstream
#'     compatibility.}
#'   \item{\code{vegdist_raw}}{\code{"dist"} object from \code{vegan::vegdist},
#'     or \code{NULL} when \code{x} was already a \code{"dist"} object.}
#'   \item{\code{permanova}}{When \code{rhs} is not \code{NULL}: a list with
#'     elements \code{table} (the full \code{adonis2} result as a
#'     \code{data.frame}), \code{R2} (numeric, proportion of variance for the
#'     grouping term), \code{F_stat} (numeric, pseudo-F), \code{p_value}
#'     (numeric, permutation p-value), \code{rhs} (character, the grouping
#'     variable name), and \code{perms} (integer, number of permutations used).
#'     \code{NULL} when PERMANOVA was not requested.}
#'   \item{\code{pairwise_adonis}}{When \code{groups} is specified: the full
#'     result table from \code{OmicFlow::pairwise_adonis}.}
#'   \item{\code{call}}{The matched function call.}
#' }
#'
#' @examples
#' ## ── Example 1: Bray–Curtis PCoA on the dune dataset ──────────────────────
#' if (requireNamespace("vegan", quietly = TRUE) &&
#'     requireNamespace("ape",   quietly = TRUE)) {
#'
#'   data(dune,     package = "vegan")
#'   data(dune.env, package = "vegan")
#'
#'   result <- run_pcoa(dune, method = "bray")
#'   print(result)
#'   summary(result)
#'   head(result$scores[, c("PC1", "PC2")])
#'   result$variance_explained[1:3]
#' }
#'
#' ## ── Example 2: PCoA + PERMANOVA ───────────────────────────────────────────
#' if (requireNamespace("vegan", quietly = TRUE) &&
#'     requireNamespace("ape",   quietly = TRUE)) {
#'
#'   data(dune,     package = "vegan")
#'   data(dune.env, package = "vegan")
#'
#'   result_perm <- run_pcoa(
#'     x        = dune,
#'     method   = "bray",
#'     metadata = dune.env,
#'     rhs      = "Management",
#'     perms    = 999
#'   )
#'
#'   # Extract PERMANOVA results individually
#'   result_perm$permanova$R2
#'   result_perm$permanova$p_value
#'   result_perm$permanova$table
#' }
#'
#' ## ── Example 3: Lingoes correction ─────────────────────────────────────────
#' if (requireNamespace("vegan", quietly = TRUE) &&
#'     requireNamespace("ape",   quietly = TRUE)) {
#'
#'   data(dune, package = "vegan")
#'   result_ling <- run_pcoa(dune, method = "bray", correction = "lingoes")
#'   result_ling$variance_explained[1:3]
#' }
#'
#' ## ── Example 4: Precomputed dist + PERMANOVA ───────────────────────────────
#' if (requireNamespace("vegan", quietly = TRUE) &&
#'     requireNamespace("ape",   quietly = TRUE)) {
#'
#'   data(dune,     package = "vegan")
#'   data(dune.env, package = "vegan")
#'   d <- vegan::vegdist(dune, method = "jaccard")
#'
#'   result_dist <- run_pcoa(
#'     x        = d,
#'     metadata = dune.env,
#'     rhs      = "Management",
#'     perms    = 999
#'   )
#'   result_dist$method          # "user-supplied dist"
#'   result_dist$permanova$R2
#' }
#'
#' @references
#' Anderson, M.J. (2001) A new method for non-parametric multivariate analysis
#'   of variance. \emph{Austral Ecology}, \strong{26}, 32–46.
#'
#' Cailliez, F. (1983) The analytical solution of the additive constant
#'   problem. \emph{Psychometrika}, \strong{48}, 305–308.
#'
#' Excoffier, L., Smouse, P.E., and Quattro, J.M. (1992) Analysis of molecular
#'   variance inferred from metric distances among DNA haplotypes: Application
#'   to human mitochondrial DNA restriction data. \emph{Genetics},
#'   \strong{131}, 479–491.
#'
#' Gower, J.C. (1966) Some distance properties of latent root and vector
#'   methods used in multivariate analysis. \emph{Biometrika}, \strong{53},
#'   325–338.
#'
#' Gower, J.C. and Legendre, P. (1986) Metric and Euclidean properties of
#'   dissimilarity coefficients. \emph{Journal of Classification}, \strong{3},
#'   5–48.
#'
#' Legendre, P. and Anderson, M.J. (1999) Distance-based redundancy analysis:
#'   Testing multispecies responses in multifactorial ecological experiments.
#'   \emph{Ecological Monographs}, \strong{69}, 1–24.
#'
#' Legendre, P. and Gallagher, E.D. (2001) Ecologically meaningful
#'   transformations for ordination of species data. \emph{Oecologia},
#'   \strong{129}, 271–280.
#'
#' Legendre, P. and Legendre, L. (1998) \emph{Numerical Ecology}, 2nd English
#'   edition. Amsterdam: Elsevier Science BV.
#'
#' Lingoes, J.C. (1971) Some boundary conditions for a monotone analysis of
#'   symmetric matrices. \emph{Psychometrika}, \strong{36}, 195–203.
#'
#' McArdle, B.H. and Anderson, M.J. (2001) Fitting multivariate models to
#'   community data: A comment on distance-based redundancy analysis.
#'   \emph{Ecology}, \strong{82}, 290–297.
#'
#' Warton, D.I., Wright, T.W., and Wang, Y. (2012) Distance-based multivariate
#'   analyses confound location and dispersion effects. \emph{Methods in
#'   Ecology and Evolution}, \strong{3}, 89–101.
#'
#' @author John Lennon L. Calorio
#' @export
run_pcoa <- function(x,
                     method     = "bray",
                     correction = "none",
                     metadata   = NULL,
                     rhs        = NULL,
                     perms      = 9999L,
                     groups     = NULL,
                     p_adjust   = "BH",
                     ...) {

  # ── 0. Capture call ──────────────────────────────────────────────────────────
  .call <- match.call()

  # ── 1. Check required packages ───────────────────────────────────────────────
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop(
      "'vegan' is required for run_pcoa() but is not installed.\n",
      "  Install it with: install.packages(\"vegan\")",
      call. = FALSE
    )
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop(
      "'ape' is required for run_pcoa() but is not installed.\n",
      "  Install it with: install.packages(\"ape\")",
      call. = FALSE
    )
  }

  # ── 2. Validate scalar arguments ─────────────────────────────────────────────
  correction <- match.arg(correction, choices = c("none", "lingoes", "cailliez"))

  if (!is.null(rhs)) {
    if (!is.character(rhs) || length(rhs) != 1L || is.na(rhs)) {
      stop(
        "'rhs' must be a single non-NA character string naming a column of ",
        "'metadata'.",
        call. = FALSE
      )
    }
    if (!is.numeric(perms) || length(perms) != 1L || perms < 1L) {
      stop(
        "'perms' must be a single positive integer (e.g., 999 or 9999).",
        call. = FALSE
      )
    }
    perms <- as.integer(perms)
  }

  # ── 3. Validate metadata when PERMANOVA is requested ─────────────────────────
  run_permanova <- !is.null(rhs)

  if (run_permanova) {
    if (is.null(metadata)) {
      stop(
        "'metadata' must be provided when 'rhs' is specified.\n",
        "  Supply a data.frame with one row per observation matching 'x'.",
        call. = FALSE
      )
    }
    if (!is.data.frame(metadata) && !inherits(metadata, "tbl_df")) {
      stop(
        "'metadata' must be a data.frame or tibble.\n",
        "  Received: ", paste(class(metadata), collapse = ", "),
        call. = FALSE
      )
    }
  }

  # ── 4. Branch on whether x is already a "dist" object ────────────────────────
  is_precomputed_dist <- inherits(x, "dist")

  if (is_precomputed_dist) {

    if (!missing(method)) {
      warning(
        "'x' is already a \"dist\" object; the 'method' argument is ignored.",
        call. = FALSE
      )
    }
    dist_matrix <- x
    vegdist_raw <- NULL
    method_used <- "user-supplied dist"

    # Validate metadata row count against dist size
    if (run_permanova) {
      n_dist <- attr(dist_matrix, "Size")
      if (nrow(metadata) != n_dist) {
        stop(
          "'metadata' has ", nrow(metadata), " row(s) but the distance matrix ",
          "has ", n_dist, " observation(s). Row counts must match.",
          call. = FALSE
        )
      }
      if (!rhs %in% colnames(metadata)) {
        stop(
          "Column '", rhs, "' specified in 'rhs' not found in 'metadata'.\n",
          "  Available columns: ",
          paste(colnames(metadata), collapse = ", "),
          call. = FALSE
        )
      }
    }

  } else {

    # ── 4b. x is a matrix / data.frame ──────────────────────────────────────────
    if (!is.data.frame(x) && !is.matrix(x) && !inherits(x, "tbl_df")) {
      stop(
        "'x' must be a numeric data.frame, tibble, matrix, or an object of ",
        "class \"dist\".\n",
        "  Received: ", paste(class(x), collapse = ", "),
        call. = FALSE
      )
    }

    x_mat <- as.matrix(x)

    if (nrow(x_mat) < 2L) {
      stop(
        "'x' must have at least 2 rows (observations).\n",
        "  'x' has ", nrow(x_mat), " row(s).",
        call. = FALSE
      )
    }
    if (ncol(x_mat) < 1L) {
      stop("'x' must have at least 1 column.", call. = FALSE)
    }
    if (!is.numeric(x_mat)) {
      stop(
        "'x' must be numeric. Non-numeric columns detected.\n",
        "  Convert or remove them before calling run_pcoa().",
        call. = FALSE
      )
    }
    if (any(is.nan(x_mat))) {
      stop(
        "'x' contains NaN values. Replace or remove them before calling ",
        "run_pcoa().",
        call. = FALSE
      )
    }
    if (any(is.infinite(x_mat))) {
      stop(
        "'x' contains infinite values. Replace or remove them before calling ",
        "run_pcoa().",
        call. = FALSE
      )
    }
    if (anyNA(x_mat)) {
      warning(
        "'x' contains NA values. vegan::vegdist() will handle them according ",
        "to its 'na.rm' argument (default FALSE). Consider imputing or ",
        "removing missing values first.",
        call. = FALSE
      )
    }
    if (method %in% c("bray", "kulczynski", "canberra", "clark",
                      "binomial", "cao", "chao") && any(x_mat < 0)) {
      stop(
        "Method \"", method, "\" requires non-negative values, but 'x' ",
        "contains negative entries.\n",
        "  Consider transforming your data or choosing a different method.",
        call. = FALSE
      )
    }
    if (method == "aitchison" && any(x_mat <= 0)) {
      stop(
        "Method \"aitchison\" requires strictly positive values, but 'x' ",
        "contains values <= 0.\n",
        "  Add a pseudocount (e.g., x + 0.5) or use method = ",
        "\"robust.aitchison\".",
        call. = FALSE
      )
    }
    if (!is.character(method) || length(method) != 1L || is.na(method)) {
      stop(
        "'method' must be a single non-NA character string.\n",
        "  See ?vegan::vegdist for available options.",
        call. = FALSE
      )
    }

    # Validate metadata dimensions when PERMANOVA requested
    if (run_permanova) {
      if (nrow(metadata) != nrow(x_mat)) {
        stop(
          "'metadata' has ", nrow(metadata), " row(s) but 'x' has ",
          nrow(x_mat), " row(s). They must match.",
          call. = FALSE
        )
      }
      if (!rhs %in% colnames(metadata)) {
        stop(
          "Column '", rhs, "' specified in 'rhs' not found in 'metadata'.\n",
          "  Available columns: ",
          paste(colnames(metadata), collapse = ", "),
          call. = FALSE
        )
      }
    }

    # ── Separate ... for vegdist / pcoa / adonis2 ──────────────────────────────
    dots      <- list(...)
    vd_args   <- .split_dots_vegdist(dots)
    pcoa_args <- .split_dots_pcoa(dots)

    known_nms <- Reduce(
      union,
      list(
        names(formals(vegan::vegdist)),
        names(formals(ape::pcoa)),
        names(formals(vegan::adonis2))
      )
    )
    unknown <- setdiff(names(dots), known_nms)
    if (length(unknown) > 0L) {
      warning(
        "The following '...' argument(s) were not recognised by vegdist(), ",
        "pcoa(), or adonis2() and will be ignored:\n  ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }

    # ── Compute dissimilarity (once) ───────────────────────────────────────────
    vegdist_raw <- tryCatch(
      do.call(vegan::vegdist, c(list(x = x_mat, method = method), vd_args)),
      error = function(e) {
        stop(
          "vegan::vegdist() failed:\n  ", conditionMessage(e), "\n",
          "  Check that 'method' is valid and 'x' meets its requirements.",
          call. = FALSE
        )
      }
    )
    dist_matrix <- vegdist_raw
    method_used <- method
  }

  # ── 5. Separate ... for pcoa when x was a dist ───────────────────────────────
  if (is_precomputed_dist) {
    dots      <- list(...)
    pcoa_args <- .split_dots_pcoa(dots)
  }

  # ── 6. Run PCoA ───────────────────────────────────────────────────────────────
  pcoa_raw <- tryCatch(
    do.call(ape::pcoa, c(list(D = dist_matrix, correction = correction), pcoa_args)),
    error = function(e) {
      stop(
        "ape::pcoa() failed:\n  ", conditionMessage(e), "\n",
        "  If negative eigenvalues are the cause, try correction = ",
        "\"lingoes\" or \"cailliez\".",
        call. = FALSE
      )
    }
  )

  # ── 7. Extract and rename scores ──────────────────────────────────────────────
  raw_scores <- if (!is.null(pcoa_raw$vectors.cor)) pcoa_raw$vectors.cor else pcoa_raw$vectors

  if (is.null(raw_scores) || ncol(raw_scores) == 0L) {
    stop(
      "ape::pcoa() returned no principal coordinate axes.\n",
      "  The dissimilarity matrix may be degenerate or all eigenvalues ",
      "non-positive.\n",
      "  Try correction = \"lingoes\" or \"cailliez\".",
      call. = FALSE
    )
  }

  n_axes <- ncol(raw_scores)
  pc_nms <- paste0("PC", seq_len(n_axes))
  scores <- as.data.frame(raw_scores, stringsAsFactors = FALSE)
  colnames(scores) <- pc_nms

  if (!is.null(attr(dist_matrix, "Labels"))) {
    rownames(scores) <- attr(dist_matrix, "Labels")
  }

  # ── 8. Variance explained ─────────────────────────────────────────────────────
  eig_values <- .extract_eigenvalues(pcoa_raw, correction)
  pos_eig    <- eig_values[eig_values > 0]
  var_pct    <- eig_values / sum(pos_eig) * 100
  var_pct    <- var_pct[seq_len(n_axes)]
  names(var_pct) <- pc_nms

  # ── 9. Warn on negative eigenvalues (no correction) ──────────────────────────
  all_eig <- pcoa_raw$values$Eigenvalues
  if (correction == "none" && any(all_eig < 0, na.rm = TRUE)) {
    n_neg <- sum(all_eig < 0, na.rm = TRUE)
    warning(
      n_neg, " negative eigenvalue(s) detected in the PCoA decomposition.\n",
      "  This is common with non-Euclidean dissimilarity indices (e.g., ",
      "Bray-Curtis).\n",
      "  Consider correction = \"lingoes\" or \"cailliez\" for a fully ",
      "Euclidean representation.",
      call. = FALSE
    )
  }

  n_obs <- nrow(scores)

  # ── 10. PERMANOVA via vegan::adonis2 ─────────────────────────────────────────
  permanova_result <- NULL

  if (run_permanova) {
    # Build a minimal data frame that adonis2 can resolve the formula against.
    # We bind the grouping column from metadata; the LHS is the dist object
    # passed directly so 'x' (the raw feature matrix) is NOT re-used here —
    # the dissimilarity is never recomputed.
    perm_data <- data.frame(
      .grp = metadata[[rhs]],
      stringsAsFactors = FALSE
    )
    names(perm_data) <- rhs   # restore original column name for formula

    perm_formula <- stats::as.formula(paste0("dist_matrix ~ `", rhs, "`"))

    # Route additional ... to adonis2 only
    adonis_args <- if (!is_precomputed_dist) .split_dots_adonis2(dots) else {
      all_dots <- list(...)
      .split_dots_adonis2(all_dots)
    }

    adonis_out <- tryCatch(
      do.call(
        vegan::adonis2,
        c(
          list(
            formula      = perm_formula,
            data         = perm_data,
            permutations = perms,
            by           = adonis_args[["by"]] %||% "margin"
          ),
          adonis_args[setdiff(names(adonis_args), "by")]
        )
      ),
      error = function(e) {
        stop(
          "vegan::adonis2() failed:\n  ", conditionMessage(e), "\n",
          "  Check that 'rhs' is a valid column in 'metadata' and that the ",
          "grouping variable has at least 2 levels.",
          call. = FALSE
        )
      }
    )

    # # Extract key scalars for the grouping term (first data row)
    # adonis_df  <- as.data.frame(adonis_out)
    # grp_row    <- adonis_df[rhs, , drop = FALSE]

    # Extract key scalars for the grouping term (first data row)
    adonis_df  <- as.data.frame(adonis_out)
    
    # The grouping variable is always the first row in the adonis2 output table.
    # Extracting by index (1) completely avoids row name mismatches caused by backticks or spaces.
    grp_row    <- adonis_df[1, , drop = FALSE]

    permanova_result <- list(
      table   = adonis_df,
      R2      = grp_row[["R2"]],
      F_stat  = grp_row[["F"]],
      p_value = grp_row[["Pr(>F)"]],
      rhs     = rhs,
      perms   = perms
    )
  }

  # ── 11. Post-hoc via OmicFlow::pairwise_adonis ───────────────────────────────
  pairwise_result <- NULL

  if (!is.null(groups)) {
    if (!is.character(groups) || length(groups) != 1L) {
      stop("'groups' must be a single character string.", call. = FALSE)
    }
    if (!groups %in% colnames(metadata)) {
      stop("Column '", groups, "' specified in 'groups' not found in 'metadata'.", call. = FALSE)
    }

    if (!requireNamespace("OmicFlow", quietly = TRUE)) {
      warning(
        "'OmicFlow' is required for post-hoc analysis but is not installed.\n",
        "  Skipping pairwise PERMANOVA.",
        call. = FALSE
      )
    } else {
      # Extract args specific to pairwise_adonis
      pa_args <- .split_dots_pairwise_adonis(dots)

      pa_out <- tryCatch({
        do.call(
          OmicFlow::pairwise_adonis,
          c(
            list(
              x               = dist_matrix,
              # pairwise_adonis strictly requires a character vector here:
              groups          = as.character(metadata[[groups]]),
              # Pass the full data.frame to metadata to satisfy perm_design requirements cleanly
              metadata        = metadata,
              p.adjust.method = p_adjust,
              perm            = perms
            ),
            pa_args[setdiff(names(pa_args), "metadata")] # Prevent duplicate arg errors
          )
        )
      }, error = function(e) {
        warning("OmicFlow::pairwise_adonis() failed:\n  ", conditionMessage(e), call. = FALSE)
        NULL
      })

      pairwise_result <- pa_out
    }
  }

  # ── 12. Assemble result ───────────────────────────────────────────────────────
  result <- list(
    scores             = scores,
    variance_explained = var_pct,
    eigenvalues        = eig_values,
    dist_matrix        = dist_matrix,
    method             = method_used,
    correction         = correction,
    n_obs              = n_obs,
    n_axes             = n_axes,
    pcoa_raw           = pcoa_raw,
    vegdist_raw        = if (is_precomputed_dist) NULL else vegdist_raw,
    permanova          = permanova_result,
    pairwise_adonis    = pairwise_result,   # New element added here
    call               = .call,
    p_adjust           = p_adjust           # Saved for summary output
  )
  class(result) <- c("run_pcoa", "list")
  result
}

# ── S3 print ──────────────────────────────────────────────────────────────────

#' @export
print.run_pcoa <- function(x, n_axes = 5L, ...) {
  cat("Principal Coordinate Analysis (PCoA)\n")
  cat(strrep("-", 40L), "\n")
  cat("Dissimilarity method :", x$method, "\n")
  cat("Correction           :", x$correction, "\n")
  cat("Observations         :", x$n_obs, "\n")
  cat("Axes retained        :", x$n_axes, "\n\n")

  show_n <- min(n_axes, x$n_axes)
  cat("Variance explained (first", show_n, "axes):\n")
  ve <- x$variance_explained[seq_len(show_n)]
  cat(sprintf("  %s: %.2f%%\n", names(ve), ve), sep = "")

  if (!is.null(x$permanova)) {
    cat("\nPERMANOVA (", x$permanova$rhs, ", ",
        x$permanova$perms, " permutations):\n", sep = "")
    cat(sprintf("  R\u00b2  = %s\n", .fmt_r2(x$permanova$R2)))
    cat(sprintf("  F    = %.4f\n", x$permanova$F_stat))
    cat(sprintf("  p    = %s\n",   .fmt_pval(x$permanova$p_value)))
  }
  cat("\n")
  invisible(x)
}


# ── S3 summary ────────────────────────────────────────────────────────────────

#' @export
summary.run_pcoa <- function(object, ...) {
  cat("Principal Coordinate Analysis Summary\n")
  cat(strrep("=", 45L), "\n")
  cat("Call:\n  ", deparse(object$call), "\n\n")
  cat("Dissimilarity method :", object$method, "\n")
  cat("Correction           :", object$correction, "\n")
  cat("Observations (n)     :", object$n_obs, "\n")
  cat("Total axes           :", object$n_axes, "\n\n")

  cat("Eigenvalues and variance explained:\n")
  df <- data.frame(
    Axis           = names(object$variance_explained),
    Eigenvalue     = round(object$eigenvalues[seq_len(object$n_axes)], 6L),
    Variance_pct   = round(object$variance_explained, 3L),
    Cumulative_pct = round(cumsum(object$variance_explained), 3L),
    stringsAsFactors = FALSE,
    check.names    = FALSE
  )
  print(df, row.names = FALSE)

  if (!is.null(object$permanova)) {
    cat("\nPERMANOVA results (", object$permanova$rhs, ", ",
        object$permanova$perms, " permutations):\n", sep = "")
    print(object$permanova$table)
  }
  invisible(object)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

.split_dots_vegdist <- function(dots) {
  dots[intersect(names(dots), names(formals(vegan::vegdist)))]
}

.split_dots_pcoa <- function(dots) {
  dots[intersect(names(dots), names(formals(ape::pcoa)))]
}

.split_dots_adonis2 <- function(dots) {
  dots[intersect(names(dots), names(formals(vegan::adonis2)))]
}

.extract_eigenvalues <- function(pcoa_obj, correction) {
  vals <- pcoa_obj$values
  if (correction != "none" && !is.null(vals$Corr_eig)) return(vals$Corr_eig)
  vals$Eigenvalues
}

# Null-coalescing operator (defined locally to avoid rlang dependency)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Format R2: expand decimal places until a non-1 representation is found,
# cap at 5 dp then switch to scientific notation with 2 dp.
.fmt_r2 <- function(r2) {
  if (is.na(r2)) return("NA")
  for (dp in 2:5) {
    s <- formatC(r2, format = "f", digits = dp)
    if (as.numeric(s) != 1) return(s)
  }
  formatC(r2, format = "e", digits = 2)
}

# Format p-value: 3 dp maximum; values < 0.001 shown as "<0.001".
.fmt_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(round(p, 3), format = "f", digits = 3)
}

.split_dots_pairwise_adonis <- function(dots) {
  if (!requireNamespace("OmicFlow", quietly = TRUE)) return(list())
  dots[intersect(names(dots), names(formals(OmicFlow::pairwise_adonis)))]
}