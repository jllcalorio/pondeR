#' Impute Missing Values in Metabolomics Data
#'
#' @description
#' Replaces missing values (NAs) using either a simple deterministic method
#' (fraction of minimum), a statistical imputation method (quantile regression
#' imputation of left-censored data, QRILC), or any method available in
#' \code{\link[mice]{mice}} (Multiple Imputation by Chained Equations).
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features
#'   in columns. Missing values should be represented as `NA`.
#' @param method Numeric or character. Controls the imputation strategy:
#'   \itemize{
#'     \item **Numeric** — fraction of the smallest observed value per feature
#'       used for deterministic imputation (e.g., `0.2` = 1/5th of the
#'       smallest value). By default (`positive_only = TRUE`), the minimum is
#'       computed from strictly positive values only, which is the standard
#'       MNAR convention in metabolomics. Set `positive_only = FALSE` to use
#'       the global minimum including zeros and negative values. Recommended
#'       range: 0.01–0.5. Values > 1 trigger a warning.
#'     \item **`"quantileregression"`** — uses the QRILC algorithm from the
#'       \pkg{imputeLCMD} package. Assumes Missing Not At Random (MNAR) data.
#'       Controlled further by `tune_sigma`.
#'     \item **Any `mice` method string** — e.g., `"pmm"`, `"norm"`,
#'       `"norm.boot"`, `"norm.nob"`, `"mean"`, `"2l.norm"`, `"rf"`,
#'       `"cart"`, etc. Delegates to \code{\link[mice]{mice}} from the
#'       \pkg{mice} package. Controlled further by `m`, `maxit`, `seed`,
#'       and `...`.
#'   }
#'   Default: `0.2` (1/5th of minimum, deterministic).
#' @param positive_only Logical. Only used when `method` is numeric. If `TRUE`
#'   (default), the per-feature minimum is computed from strictly positive
#'   values only (zeros and negatives are excluded), which is standard practice
#'   for MNAR imputation in metabolomics where missing values represent
#'   abundances below the detection limit. If `FALSE`, the global minimum
#'   across all observed (non-NA) values is used, which may be appropriate
#'   when features can legitimately be zero or negative (e.g., log-transformed
#'   or mean-centred data). Has no effect when `method` is
#'   `"quantileregression"` or a \pkg{mice} method string.
#' @param tune_sigma Numeric. Only used when `method = "quantileregression"`.
#'   Controls the standard deviation of the left-censored missing value
#'   distribution. Smaller values (e.g., `0.5`) assume a tighter distribution;
#'   larger values (e.g., `2`) allow more spread. Default: `1`.
#' @param m Integer. Only used when `method` routes to \pkg{mice}. Number of
#'   multiple imputations. The first completed dataset is used as the imputed
#'   result. Default: `5`.
#' @param maxit Integer. Only used when `method` routes to \pkg{mice}. Number
#'   of iterations for the chained equations algorithm. Default: `5`.
#' @param seed Integer or `NA`. Only used when `method` routes to \pkg{mice}.
#'   Random seed passed to \code{\link[mice]{mice}} for reproducibility.
#'   Default: `NA` (no seed set).
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#' @param ... Additional arguments passed to \code{\link[mice]{mice}} when
#'   `method` routes to the \pkg{mice} engine. Has no effect for numeric or
#'   `"quantileregression"` methods. A notable argument is `predictorMatrix`:
#'   by default, `run_mvimpute` builds one automatically via
#'   \code{\link[mice]{quickpred}} to avoid failures caused by constant or
#'   collinear features, which are common in wide metabolomics/metagenomics
#'   data. Supply your own `predictorMatrix` via `...` to override this.
#'
#' @return A list with the following elements:
#'   \item{data}{Matrix or data frame with imputed values (same class as input `x`).}
#'   \item{n_missing_before}{Integer. Number of missing values before imputation.}
#'   \item{n_missing_after}{Integer. Number of missing values after imputation
#'     (should be 0).}
#'   \item{imputed_summary}{Data frame with per-feature imputation statistics.}
#'   \item{parameters}{List of parameters used, including method, positive_only,
#'     tune_sigma (QRILC only), and m/maxit/seed (mice only).}
#'   \item{method_used}{Character string describing the imputation method applied.}
#'
#' @details
#' **Missing Value Imputation Strategies:**
#'
#' 1. **Deterministic (fraction of minimum):** When `method` is numeric, each
#'    missing value in a feature is replaced by a fraction of the per-feature
#'    minimum, multiplied by the specified fraction. The minimum is determined
#'    by `positive_only`:
#'    \itemize{
#'      \item `positive_only = TRUE` (default): uses the smallest strictly
#'        positive observed value per feature, excluding zeros and negatives.
#'        This is standard in metabolomics for values below the detection limit
#'        (MNAR — Missing Not At Random).
#'      \item `positive_only = FALSE`: uses the global minimum across all
#'        observed (non-NA) values, including zeros and negatives. Suitable for
#'        data that has been log-transformed or mean-centred.
#'    }
#'    If no valid minimum can be found (e.g., a feature is entirely NA or, when
#'    `positive_only = TRUE`, has no positive values), a fallback of `1e-9` is
#'    used.
#'
#' 2. **Quantile Regression (QRILC):** When `method = "quantileregression"`,
#'    uses the Quantile Regression Imputation of Left-Censored data algorithm
#'    from the \pkg{imputeLCMD} package. Models missing values from a
#'    left-censored distribution, controlled by `tune_sigma`. More
#'    statistically principled than simple fraction-of-minimum for MNAR data.
#'    Requires \pkg{imputeLCMD}: `BiocManager::install("imputeLCMD")`.
#'
#' 3. **mice (Multiple Imputation by Chained Equations):** When `method` is
#'    a recognised \pkg{mice} method string, \code{\link[mice]{mice}} is
#'    called internally on the feature matrix. The first of the `m` completed
#'    datasets is extracted via `mice::complete()` and returned. Supported
#'    `mice` methods include (non-exhaustive): `"pmm"` (predictive mean
#'    matching), `"norm"`, `"norm.boot"`, `"norm.nob"`, `"mean"`, `"rf"`
#'    (random forest), `"cart"`, `"2l.norm"`, `"midastouch"`, `"sample"`.
#'    See \code{\link[mice]{mice}} for the full list. Requires \pkg{mice}:
#'    `install.packages("mice")`. By default, `run_mvimpute` builds a
#'    predictor matrix automatically via \code{\link[mice]{quickpred}} to
#'    guard against the common failure where `mice` removes all predictors due
#'    to constant or collinear columns (frequent in wide, sparse
#'    metabolomics/metagenomics data). Supply a custom `predictorMatrix` via
#'    `...` to override. Additional arguments for \pkg{mice} can also be
#'    passed via `...`.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
#' Imputation by Chained Equations in R. *Journal of Statistical Software*,
#' 45(3), 1–67. \doi{10.18637/jss.v045.i03}
#'
#' Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
#' Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics
#' Data. *Scientific Reports*, 8(1), 663.
#' \doi{10.1038/s41598-017-19120-0}
#'
#' Lazar, C., Gatto, L., Ferro, M., Bruley, C., & Burger, T. (2016).
#' Accounting for the Multiple Natures of Missing Values in Label-Free
#' Quantitative Proteomics Data Sets to Compare Imputation Strategies.
#' *Journal of Proteome Research*, 15(4), 1116–1125.
#' \doi{10.1021/acs.jproteome.5b00981}
#'
#' @seealso \code{\link[mice]{mice}}, \code{\link[imputeLCMD]{impute.QRILC}}
#'
#' @importFrom mice mice complete
#' @importFrom imputeLCMD impute.QRILC
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data with missing values
#' set.seed(123)
#' x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
#' x[sample(length(x), 200)] <- NA
#' colnames(x) <- paste0("Feature", 1:50)
#'
#' # 1. Deterministic imputation (1/5th of minimum positive value per feature)
#' result1 <- run_mvimpute(x, method = 0.2)
#'
#' # 2. Deterministic imputation using global minimum (including zeros/negatives)
#' result2 <- run_mvimpute(x, method = 0.2, positive_only = FALSE)
#'
#' # 3. Quantile regression imputation (QRILC)
#' result3 <- run_mvimpute(x, method = "quantileregression", tune_sigma = 1)
#'
#' # 4. Predictive mean matching via mice
#' result4 <- run_mvimpute(x, method = "pmm", m = 5, maxit = 5, seed = 42)
#'
#' # 5. Random forest imputation via mice
#' result5 <- run_mvimpute(x, method = "rf", m = 5, maxit = 5, seed = 42)
#'
#' # 6. Passing additional mice arguments via ...
#' result6 <- run_mvimpute(x, method = "norm.boot", m = 10, printFlag = FALSE)
#' }
run_mvimpute <- function(
    x,
    method       = 0.2,
    positive_only = TRUE,
    tune_sigma   = 1,
    m            = 5,
    maxit        = 5,
    seed         = NA,
    verbose      = TRUE,
    ...
) {

  msg <- function(...) if (verbose) message(...)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  if (!is.logical(positive_only) || length(positive_only) != 1L) {
    stop("'positive_only' must be a single logical value (TRUE or FALSE).")
  }

  was_matrix <- is.matrix(x)
  x_matrix   <- as.matrix(x)

  n_missing_before <- sum(is.na(x_matrix))

  if (n_missing_before == 0L) {
    msg("No missing values detected. Returning original data.")
    return(list(
      data             = if (was_matrix) x_matrix else as.data.frame(x_matrix),
      n_missing_before = 0L,
      n_missing_after  = 0L,
      imputed_summary  = data.frame(
        Feature          = colnames(x_matrix),
        N_Missing        = 0L,
        Imputation_Value = NA_real_
      ),
      parameters  = list(method = method, positive_only = positive_only,
                         tune_sigma = tune_sigma, m = m, maxit = maxit,
                         seed = seed),
      method_used = "None (no missing values)"
    ))
  }

  msg(sprintf("Starting missing value imputation (%d missing values)...",
              n_missing_before))

  # ---------------------------------------------------------------------------
  # Dispatch: numeric -> deterministic | "quantileregression" -> QRILC
  #           anything else -> mice
  # ---------------------------------------------------------------------------

  if (is.numeric(method)) {

    # -- Deterministic (fraction of minimum) ----------------------------------
    if (length(method) != 1L) {
      stop("'method' when numeric must be a single scalar value.")
    }
    if (method <= 0) {
      stop("'method' when numeric must be a positive value.")
    }
    if (method > 1) {
      warning(
        "'method' > 1 means imputed values will exceed the observed minimum. ",
        "This is unusual for MNAR imputation. Consider using a value < 1."
      )
    }

    if (positive_only) {
      msg(sprintf(
        "Using deterministic imputation: %.4f x minimum positive value per feature.",
        method
      ))
      get_min <- function(col) {
        v <- col[!is.na(col) & col > 0]
        if (length(v) == 0L) 1e-9 else min(v, na.rm = TRUE)
      }
    } else {
      msg(sprintf(
        "Using deterministic imputation: %.4f x minimum observed value per feature (including zeros/negatives).",
        method
      ))
      get_min <- function(col) {
        v <- col[!is.na(col)]
        if (length(v) == 0L) 1e-9 else min(v, na.rm = TRUE)
      }
    }

    na_indices    <- is.na(x_matrix)
    min_vals      <- apply(x_matrix, 2L, get_min)
    imputation_vals <- min_vals * method
    x_matrix[na_indices] <- imputation_vals[col(x_matrix)][na_indices]

    method_description <- sprintf(
      "Deterministic: %.4f x min %s",
      method,
      if (positive_only) "positive" else "observed (all)"
    )

    imputed_summary <- data.frame(
      Feature          = colnames(x_matrix),
      N_Missing        = colSums(na_indices),
      Min_Value        = min_vals,
      Imputation_Value = imputation_vals
    )

    params_out <- list(method = method, positive_only = positive_only,
                       tune_sigma = NA, m = NA, maxit = NA, seed = NA)

  } else if (is.character(method)) {

    method_lc <- tolower(trimws(method))

    if (method_lc == "quantileregression") {

      # -- QRILC --------------------------------------------------------------
      if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
        stop(
          "Package 'imputeLCMD' is required for quantile regression imputation. ",
          "Install it with: BiocManager::install('imputeLCMD')"
        )
      }
      if (!is.numeric(tune_sigma) || length(tune_sigma) != 1L || tune_sigma <= 0) {
        stop("'tune_sigma' must be a single positive numeric value.")
      }

      msg(sprintf(
        "Using quantile regression imputation (QRILC) with tune_sigma = %.2f.",
        tune_sigma
      ))

      na_indices <- is.na(x_matrix)

      tryCatch({
        x_imputed <- imputeLCMD::impute.QRILC(t(x_matrix), tune.sigma = tune_sigma)
        x_matrix  <- t(x_imputed[[1]])
      }, error = function(e) {
        stop(
          "Quantile regression imputation failed: ", e$message,
          "\nConsider using a deterministic method instead."
        )
      })

      method_description <- sprintf("QRILC (tune_sigma = %.2f)", tune_sigma)

      imputed_summary <- data.frame(
        Feature   = colnames(x_matrix),
        N_Missing = colSums(na_indices),
        Method    = "QRILC"
      )

      params_out <- list(method = method, positive_only = NA,
                         tune_sigma = tune_sigma, m = NA, maxit = NA, seed = NA)

    } else {

      # -- mice ---------------------------------------------------------------
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop(
          "Package 'mice' is required for method = '", method, "'. ",
          "Install it with: install.packages('mice')"
        )
      }

      if (!is.numeric(m) || length(m) != 1L || m < 1L) {
        stop("'m' must be a single positive integer.")
      }
      if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1L) {
        stop("'maxit' must be a single positive integer.")
      }
      if (!is.na(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
        stop("'seed' must be a single integer or NA.")
      }

      msg(sprintf(
        "Using mice imputation (method = '%s', m = %d, maxit = %d, seed = %s).",
        method, as.integer(m), as.integer(maxit),
        if (is.na(seed)) "NA" else as.character(as.integer(seed))
      ))

      na_indices <- is.na(x_matrix)

      dots <- list(...)

      if (is.null(dots$predictorMatrix)) {
        pred_matrix <- tryCatch(
          mice::quickpred(as.data.frame(x_matrix)),
          error = function(e) {
            stop(
              "mice::quickpred() failed to build a predictor matrix: ", e$message,
              "\nYour data may have too few samples or too many constant/collinear ",
              "features for mice-based imputation. ",
              "Consider using method = 0.2 or method = 'quantileregression' instead."
            )
          }
        )
        dots$predictorMatrix <- pred_matrix
        msg("Predictor matrix built automatically via mice::quickpred().")
      }

      tryCatch({
        mice_obj <- do.call(
          mice::mice,
          c(
            list(
              data      = as.data.frame(x_matrix),
              method    = method,
              m         = as.integer(m),
              maxit     = as.integer(maxit),
              seed      = if (is.na(seed)) NA_integer_ else as.integer(seed),
              printFlag = FALSE
            ),
            dots
          )
        )
        x_matrix <- as.matrix(mice::complete(mice_obj, action = 1L))
      }, error = function(e) {
        stop(
          "mice imputation failed (method = '", method, "'): ", e$message,
          "\nIf the error mentions constant or collinear variables, your data may ",
          "be too sparse or wide for mice-based imputation. ",
          "Consider using method = 0.2 or method = 'quantileregression' instead."
        )
      })

      method_description <- sprintf(
        "mice: %s (m = %d, maxit = %d)", method,
        as.integer(m), as.integer(maxit)
      )

      imputed_summary <- data.frame(
        Feature   = colnames(x_matrix),
        N_Missing = colSums(na_indices),
        Method    = method
      )

      params_out <- list(method = method, positive_only = NA,
                         tune_sigma = NA, m = as.integer(m),
                         maxit = as.integer(maxit),
                         seed = if (is.na(seed)) NA_integer_ else as.integer(seed))
    }

  } else {
    stop(
      "'method' must be numeric (fraction of minimum), ",
      "'quantileregression', or a mice method string (e.g., 'pmm', 'rf')."
    )
  }

  # ---------------------------------------------------------------------------
  # Finalise and return
  # ---------------------------------------------------------------------------
  n_missing_after <- sum(is.na(x_matrix))

  msg(sprintf("Imputation complete. Missing values: %d -> %d.",
              n_missing_before, n_missing_after))

  list(
    data             = if (was_matrix) x_matrix else as.data.frame(x_matrix),
    n_missing_before = n_missing_before,
    n_missing_after  = n_missing_after,
    imputed_summary  = imputed_summary,
    parameters       = params_out,
    method_used      = method_description
  )
}

# S3 print method for cleaner output
#' @export
print.run_mvimpute <- function(x, ...) {
  cat("=== Missing Value Imputation Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  cat(sprintf("Missing values before: %d\n", x$n_missing_before))
  cat(sprintf("Missing values after:  %d\n", x$n_missing_after))
  invisible(x)
}