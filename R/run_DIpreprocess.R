#' Direct-Injection Metabolomics Preprocessing Pipeline
#'
#' @description
#' Runs a complete, sequential preprocessing pipeline for metabolomics data by
#' orchestrating the individual \code{run_*} functions of \pkg{pondeR}. Steps
#' execute in a fixed, reproducible order regardless of argument arrangement.
#'
#' @param x Matrix or data frame. Numeric feature data with \strong{samples in rows}
#'   and \strong{features in columns}. All values should be non-negative raw intensities;
#'   zeros are treated as missing internally.
#' @param metadata Data frame. Sample metadata with \code{nrow(metadata) == nrow(x)}.
#'   At minimum, the group column (default \code{"Group"}) must be present.
#' @param sample_id_col Character. Column in \code{metadata} containing unique
#'   sample identifiers that match \code{rownames(x)}. Default: \code{"Sample"}.
#' @param group_col Character. Column in \code{metadata} containing group labels.
#'   Default: \code{"Group"}.
#' @param qc_types Character vector. Group labels that identify QC samples. 
#'   Default: \code{c("QC", "SQC", "EQC")}.
#' @param batch_col Character. Column in \code{metadata} containing batch numbers.
#'   Required for drift correction. Default: \code{"Batch"}.
#' @param injection_col Character. Column in \code{metadata} containing injection
#'   order (integer). Required for drift correction. Default: \code{"InjectionSequence"}.
#' @param norm_factor_col Character. Column in \code{metadata} containing
#'   external normalization factors. Default: \code{"Normalization"}.
#' @param subject_id_col Character. Column in \code{metadata} containing
#'   subject identifiers for technical-replicate merging. Default: \code{"SubjectID"}.
#' @param outliers Character vector or \code{NULL}. Row names of \code{x} 
#'   to remove before processing. Default: \code{NULL}.
#' @param missing_threshold Numeric \[0, 1\]. Missing value threshold. Default: \code{0.2}.
#' @param missing_by_group Logical. If \code{TRUE}, assess per group. Default: \code{TRUE}.
#' @param missing_include_qc Logical. If \code{TRUE}, include QC samples. Default: \code{FALSE}.
#' @param impute_fraction Numeric > 0. Fraction of the smallest positive value 
#'   per feature used for imputation. Default: \code{0.2}.
#' @param positive_only Logical. If \code{TRUE}, imputation only considers positive values. Default: \code{TRUE}.
#' @param correct_drift Logical. If \code{TRUE}, apply QCRSC correction. Default: \code{TRUE}.
#' @param remove_uncorrected Logical. If \code{TRUE}, remove features QCRSC 
#'   could not correct. Default: \code{FALSE}.
#' @param spline_smooth Numeric \[0, 1\]. Smoothing parameter. Default: \code{0}.
#' @param spline_spar_limit Numeric vector. Lower/upper bounds for spline. Default: \code{c(-1.5, 1.5)}.
#' @param correct_on_log Logical. If \code{TRUE}, fit on log-transformed data. Default: \code{TRUE}.
#' @param min_qc_per_batch Integer. Minimum QC samples per batch. Default: \code{5}.
#' @param normalize_method Character. Normalization method.
#' \itemize{
#'   \item \code{"sum"}: Normalizes by sum.
#'   \item \code{"median"}: Normalizes by median.
#'   \item \code{"specific_factor"}: Uses external factors in \code{norm_factor_col}.
#'   \item \code{"pqn_global"}: PQN using a global reference.
#'   \item \code{"pqn_reference"}: PQN using \code{normalize_ref_sample}.
#'   \item \code{"pqn_group"}: PQN grouping.
#'   \item \code{"quantile"}: Quantile normalization.
#'   \item \code{"none"}: No normalization.
#' }
#' Default: \code{"sum"}.
#' @param normalize_ref_sample Character. Reference sample name for PQN. Default: \code{NULL}.
#' @param normalize_qc_method Character. QC normalization when using factors.
#'   One of \code{"mean"}, \code{"median"}, \code{"none"}. Default: \code{"median"}.
#' @param transform_method Character. Transformation method. 
#' \itemize{
#'   \item \code{"log2"}, \code{"log10"}, \code{"sqrt"}, \code{"cbrt"}, 
#'   \code{"clr"}, \code{"vsn"}, \code{"glog"}, or \code{"none"}.
#' }
#' Default: \code{"log2"}.
#' @param vsn_cores Integer or \code{"max"}. Cores for VSN. Default: \code{1}.
#' @param scale_nonpls Character. Scaling for NONPLS branch (\code{"auto"}, \code{"pareto"}, \code{"mean"}).
#' @param scale_pls Character. Scaling for PLS branch (\code{"pareto"}, \code{"auto"}, \code{"mean"}).
#' @param rsd_threshold Numeric \[0, 1\]. QC-RSD threshold. Default: \code{0.3}.
#' @param rsd_qc_type Character. QC type for RSD (\code{"EQC"}, \code{"SQC"}, \code{"QC"}).
#' @param variance_percentile Numeric \[0, 100\]. Percentile to filter. Default: \code{10}.
#' @param scale_filter_ref Character. Harmonisation strategy.
#' \itemize{
#'   \item \code{"auto"}: Keep intersection of both branches.
#'   \item \code{"NONPLS"}: Apply NONPLS filters to both.
#'   \item \code{"PLS"}: Apply PLS filters to both.
#' }
#' @param merge_replicates Logical. Average technical replicates. Default: \code{FALSE}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return A named list of class \code{run_DIpreprocess} containing:
#' \itemize{
#'   \item \code{metadata}: Processed sample metadata.
#'   \item \code{data_raw}: Original feature matrix.
#'   \item \code{data_nonpls}: Final NONPLS-scaled, filtered data.
#'   \item \code{data_pls}: Final PLS-scaled, filtered data.
#'   \item \code{features_final}: Character vector of retained features.
#'   \item \code{dimensions}: Sample and feature counts at each step.
#'   \item \code{elapsed_seconds}: Execution time.
#' }
#'
#' @details
#' \strong{Preprocessing Workflow (Fixed Order):}
#' \enumerate{
#'   \item Input validation
#'   \item Outlier removal
#'   \item Missing-value filtering (\code{\link{run_filtermissing}})
#'   \item Missing-value imputation (\code{\link{run_mvimpute}})
#'   \item Signal drift correction (\code{\link{run_driftBatchCorrect}})
#'   \item Normalization (\code{\link{run_normalize}})
#'   \item Transformation (\code{\link{run_transform}})
#'   \item Scaling (\code{\link{run_scale}})
#'   \item Quality filtering (\code{\link{run_filterRSD}}, \code{\link{run_filtervariance}})
#'   \item Common-feature harmonisation
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @seealso
#' \code{\link{run_filtermissing}}, \code{\link{run_mvimpute}},
#' \code{\link{run_driftBatchCorrect}}, \code{\link{run_normalize}}
#'
#' @references
#' Kirwan, J.A., et al. (2013). \emph{Analytical and Bioanalytical Chemistry}, 405, 5147–5157.
#'
#' @importFrom matrixStats colMaxs
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n_s <- 80; n_f <- 100
#' x <- matrix(abs(rnorm(n_s * n_f, 500, 150)), nrow = n_s, ncol = n_f)
#' colnames(x) <- paste0("Feature", seq_len(n_f))
#' x[sample(length(x), 400)] <- 0
#'
#' meta <- data.frame(
#'   Sample          = paste0("S", seq_len(n_s)),
#'   Group           = rep(c("Control", "Treatment", "QC"), c(30, 30, 20)),
#'   Batch           = rep(1:2, each = 40),
#'   InjectionSequence = seq_len(n_s),
#'   SubjectID       = c(paste0("BIO", seq_len(60)), rep(NA, 20)),
#'   stringsAsFactors = FALSE
#' )
#' rownames(x) <- meta$Sample
#'
#' result <- run_DIpreprocess(
#'   x              = x,
#'   metadata       = meta,
#'   normalize_method = "pqn_group",
#'   transform_method = "log2",
#'   scale_nonpls   = "auto",
#'   scale_pls      = "pareto",
#'   correct_drift  = FALSE
#' )
#'
#' dim(result$data_nonpls)
#' result$dimensions
#' }
run_DIpreprocess <- function(
    x,
    metadata,

    # — Column names ——————————————————————————————————————————————————————————
    sample_id_col    = "Sample",
    group_col        = "Group",
    qc_types         = c("QC", "SQC", "EQC"),
    batch_col        = "Batch",
    injection_col    = "InjectionSequence",
    norm_factor_col  = "Normalization",
    subject_id_col   = "SubjectID",

    # — Outlier removal ———————————————————————————————————————————————————————
    outliers         = NULL,

    # — Missing-value filter ——————————————————————————————————————————————————
    missing_threshold  = 0.2,
    missing_by_group   = TRUE,
    missing_include_qc = FALSE,

    # — Missing-value imputation ——————————————————————————————————————————————
    impute_fraction  = 0.2,
    positive_only    = TRUE,

    # — Drift / batch correction ——————————————————————————————————————————————
    correct_drift        = TRUE,
    remove_uncorrected   = FALSE,
    spline_smooth        = 0,
    spline_spar_limit    = c(-1.5, 1.5),
    correct_on_log       = TRUE,
    min_qc_per_batch     = 5L,

    # — Normalization —————————————————————————————————————————————————————————
    normalize_method     = "sum",
    normalize_ref_sample = NULL,
    normalize_qc_method  = "median",

    # — Transformation ————————————————————————————————————————————————————————
    transform_method = "log2",
    vsn_cores        = 1L,

    # — Scaling ———————————————————————————————————————————————————————————————
    scale_nonpls     = "auto",
    scale_pls        = "pareto",

    # — Quality filtering —————————————————————————————————————————————————————
    rsd_threshold      = 0.3,
    rsd_qc_type        = "EQC",
    variance_percentile = 10,
    scale_filter_ref   = "auto",

    # — Replicate merging —————————————————————————————————————————————————————
    merge_replicates = FALSE,

    verbose          = TRUE
) {

  t_start <- proc.time()[["elapsed"]]
  msg <- function(...) if (verbose) message(...)

  # ---------------------------------------------------------------------------
  # Internal helpers
  # ---------------------------------------------------------------------------

  .record_dim <- function(dims_df, label, data) {
    rbind(dims_df, data.frame(
      Step     = label,
      Samples  = nrow(data),
      Features = ncol(data),
      stringsAsFactors = FALSE
    ))
  }

  .safe_impute <- function(mat, fraction) {
    # Deterministic: fraction * column minimum (positive values only)
    pos_min <- function(col) {
      v <- col[!is.na(col) & col > 0]
      if (length(v) == 0L) 1e-9 else min(v)
    }
    min_vals <- apply(mat, 2L, pos_min) * fraction
    na_idx   <- is.na(mat)
    mat[na_idx] <- min_vals[col(mat)][na_idx]
    mat
  }

  # ---------------------------------------------------------------------------
  # Result scaffold
  # ---------------------------------------------------------------------------

  out <- list(
    metadata         = NULL,
    data_raw         = NULL,
    data_missing_filtered = NULL,
    data_imputed     = NULL,
    data_corrected   = NULL,
    data_normalized  = NULL,
    data_transformed = NULL,
    data_nonpls      = NULL,
    data_pls         = NULL,
    data_nonpls_merged = NULL,
    data_pls_merged    = NULL,
    metadata_merged    = NULL,
    features_final   = character(0L),
    uncorrected_features = character(0L),
    dimensions       = data.frame(
      Step = character(), Samples = integer(), Features = integer(),
      stringsAsFactors = FALSE
    ),
    parameters       = as.list(match.call()[-1L]),
    elapsed_seconds  = NA_real_,
    error            = NULL,
    other_data = NULL # REMOVE AFTER CHECKING
  )

  tryCatch({

    # =========================================================================
    # 1. INPUT VALIDATION
    # =========================================================================

    msg("Step 1/10: Validating inputs...")

    if (!is.matrix(x) && !is.data.frame(x))
      stop("'x' must be a matrix or data frame.")
    if (!is.data.frame(metadata))
      stop("'metadata' must be a data frame.")
    if (nrow(x) != nrow(metadata))
      stop("nrow(x) != nrow(metadata).")
    if (!is.numeric(missing_threshold) || length(missing_threshold) != 1L ||
        missing_threshold < 0 || missing_threshold > 1)
      stop("'missing_threshold' must be a single numeric value in [0, 1].")
    if (!is.numeric(impute_fraction) || length(impute_fraction) != 1L ||
        impute_fraction <= 0)
      stop("'impute_fraction' must be a single positive numeric value.")
    if (!scale_filter_ref %in% c("auto", "NONPLS", "PLS"))
      stop("'scale_filter_ref' must be one of: 'auto', 'NONPLS', 'PLS'.")
    if (!rsd_qc_type %in% c("QC", "SQC", "EQC"))
      stop("'rsd_qc_type' must be one of: 'QC', 'SQC', 'EQC'.")
    if (!is.numeric(variance_percentile) || length(variance_percentile) != 1L ||
        variance_percentile < 0 || variance_percentile > 100)
      stop("'variance_percentile' must be a numeric value in [0, 100].")
    if (!group_col %in% colnames(metadata))
      stop(sprintf("group_col '%s' not found in metadata.", group_col))

    was_matrix <- is.matrix(x)
    df <- as.data.frame(as.matrix(x))

    # Synchronise row names: prefer sample_id_col from metadata
    if (sample_id_col %in% colnames(metadata)) {
      rownames(df)       <- as.character(metadata[[sample_id_col]])
      rownames(metadata) <- as.character(metadata[[sample_id_col]])
    } else if (!is.null(rownames(x))) {
      rownames(metadata) <- rownames(x)
    }

    # Ensure all feature columns are numeric
    df[] <- lapply(df, function(col) suppressWarnings(as.numeric(as.character(col))))

    out$dimensions <- .record_dim(out$dimensions, "Original", df)

    # =========================================================================
    # 2. OUTLIER REMOVAL
    # =========================================================================

    msg("Step 2/10: Checking for outliers...")

    if (!is.null(outliers) && length(outliers) > 0L) {
      present <- intersect(outliers, rownames(df))
      missing_o <- setdiff(outliers, rownames(df))

      if (length(present) > 0L) {
        df       <- df[!rownames(df) %in% present, , drop = FALSE]
        metadata <- metadata[!rownames(metadata) %in% present, , drop = FALSE]
        msg(sprintf("  Removed %d outlier(s): %s",
                    length(present), paste(present, collapse = ", ")))
      }
      if (length(missing_o) > 0L)
        warning("Outliers not found in data: ", paste(missing_o, collapse = ", "),
                call. = FALSE)

      out$dimensions <- .record_dim(out$dimensions, "After outlier removal", df)
    }

    # =========================================================================
    # METADATA SNAPSHOT (after outlier removal)
    # =========================================================================

    # Helper: safely pull a metadata column or return NA vector
    .meta_col <- function(col) {
      if (col %in% colnames(metadata)) metadata[[col]]
      else rep(NA, nrow(metadata))
    }

    meta_snap <- data.frame(
      Sample            = rownames(metadata),
      Group             = .meta_col(group_col),
      Group_            = ifelse(.meta_col(group_col) %in% qc_types,
                                 "QC", .meta_col(group_col)),
      Batch             = suppressWarnings(as.integer(.meta_col(batch_col))),
      InjectionSequence = suppressWarnings(as.numeric(.meta_col(injection_col))),
      Normalization     = suppressWarnings(as.numeric(.meta_col(norm_factor_col))),
      SubjectID         = as.character(.meta_col(subject_id_col)),
      stringsAsFactors  = FALSE,
      row.names         = rownames(metadata)
    )
    out$metadata <- meta_snap

    qc_rows     <- meta_snap$Group_ == "QC"
    non_qc_rows <- !qc_rows

    # =========================================================================
    # ZEROS → NA  (always; not a user parameter)
    # =========================================================================

    n_zeros <- sum(df == 0L, na.rm = TRUE)
    if (n_zeros > 0L) {
      df[df == 0L] <- NA
      msg(sprintf("  Converted %d zero(s) to NA (structural missingness).", n_zeros))
    }

    # Remove all-NA columns before filtering
    all_na <- colSums(is.na(df), na.rm = TRUE) == nrow(df)
    if (any(all_na)) {
      df <- df[, !all_na, drop = FALSE]
      msg(sprintf("  Removed %d all-NA feature(s).", sum(all_na)))
    }

    out$data_raw   <- df
    out$dimensions <- .record_dim(out$dimensions, "After zero → NA conversion", df)

    # =========================================================================
    # 3. MISSING-VALUE FILTERING
    # =========================================================================

    msg(sprintf("Step 3/10: Missing-value filtering (threshold = %.0f%%)...",
                missing_threshold * 100))

    filt_miss <- run_filtermissing(
      x              = df,
      metadata       = metadata,
      threshold      = missing_threshold,
      filter_by_group = missing_by_group,
      include_QC     = missing_include_qc,
      group_col      = group_col,
      qc_types       = qc_types,
      zero_as_missing = FALSE,   # already converted above
      verbose        = FALSE
    )

    df             <- filt_miss$data
    out$data_missing_filtered <- df
    out$dimensions <- .record_dim(out$dimensions, "After missing filter", df)

    msg(sprintf("  Removed %d feature(s) (%.1f%%).",
                filt_miss$n_features_removed,
                filt_miss$n_features_removed / filt_miss$n_features_before * 100))

    # =========================================================================
    # 4. MISSING-VALUE IMPUTATION
    # =========================================================================

    msg(sprintf("Step 4/10: Imputing missing values (fraction = %.4f)...",
                impute_fraction))

    imp <- run_mvimpute(
      x             = df,
      method        = impute_fraction,
      positive_only = positive_only,
      verbose       = FALSE
    )

    df           <- imp$data
    out$data_imputed <- df

    msg(sprintf("  Imputed %d missing value(s).", imp$n_missing_before))

    # =========================================================================
    # 5. SIGNAL DRIFT AND BATCH CORRECTION
    # =========================================================================
    
    msg("Step 5/10: Drift/batch correction...")

    uncorrected_features <- character(0L)

    if (correct_drift) {

      required_meta <- c(batch_col, injection_col, group_col)
      missing_meta  <- setdiff(required_meta, colnames(metadata))
      if (length(missing_meta) > 0L) {
        warning("Skipping drift correction: metadata column(s) not found: ",
                paste(missing_meta, collapse = ", "), call. = FALSE)
        correct_drift <- FALSE
      }
    }

    if (correct_drift) {

      corr_result <- run_driftBatchCorrect(
        x                  = df,
        metadata           = metadata,
        perform_correction = TRUE,
        batch_corr_only    = FALSE,
        injection_sequence = injection_col,
        batch_numbers      = batch_col,
        groups             = group_col,
        qc_label           = "QC",
        qc_types           = qc_types,
        spline_smooth_param = spline_smooth,
        min_QC             = min_qc_per_batch,
        spar_limit         = spline_spar_limit,
        log_scale          = correct_on_log,
        verbose            = FALSE
      )

      uncorrected_features <- corr_result$uncorrected_features

      if (remove_uncorrected && length(uncorrected_features) > 0L) {
        corr_df <- corr_result$data
        keep    <- setdiff(colnames(corr_df), uncorrected_features)
        df      <- corr_df[, keep, drop = FALSE]
        msg(sprintf("  Removed %d uncorrected feature(s).",
                    length(uncorrected_features)))
      } else {
        df <- corr_result$data
        if (length(uncorrected_features) > 0L)
          msg(sprintf("  %d feature(s) could not be corrected (retained).",
                      length(uncorrected_features)))
      }

      # Re-impute any NAs introduced by QCRSC.
      # Compare against the pre-correction NA count (should be 0 after step 4,
      # but diffing is safer than assuming).
      na_before_corr <- sum(is.na(as.matrix(out$data_imputed)))
      na_after_corr  <- sum(is.na(as.matrix(df)))
      new_na         <- na_after_corr - na_before_corr

      if (na_after_corr > 0L) {
        # df <- as.data.frame(.safe_impute(as.matrix(df), impute_fraction))
        imp <- run_mvimpute(
          x       = df,
          method  = impute_fraction,
          verbose = FALSE
        )

        df <- imp$data

        msg(sprintf(
          "  Re-imputed %d NA(s) after QCRSC (%d newly introduced, %d pre-existing).",
          na_after_corr, max(0L, new_na), max(0L, -new_na)
        ))
      }

      msg(sprintf("  QCRSC applied (%s).", corr_result$method_used))

    } else {
      msg("  Skipped (correct_drift = FALSE).")
    }

    out$data_corrected       <- df
    out$uncorrected_features <- uncorrected_features
    out$dimensions           <- .record_dim(out$dimensions, "After drift correction", df)

    # =========================================================================
    # 6. NORMALIZATION
    # =========================================================================

    msg(sprintf("Step 6/10: Normalizing (%s)...", normalize_method))

    norm_method_arg <- if (normalize_method == "none") "none" else normalize_method

    if (normalize_method == "none") {
      msg("  Skipped (normalize_method = 'none').")
      df_norm <- df
    } else {
      norm_result <- run_normalize(
        x              = df,
        metadata       = metadata,
        method         = norm_method_arg,
        factor_col     = norm_factor_col,
        ref_sample     = normalize_ref_sample,
        qc_normalize   = normalize_qc_method,
        groups         = group_col,
        qc_types       = qc_types,
        sample_id_col  = sample_id_col,
        verbose        = FALSE
      )
      df_norm <- norm_result$data
    }

    out$data_normalized_org <- df_norm # REMOVE AFTER CHECKING

    # =========================================================================
    # 7. TRANSFORMATION
    # =========================================================================

    msg(sprintf("Step 7/10: Transforming (%s)...", transform_method))

    if (transform_method == "none") {
      msg("  Skipped (transform_method = 'none').")
      df_trans <- df_norm
    } else {
      trans_result <- run_transform(
        x         = df_norm,
        method    = transform_method,
        metadata  = NULL,
        groups    = group_col,
        qc_types  = qc_types,
        num_cores = vsn_cores,
        verbose   = FALSE
      )
      df_trans <- trans_result$data
    }

    out$data_transformed_org <- df_trans # REMOVE AFTER CHECKING

    # =========================================================================
    # 8. SCALING (two parallel branches)
    # =========================================================================

    msg(sprintf("Step 8/10: Scaling (NONPLS = '%s', PLS = '%s')...",
                scale_nonpls, scale_pls))

    .scale_branch <- function(data, method) {
      if (method == "none") return(as.data.frame(data))
      run_scale(x = data, method = method, verbose = FALSE)$data
    }

    df_nonpls <- .scale_branch(df_trans, scale_nonpls)
    df_pls    <- .scale_branch(df_trans, scale_pls)

    out$data_nonpls_scaled_org <- df_nonpls # REMOVE AFTER CHECKING
    out$data_pls_scaled_org    <- df_pls    # REMOVE AFTER CHECKING

    # =========================================================================
    # 9. QUALITY FILTERING PER BRANCH
    # =========================================================================

    msg(sprintf(
      "Step 9/10: Quality filtering (RSD <= %.0f%%, bottom %g%% variance, ref = '%s')...",
      rsd_threshold * 100, variance_percentile, scale_filter_ref
    ))

    .filter_branch <- function(data, branch_name) {

      # RSD filter
      rsd_res  <- run_filterRSD(
        x        = data,
        metadata = metadata,
        max_rsd  = rsd_threshold,
        groups   = group_col,
        qc_type  = rsd_qc_type,
        qc_types = qc_types,
        verbose  = FALSE
      )
      data_rsd <- rsd_res$data

      msg(sprintf("  [%s] RSD: removed %d feature(s).",
                  branch_name, rsd_res$n_features_removed))

      # Variance filter
      var_res  <- run_filtervariance(
        x          = data_rsd,
        percentile = variance_percentile,
        verbose    = FALSE
      )
      data_var <- var_res$data

      msg(sprintf("  [%s] Variance: removed %d feature(s).",
                  branch_name, var_res$n_features_removed))

      list(rsd = data_rsd, final = data_var)
    }

    if (scale_filter_ref == "auto") {

      filt_nonpls <- .filter_branch(df_nonpls, "NONPLS")
      filt_pls    <- .filter_branch(df_pls,    "PLS")

      features_final <- intersect(
        colnames(filt_nonpls$final),
        colnames(filt_pls$final)
      )

      if (length(features_final) == 0L)
        stop("No features passed quality filtering in both branches. ",
             "Consider relaxing 'rsd_threshold' or 'variance_percentile'.")

      out$dimensions <- .record_dim(
        out$dimensions, "After NONPLS RSD filter", filt_nonpls$rsd
      )
      out$dimensions <- .record_dim(
        out$dimensions, "After NONPLS variance filter", filt_nonpls$final
      )
      out$dimensions <- .record_dim(
        out$dimensions, "After PLS RSD filter", filt_pls$rsd
      )
      out$dimensions <- .record_dim(
        out$dimensions, "After PLS variance filter", filt_pls$final
      )

    } else if (scale_filter_ref == "NONPLS") {

      filt_nonpls    <- .filter_branch(df_nonpls, "NONPLS")
      features_final <- colnames(filt_nonpls$final)

      out$dimensions <- .record_dim(
        out$dimensions, "After NONPLS RSD filter", filt_nonpls$rsd
      )
      out$dimensions <- .record_dim(
        out$dimensions, "After NONPLS variance filter (reference)", filt_nonpls$final
      )

    } else {  # "PLS"

      filt_pls       <- .filter_branch(df_pls, "PLS")
      features_final <- colnames(filt_pls$final)

      out$dimensions <- .record_dim(
        out$dimensions, "After PLS RSD filter", filt_pls$rsd
      )
      out$dimensions <- .record_dim(
        out$dimensions, "After PLS variance filter (reference)", filt_pls$final
      )
    }

    msg(sprintf("  Final feature set: %d feature(s).", length(features_final)))

    # Subset both scaled branches to the common feature set
    out$data_nonpls      <- df_nonpls[, features_final, drop = FALSE]
    out$data_pls         <- df_pls[,    features_final, drop = FALSE]

    # Subset pre-scale data to the same feature set (fold-change ready)
    out$data_normalized  <- as.data.frame(df_norm)[,  features_final, drop = FALSE]
    out$data_transformed <- as.data.frame(df_trans)[, features_final, drop = FALSE]

    out$features_final   <- features_final

    out$dimensions <- .record_dim(out$dimensions, "Final (NONPLS)", out$data_nonpls)
    out$dimensions <- .record_dim(out$dimensions, "Final (PLS)",    out$data_pls)

    # =========================================================================
    # 10. TECHNICAL-REPLICATE MERGING (base R)
    # =========================================================================

    msg("Step 10/10: Replicate merging...")

    do_merge <- merge_replicates &&
      subject_id_col %in% colnames(metadata) &&
      !all(is.na(meta_snap$SubjectID) | meta_snap$SubjectID == "" |
             meta_snap$SubjectID == "NA")

    if (do_merge) {

      # Identify rows: QC (never merged), biological with SubjectID, others
      sid        <- meta_snap$SubjectID
      is_qc      <- qc_rows
      has_sid    <- !is.na(sid) & sid != "" & sid != "NA"
      mergeable  <- !is_qc & has_sid

      unique_sids <- unique(sid[mergeable])
      has_dups    <- any(tabulate(match(sid[mergeable], unique_sids)) > 1L)

      if (!has_dups) {
        msg("  No technical replicates detected; skipping merge.")
        out$data_nonpls_merged <- NULL
        out$data_pls_merged    <- NULL
        out$metadata_merged    <- NULL
      } else {

        .merge_branch <- function(data) {
          bio_data  <- data[mergeable,  , drop = FALSE]
          keep_data <- data[!mergeable, , drop = FALSE]
          bio_sid   <- sid[mergeable]

          # Average features per SubjectID, preserving column types
          grp_idx   <- split(seq_len(nrow(bio_data)), bio_sid)
          merged_rows <- do.call(rbind, lapply(grp_idx, function(idx) {
            colMeans(bio_data[idx, , drop = FALSE], na.rm = TRUE)
          }))
          merged_df          <- as.data.frame(merged_rows)
          rownames(merged_df) <- names(grp_idx)

          rbind(merged_df, keep_data)
        }

        nonpls_merged <- .merge_branch(out$data_nonpls)
        pls_merged    <- .merge_branch(out$data_pls)

        # Build merged metadata: one row per merged subject + all unmerged rows
        mergeable_meta <- meta_snap[mergeable, , drop = FALSE]
        keep_meta      <- meta_snap[!mergeable, , drop = FALSE]

        grp_idx_meta <- split(seq_len(nrow(mergeable_meta)), sid[mergeable])
        merged_meta_rows <- do.call(rbind, lapply(grp_idx_meta, function(idx) {
          r <- mergeable_meta[idx[1L], , drop = FALSE]
          # Use the earliest injection order as the representative
          r$InjectionSequence <- min(
            mergeable_meta$InjectionSequence[idx], na.rm = TRUE
          )
          r
        }))
        rownames(merged_meta_rows) <- names(grp_idx_meta)

        merged_meta <- rbind(merged_meta_rows, keep_meta)

        out$data_nonpls_merged <- nonpls_merged
        out$data_pls_merged    <- pls_merged
        out$metadata_merged    <- merged_meta

        n_before <- sum(mergeable)
        n_after  <- nrow(merged_meta_rows)
        msg(sprintf("  Merged %d rows into %d unique subject(s).",
                    n_before, n_after))

        out$dimensions <- .record_dim(
          out$dimensions, "After replicate merge (NONPLS)", nonpls_merged
        )
        out$dimensions <- .record_dim(
          out$dimensions, "After replicate merge (PLS)", pls_merged
        )
      }

    } else {
      if (merge_replicates && !do_merge)
        msg("  merge_replicates = TRUE but no valid SubjectID data found; skipping.")
      else
        msg("  Skipped (merge_replicates = FALSE).")

      out$data_nonpls_merged <- NULL
      out$data_pls_merged    <- NULL
      out$metadata_merged    <- NULL
    }

    msg("Pipeline completed successfully.")

  }, error = function(e) {
    message("run_DIpreprocess ERROR: ", e$message)
    out$error <<- e$message
  })

  out$elapsed_seconds <- proc.time()[["elapsed"]] - t_start
  msg(sprintf("Elapsed: %.1f second(s).", out$elapsed_seconds))

  class(out) <- c("run_DIpreprocess", "list")
  out
}


# ==============================================================================
# S3 METHODS
# ==============================================================================

#' @export
print.run_DIpreprocess <- function(x, ...) {
  cat("=== run_DIpreprocess pipeline result ===\n")

  if (!is.null(x$error)) {
    cat(sprintf("STATUS : FAILED — %s\n", x$error))
    return(invisible(x))
  }

  cat("STATUS : OK\n")

  n_samp <- if (!is.null(x$data_nonpls)) nrow(x$data_nonpls) else NA_integer_
  n_feat <- length(x$features_final)
  cat(sprintf("Data   : %d samples × %d features (final)\n", n_samp, n_feat))

  if (!is.null(x$data_nonpls_merged))
    cat(sprintf("Merged : %d samples × %d features (replicate-merged)\n",
                nrow(x$data_nonpls_merged), ncol(x$data_nonpls_merged)))

  if (length(x$uncorrected_features) > 0L)
    cat(sprintf("QCRSC  : %d uncorrected feature(s) retained\n",
                length(x$uncorrected_features)))

  cat(sprintf("Time   : %.1f second(s)\n", x$elapsed_seconds))
  invisible(x)
}


#' @export
summary.run_DIpreprocess <- function(object, ...) {
  cat("=== run_DIpreprocess — step-by-step dimensions ===\n")
  print(object$dimensions, row.names = FALSE)

  if (!is.null(object$data_nonpls) && !is.null(object$data_pls)) {
    cat("\nFinal datasets:\n")
    cat(sprintf("  data_nonpls      : %d × %d\n",
                nrow(object$data_nonpls), ncol(object$data_nonpls)))
    cat(sprintf("  data_pls         : %d × %d\n",
                nrow(object$data_pls), ncol(object$data_pls)))
    cat(sprintf("  data_normalized  : %d × %d (pre-scale)\n",
                nrow(object$data_normalized), ncol(object$data_normalized)))
    cat(sprintf("  data_transformed : %d × %d (pre-scale)\n",
                nrow(object$data_transformed), ncol(object$data_transformed)))

    if (!is.null(object$data_nonpls_merged))
      cat(sprintf("  data_nonpls_merged: %d × %d\n",
                  nrow(object$data_nonpls_merged), ncol(object$data_nonpls_merged)))
  }

  invisible(object)
}