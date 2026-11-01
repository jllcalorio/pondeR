# Extracted from test-render_metabolomics.R:65

# prequel ----------------------------------------------------------------------
skip_if_not_installed("dimsprepr")
make_files <- function(dir) {
  set.seed(2741)
  n <- 40
  meta <- data.frame(Sample = paste0("S", 1:n),
                     Group = c(rep("QC", 8), rep(c("A", "B"), each = 16)),
                     Batch = rep(1:2, each = 20), InjectionSequence = 1:n)
  feat <- as.data.frame(matrix(rlnorm(n * 10, 11, 0.3), n, 10,
                               dimnames = list(NULL, paste0("m", 1:10))))
  utils::write.csv(meta, file.path(dir, "meta.csv"), row.names = FALSE)
  utils::write.csv(cbind(Sample = meta$Sample, feat), file.path(dir, "feat.csv"), row.names = FALSE)
  utils::write.csv(feat, file.path(dir, "feat_noid.csv"), row.names = FALSE)
  utils::write.csv(cbind(meta, feat), file.path(dir, "combined.csv"), row.names = FALSE)
  invisible(list(meta = meta, feat = feat))
}

# test -------------------------------------------------------------------------
skip_on_cran()
skip_if_not_installed("quarto")
skip_if(is.null(quarto::quarto_path()))
d <- withr::local_tempdir()
make_files(d)
out <- suppressMessages(render_metabolomics(
    file.path(d, "feat.csv"), file.path(d, "meta.csv"), qc_types = "QC",
    run_DIpreprocess = list(min_qc_per_batch = 3)))
