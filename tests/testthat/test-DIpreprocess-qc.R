skip_if_not_installed("dimsprepr")

make_di <- function() {
  set.seed(3917)
  n <- 80; p <- 30
  grp <- rep("QC", n)
  grp[-seq(1, n, by = 5)] <- sample(rep(c("Control", "Case"), each = 32))
  meta <- data.frame(Sample = paste0("S", 1:n), Group = grp,
                     Batch = sort(rep(1:2, length.out = n)), InjectionSequence = 1:n)
  x <- matrix(rlnorm(n * p, log(1e5), 0.3), n, p,
              dimnames = list(meta$Sample, paste0("m", 1:p)))
  x[grp == "Case", seq(1, p, by = 6)] <- x[grp == "Case", seq(1, p, by = 6)] * 4
  rownames(meta) <- meta$Sample
  suppressMessages(dimsprepr::run_DIpreprocess(
    as.data.frame(x), meta, qc_types = "QC", verbose = FALSE))
}

test_that("run_pls drops QC samples by default for run_DIpreprocess input", {
  di <- make_di()
  p <- suppressWarnings(suppressMessages(
    run_pls(di, group = "Group", method = "plsda", verbose = FALSE)))
  expect_equal(p$n_samples, sum(di$metadata$Group != "QC"))
  expect_false("QC" %in% p$metadata$Group)
})

test_that("run_auc drops QC samples for a bare run_DIpreprocess object", {
  di <- make_di()
  r <- suppressWarnings(suppressMessages(run_auc(di)))
  expect_equal(r$n_cases + r$n_controls, sum(di$metadata$Group != "QC"))
})

test_that("plot_volcano falls back to mean abundance without run_diff results", {
  di <- make_di()
  fc <- suppressMessages(run_foldchange(di, from = "normalized"))
  expect_message(v <- plot_volcano(list(fc), y = "fold_change", z = "p_value"),
                 "mean abundance")
  expect_s3_class(v, "plot_volcano")
  expect_named(v$plots, "Mean Abundance")
})

test_that("DI-aware run_* functions drop QC samples by default", {
  di <- make_di()
  n_study <- sum(di$metadata$Group != "QC")
  f1 <- colnames(di$data_nonpls)[1]
  q <- function(e) suppressWarnings(suppressMessages(e))

  fc <- q(run_foldchange(di, from = "normalized"))
  expect_equal(sum(unique(fc$summary_table[c("n_num", "n_den")])), n_study)
  expect_equal(sum(q(run_diff(di, outcome = f1))$group_ns), n_study)
  # logreg would error on 3 outcome levels if QC were kept
  expect_s3_class(q(run_logreg(di, y = "Group", indep = f1, ref = "Control")), "data.frame")
  mm <- q(run_mixedmodel(di, formula = ".feature ~ Batch + (1|Group)"))
  expect_equal(mm[[1]]$n_obs, n_study)
})
