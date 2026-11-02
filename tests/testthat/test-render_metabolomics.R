skip_if_not_installed("dimsprepr")

make_files <- function(dir) {
  set.seed(2741)
  n <- 120
  grp <- rep("QC", n)  # a QC every 5th injection: enough QCs per batch for drift correction
  grp[-seq(1, n, by = 5)] <- sample(rep(c("A", "B"), each = 48))
  meta <- data.frame(Sample = paste0("S", 1:n), Group = grp,
                     Batch = sort(rep(1:3, length.out = n)), InjectionSequence = 1:n)
  feat <- as.data.frame(matrix(rlnorm(n * 10, 11, 0.3), n, 10,
                               dimnames = list(NULL, paste0("m", 1:10))))
  feat[grp == "B", 1:2] <- feat[grp == "B", 1:2] * 3
  utils::write.csv(meta, file.path(dir, "meta.csv"), row.names = FALSE)
  utils::write.csv(cbind(Sample = meta$Sample, feat), file.path(dir, "feat.csv"), row.names = FALSE)
  utils::write.csv(feat, file.path(dir, "feat_noid.csv"), row.names = FALSE)
  utils::write.csv(cbind(meta, feat), file.path(dir, "combined.csv"), row.names = FALSE)
  invisible(list(meta = meta, feat = feat))
}

test_that(".wf_prepare aligns separate, id-less, transposed and combined inputs", {
  d <- withr::local_tempdir(); m <- make_files(d)
  p <- function(...) suppressMessages(pondeR:::.wf_prepare(
    ..., sample_id_col = "Sample", group_col = "Group", batch_col = "Batch",
    injection_col = "InjectionSequence"))

  a <- p(file.path(d, "feat.csv"), file.path(d, "meta.csv"), NULL, NULL, NULL)
  expect_equal(rownames(a$features), m$meta$Sample)
  expect_equal(ncol(a$features), 10)

  b <- p(file.path(d, "feat_noid.csv"), file.path(d, "meta.csv"), NULL, NULL, NULL)
  expect_equal(unname(as.matrix(b$features)), unname(as.matrix(a$features)))

  cmb <- p(file.path(d, "combined.csv"), NULL, NULL, NULL, NULL)
  expect_equal(ncol(cmb$features), 10)

  tr <- data.frame(feature = names(m$feat), t(m$feat), check.names = FALSE)
  names(tr)[-1] <- m$meta$Sample
  t_ <- p(tr, file.path(d, "meta.csv"), NULL, NULL, NULL)
  expect_equal(unname(as.matrix(t_$features)), unname(as.matrix(m$feat)))
})

test_that(".wf_prepare gives clear errors", {
  d <- withr::local_tempdir(); make_files(d)
  p <- function(f, md, ...) pondeR:::.wf_prepare(f, md, NULL, NULL, NULL, "Sample", "Group",
                                                 "Batch", "InjectionSequence")
  expect_error(p(file.path(d, "nope.csv"), file.path(d, "meta.csv")), "not found")
  expect_error(p(file.path(d, "feat.csv"), data.frame(Id = 1:3)), "no sample ID column")
  expect_error(p(file.path(d, "feat_noid.csv")[1], data.frame(Sample = "x", Group = "A",
                                                              Batch = 1, InjectionSequence = 1)),
               "Cannot match")
})

test_that("function overrides are validated", {
  expect_equal(pondeR:::.wf_check_overrides(list(run_foldchange = list(up = 2)))$run_foldchange$up, 2)
  expect_error(pondeR:::.wf_check_overrides(list(run_foldchange = list(upp = 2))), "not allowed")
  expect_error(pondeR:::.wf_check_overrides(list(run_foo = list(a = 1))), "Unknown function")
  expect_error(pondeR:::.wf_check_overrides(list(run_pca = c(1, 2))), "named list")
  expect_error(pondeR:::.wf_check_overrides(list(run_pca = list(x = 1))), "not allowed")
  expect_equal(pondeR:::.wf_merge(list(a = 1, b = 2), list(b = NULL))$b, NULL)
})

test_that("render_metabolomics writes one standalone report into results/", {
  skip_on_cran(); skip_if_not_installed("quarto"); skip_if(is.null(quarto::quarto_path()))
  # Quarto renders in a fresh R process, which needs the *installed* current pondeR
  inst <- find.package("pondeR", lib.loc = .libPaths(), quiet = TRUE)
  same <- length(inst) == 1 && file.exists(file.path(inst, "workflows", "metabolomics_report.qmd")) &&
    identical(readLines(file.path(inst, "workflows", "metabolomics_report.qmd")),
              readLines(system.file("workflows", "metabolomics_report.qmd", package = "pondeR")))
  skip_if_not(same,
              "install the current pondeR (devtools::install()) to run this test")
  d <- withr::local_tempdir(); make_files(d)
  out <- suppressMessages(render_metabolomics(
    file.path(d, "feat.csv"), file.path(d, "meta.csv"), qc_types = "QC"))
  expect_true(file.exists(out))
  expect_equal(normalizePath(dirname(out)), normalizePath(file.path(d, "results")))
  expect_equal(sort(list.files(file.path(d, "results"))), c("feat_report.html", "feat_results.rds"))
  expect_true(any(grepl("data:image/png;base64", readLines(out, warn = FALSE), fixed = TRUE)))
})

test_that("output names are stamped instead of overwriting unless asked", {
  d <- withr::local_tempdir()
  f <- pondeR:::.wf_output_name
  expect_equal(f(d, "x", FALSE, FALSE), "x")
  expect_match(f(d, "x", FALSE, TRUE), "^x_[0-9]{8}-[0-9]{6}$")
  file.create(file.path(d, "x_report.html"))
  expect_match(f(d, "x", FALSE, FALSE), "^x_[0-9]{8}-[0-9]{6}$")
  expect_equal(f(d, "x", TRUE, FALSE), "x")
})

test_that("plot_score.pca / plot_score.pls overrides are accepted and layered over plot_score", {
  ok <- pondeR:::.wf_check_overrides(list(plot_score = list(title = "a"),
                                          plot_score.pca = list(color_by = "Batch")))
  expect_named(ok, c("plot_score", "plot_score.pca"))
  expect_error(pondeR:::.wf_check_overrides(list(plot_score.pca = list(res = 1))), "not allowed")
  # same underlying function; the specific override wins over the shared one
  testthat::local_mocked_bindings(getExportedValue = function(ns, name) function(...) list(...),
                                  .package = "base")
  got <- pondeR:::.wf_call("plot_score.pca", list(a = 1, title = "default"),
                           list(plot_score = list(title = "shared", b = 2),
                                plot_score.pca = list(title = "pca")))
  expect_equal(got$title, "pca"); expect_equal(got$b, 2); expect_equal(got$a, 1)
})
