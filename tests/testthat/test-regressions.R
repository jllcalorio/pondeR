# One check per bug fixed during the CRAN review

iv <- droplevels(iris[iris$Species != "setosa", ])

test_that("every text size follows base_size, and explicit sizes win", {
  s <- .text_sizes(10)
  expect_equal(c(s$title, s$xlab, s$axis_text, s$label), c(12, 10, 9, 8))
  expect_equal(.text_sizes(10, title = 20)$title, 20)
  expect_error(.text_sizes(-1), "base_size")
  expect_error(.ponder_theme("nope"), "theme")

  th <- plot_bar(mtcars, x_axis = "cyl", base_size = 20, title = "t")$theme
  expect_equal(th$plot.title$size, 24)
  expect_equal(th$axis.title.x$size, 20)
})

test_that("no plot function keeps the old size arguments", {
  plots <- c("plot_bar", "plot_dend", "plot_diff", "plot_fcmetric", "plot_heatmap",
             "plot_meanmap", "plot_mixedmodel", "plot_multi_y", "plot_oddsratio",
             "plot_relabund", "plot_scree", "plot_volcano", "run_auc")
  for (fn in plots) {
    a <- formals(get(fn))
    expect_false(any(c("global_font_size", "zoom") %in% names(a)), info = fn)
    expect_equal(a$base_size, 11, info = fn)
  }
  expect_equal(formals(plot_score.run_pca)$base_size, 11)
})

test_that("point labels are sized in points, not millimetres", {
  d <- data.frame(f = c("a", "b"), fc = c(4, 0.2), p = c(0.001, 0.001))
  v <- plot_volcano(d, y = "fc", z = "p", features = "f", annotate = "a", use_ggrepel = FALSE)
  lab <- ggplot2::ggplot_build(v$plots[[1]])$data[[5]]
  expect_equal(lab$size[1], 11 * 0.8 / ggplot2::.pt)
})

test_that("plot_diff reads formatted p-values as numbers", {
  r <- run_diff(iv[, 1:2], iv["Species"], names(iv)[1:2], "Species",
                summary_table = TRUE, verbose = FALSE)
  expect_type(r$summary_table[["p-value"]], "character")
  p <- plot_diff(r, outcome = "Sepal.Width")
  expect_named(p$plots, "Sepal.Width")
  sw <- t.test(Sepal.Width ~ Species, iv)$p.value
  expect_equal(length(p$significant_plots) == 1L, p.adjust(c(sw, sw), "BH")[1] < 0.05)
})

test_that("plot_meanmap rows line up with the dendrogram leaves", {
  skip_if_not_installed("ggdendro")
  p <- plot_meanmap(iris[, 1:4], data.frame(Species = iris$Species), group = "Species")
  hc <- attr(p, "hclust_obj")
  heat <- p[[2]]$data
  expect_equal(levels(heat$variable), hc$labels[hc$order])
})

test_that("printing and summarising a multi-outcome run_diff works", {
  r <- run_diff(iv[, 1:2], iv["Species"], names(iv)[1:2], "Species",
                summary_table = TRUE, verbose = FALSE)
  expect_output(print(r), "Adjusted p-value")
  expect_output(summary(r), "Variable")
  r2 <- run_diff(iv[, 1:2], iv["Species"], names(iv)[1:2], "Species", verbose = FALSE)
  expect_output(print(r2), "2 outcomes")
})

test_that("run_diff tables carry report headers and still feed plot_diff", {
  r <- run_diff(iris, outcome = c("Sepal.Length", "Sepal.Width"), group = "Species",
                summary_table = TRUE)
  one <- r$Sepal.Length
  expect_true(all(c("Group", "Mean", "SD") %in% names(one$data_summary)))
  expect_true(all(c("Group 1", "Group 2", "Adjusted p-value") %in% names(one$posthoc_result)))
  expect_false(any(c("p.format", "p.adj.format") %in% names(one$posthoc_result)))
  layers <- plot_diff(r, outcome = "Sepal.Length")$plots[[1]]$layers
  expect_true(any(vapply(layers, function(l) inherits(l$geom, "GeomBracket"), TRUE)))
})

test_that("run_diff honours calculate_effect_size = FALSE", {
  r <- run_diff(iris, outcome = "Sepal.Length", group = "Species",
                calculate_effect_size = FALSE)
  expect_null(r$effect_size)
  expect_output(summary(r), "Post-Hoc")
  m <- run_diff(iris, outcome = c("Sepal.Length", "Sepal.Width"), group = "Species",
                calculate_effect_size = FALSE, summary_table = TRUE)
  expect_false("ES (95% CI)" %in% names(m$summary_table))
})

test_that("run_pls fits without drawing to a graphics device", {
  skip_if_not_installed("ropls")
  grDevices::graphics.off()
  m <- data.frame(Sample = rownames(iv), Species = iv$Species)
  run_pls(scale(iv[1:4]), m, method = "plsda", group = "Species", permI = 0)
  expect_null(grDevices::dev.list())
})

test_that("auto correlation picks Pearson for normal variables on different scales", {
  set.seed(1)
  a <- rnorm(40)
  d <- data.frame(a = a, b = 100 * a + rnorm(40, sd = 50))
  r <- run_correl(d, summary_table_only = TRUE)
  expect_equal(r$Method, "Pearson")
})

test_that("rank-correlation CIs use the Fieller standard errors", {
  r <- run_correl(mtcars[c("mpg", "wt")], method = "spearman", summary_table_only = TRUE)
  rho <- cor(mtcars$mpg, mtcars$wt, method = "spearman")
  se <- sqrt(1.06 / (nrow(mtcars) - 3))
  expect_equal(r[["95% CI lower"]], tanh(atanh(rho) - qnorm(0.975) * se), tolerance = 1e-8)
})

test_that("default calls do not touch the random seed", {
  set.seed(1); before <- .Random.seed
  run_subset(iris, n = 10, random = FALSE)
  expect_identical(.Random.seed, before)
  expect_null(formals(run_regreg)$seed)
  expect_null(formals(run_randomize)$seed)
  expect_null(formals(run_subset)$seed)
})

test_that("save functions refuse to write without a destination", {
  expect_error(save_excel(mtcars), "where to save")
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  expect_error(save_plots(p), "where to save")
})

test_that("run_relabel recodes simultaneously and merges labels", {
  d <- data.frame(g = c("A", "B", "C"))
  expect_equal(as.character(run_relabel(d, "g", c(A = "B", B = "C"))$g), c("B", "C", "C"))
  expect_equal(levels(run_relabel(d, "g", c(A = "X", B = "X"))$g), c("X", "C"))
})

test_that("run_reduce results print with their own method", {
  a <- data.frame(x = 1, y = 2); b <- data.frame(y = 3, z = 4)
  r <- run_reduce(a, b)
  expect_s3_class(r, "run_reduce")
  expect_output(print(r), "run_reduce result")
})

test_that("run_lefse tolerates missing values and warns cleanly", {
  x <- data.frame(a = c(NA, 2:10, 31:40), b = c(5, 3, 8, 1, 9, 2, 7, 4, 6, 10, 3, 8, 1, 9, 2, 7, 4, 6, 10, 5))
  meta <- data.frame(Group = rep(c("A", "B"), each = 10))
  r <- suppressWarnings(run_lefse(x, meta))
  expect_true("a" %in% r$lda$Feature)
  w <- tryCatch(run_lefse(data.frame(a = c(1, 2, 1, 2), b = c(2, 1, 2, 1)),
                          data.frame(Group = c("A", "A", "B", "B"))),
                warning = conditionMessage)
  expect_false(grepl("FALSE", w))
})

test_that("get_groupsizes can return a summary table", {
  skip_if_not_installed("gt")
  meta <- data.frame(G = rep(c("a", "b"), 10), A = 1:20)
  r <- get_groupsizes(meta, pred = "G", show_summary_table = TRUE)
  expect_named(r, c("group_sizes", "summary_table"))
  expect_s3_class(r$summary_table, "gt_tbl")
})

test_that("plot_relabund honours ignore_normalization and does not print", {
  counts <- data.frame(t1 = c(1, 3), t2 = c(3, 1))
  meta <- data.frame(SampleID = c("s1", "s2"))
  p <- plot_relabund(counts, metadata = meta, ignore_normalization = TRUE)
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(max(b$ymax), 4)
  expect_true(withVisible(plot_relabund(counts, metadata = meta))$visible)
})

test_that("subgroup summary tables get a clean automatic title", {
  skip_if_not_installed("gt")
  d <- data.frame(s = rep(c("m", "f"), 10), g = rep(c("a", "b"), each = 10), v = rnorm(20))
  tb <- run_summarytable(d, split_by = "g", subgroup = "s")
  expect_false(grepl("auto", tb[[1]][["_heading"]]$title))
})

# dimsprepr ------------------------------------------------------------------

di_data <- function(...) {
  set.seed(3917)
  n <- 60; p <- 20
  grp <- rep("QC", n)
  grp[-seq(1, n, by = 5)] <- sample(rep(c("Control", "Case"), each = 24))
  meta <- data.frame(Sample = paste0("S", 1:n), Class = grp,
                     Batch = sort(rep(1:2, length.out = n)), InjectionSequence = 1:n,
                     SubjectID = paste0("B", ceiling(seq_len(n) / 2)))
  rownames(meta) <- meta$Sample
  x <- as.data.frame(matrix(rlnorm(n * p, log(1e5), 0.3), n, p,
                            dimnames = list(meta$Sample, paste0("m", 1:p))))
  suppressWarnings(suppressMessages(dimsprepr::run_DIpreprocess(
    x, meta, group_col = "Class", qc_types = "QC", verbose = FALSE, ...)))
}

test_that("dimsprepr results work when the group column is not named Group", {
  skip_if_not_installed("dimsprepr")
  di <- di_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  expect_s3_class(q(run_foldchange(di)), "run_foldchange")
  expect_s3_class(q(run_diff(di, outcome = colnames(di$data_nonpls)[1])), "run_diff")
  expect_s3_class(q(run_summarytable(di)), "gt_tbl")
  expect_s3_class(q(run_auc(di)), "run_auc")
})

test_that("run_foldchange pairs merged replicates with the right metadata", {
  skip_if_not_installed("dimsprepr")
  di <- di_data(merge_replicates = TRUE)
  skip_if(is.null(di$metadata_merged))
  expect_s3_class(run_foldchange(di, from = "normalized"), "run_foldchange")
  expect_s3_class(run_foldchange(di, from = "scaled"), "run_foldchange")
})

test_that("multi-method dimsprepr results give a clear error", {
  skip_if_not_installed("dimsprepr")
  di <- di_data(normalize_method = c("median", "sum"), drift_method = "none")
  expect_s3_class(di, "run_DIpreprocess_multi")
  expect_error(run_pca(di), "Pick one")
  expect_error(run_foldchange(di), "Pick one")
})
