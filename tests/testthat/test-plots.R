# Each plot must construct AND build (ggplot_build catches rendering errors)

expect_built <- function(p) {
  expect_true(ggplot2::is_ggplot(p))
  expect_no_error(ggplot2::ggplot_build(p))
}

iv <- droplevels(iris[iris$Species != "setosa", ])

test_that("basic plotting functions return buildable ggplots", {
  expect_built(plot_bar(mtcars, x_axis = "cyl"))
  expect_built(plot_bar(as.data.frame(Titanic), x_axis = "Class", by = "Sex",
                        plot_type = "dodged"))
  expect_built(plot_oddsratio(data.frame(OR = c(1.4, 0.8), row.names = c("a", "b")),
                              or = "OR"))
  expect_built(plot_multi_y(data.frame(w = 1:5, a = 1:5, b = 5:1), x_axis = "w",
                            main_as_area = "a", main_as_line = "b"))
  expect_built(plot_heatmap(mtcars[, 1:5]))
})

test_that("plot_dend returns one or two dendrograms", {
  skip_if_not_installed("ggdendro")
  expect_built(plot_dend(mtcars))
  both <- plot_dend(mtcars, orientation = "both")
  expect_named(both, c("rows", "cols"))
  expect_built(both$cols)
})

test_that("plot_meanmap returns a patchwork", {
  skip_if_not_installed("ggdendro")
  p <- plot_meanmap(iris[, 1:4], data.frame(Species = iris$Species), group = "Species")
  expect_s3_class(p, "patchwork")
})

test_that("plot_diff returns one plot per outcome", {
  r <- suppressMessages(plot_diff(PlantGrowth, outcome = "weight", group = "group"))
  expect_named(r$plots, "weight")
  expect_built(r$plots$weight)
})

test_that("plot_volcano classifies features against thresholds", {
  d <- data.frame(f = paste0("f", 1:4), fc = c(4, 0.2, 1, 3), p = c(0.001, 0.001, 0.5, 0.5))
  v <- suppressMessages(plot_volcano(d, y = "fc", z = "p", features = "f",
                                     up = 2, down = 0.5))
  expect_s3_class(v, "plot_volcano")
  expect_built(v$plots[[1]])
  cl <- v$classified
  expect_true(is.data.frame(cl))
  expect_equal(cl$Regulation, c("Up", "Down", "NS", "NS"))
})

test_that("plot_fcmetric builds from run_foldchange", {
  fc <- run_foldchange(iv[, 1:4], metadata = iv["Species"], group = "Species")
  r  <- suppressMessages(plot_fcmetric(fc))
  expect_s3_class(r, "plot_fcmetric")
  expect_built(r$plots[[1]])
})

test_that("plot_score and plot_scree build from run_pca", {
  meta <- data.frame(Sample = rownames(iris), Species = iris$Species)
  pca <- run_pca(scale(iris[, 1:4]), meta, group = "Species", verbose = FALSE)
  expect_built(suppressMessages(plot_score(pca, color_by = "Species",
                                           points_from = "Sample")))
  expect_built(plot_scree(pca, verbose = FALSE))
})

test_that("plot_score builds from run_pcoa", {
  skip_if_not_installed("ape")
  skip_if_not_installed("vegan")
  p <- run_pcoa(as.matrix(iris[1:30, 1:4]), method = "euclidean")
  m <- data.frame(Sample = rownames(iris)[1:30], Species = iris$Species[1:30])
  expect_built(suppressMessages(plot_score(p, metadata = m, color_by = "Species",
                                           points_from = "Sample")))
})

test_that("plot_mixedmodel builds from run_mixedmodel", {
  skip_if_not_installed("lmerTest")
  set.seed(4)
  meta <- data.frame(subject_id = factor(rep(1:15, each = 2)),
                     time = rep(c("pre", "post"), 15))
  x <- data.frame(g = rnorm(30), i = rnorm(30))
  m <- suppressMessages(run_mixedmodel(x, meta, verbose = FALSE,
                        formula = ".feature ~ time + (1 | subject_id)"))
  expect_built(plot_mixedmodel(m))
})

test_that("plot_relabund normalises counts to 100%", {
  counts <- data.frame(t1 = c(1, 3, 2, 2), t2 = c(3, 1, 2, 2))
  meta <- data.frame(SampleID = paste0("s", 1:4))
  p <- suppressMessages(plot_relabund(counts, metadata = meta))
  expect_built(p)
  b <- ggplot2::ggplot_build(p)$data[[1]]
  tot <- tapply(b$ymax - b$ymin, b$x, sum)
  expect_true(all(abs(tot - max(tot)) < 1e-6))        # every bar has the same total
})

test_that("plot_relabund group bars average per-sample proportions", {
  # Unequal sample totals: summing raw counts would give A = 0.48/0.52
  counts <- data.frame(t1 = c(10, 200, 5, 1), t2 = c(30, 200, 5, 9))
  meta <- data.frame(SampleID = paste0("s", 1:4), Group = c("A", "A", "B", "B"))
  p <- suppressWarnings(suppressMessages(
    plot_relabund(counts, metadata = meta, show_x_as = "Group", group_aggregate = TRUE)))
  b <- ggplot2::ggplot_build(p)$data[[1]]
  h <- b$ymax - b$ymin
  props <- unlist(tapply(h, b$x, function(v) sort(v / sum(v))))
  expect_equal(unname(sort(props)), sort(c(0.375, 0.625, 0.3, 0.7)), tolerance = 1e-6)
})
