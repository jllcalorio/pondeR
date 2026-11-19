# Ordination and dimension reduction

meta_iris <- data.frame(Sample = rownames(iris), Species = iris$Species)

test_that("run_pca matches prcomp and works without dimsprepr", {
  pca <- run_pca(scale(iris[, 1:4]), meta_iris, group = "Species", verbose = FALSE)
  expect_s3_class(pca, "run_pca")
  ref <- summary(prcomp(iris[, 1:4], scale. = TRUE))$importance[2, ]
  expect_equal(unname(pca$variance_explained[1:4]) / 100, unname(ref), tolerance = 1e-3)
  expect_equal(pca$n_samples, 150)

  ex <- run_pca(scale(iris[, 1:4]), meta_iris, group = "Species",
                exclude = "setosa", verbose = FALSE)
  expect_equal(ex$n_samples, 100)
})

test_that("run_pcoa on Euclidean distances equals PCA", {
  skip_if_not_installed("ape")
  skip_if_not_installed("vegan")
  p <- run_pcoa(as.matrix(iris[1:30, 1:4]), method = "euclidean")
  expect_s3_class(p, "run_pcoa")
  ref <- summary(prcomp(iris[1:30, 1:4]))$importance[2, 1:3] * 100
  expect_equal(unname(p$variance_explained[1:3]), unname(ref), tolerance = 1e-3)
})

test_that("run_pcoa PERMANOVA detects group separation", {
  skip_on_cran()
  skip_if_not_installed("ape")
  skip_if_not_installed("vegan")
  set.seed(1)
  m <- data.frame(Species = iris$Species)
  p <- suppressMessages(run_pcoa(as.matrix(iris[, 1:4]), method = "euclidean",
                                 metadata = m, rhs = "Species", perms = 199))
  expect_lt(p$permanova$p_value, 0.05)
  expect_gt(p$permanova$R2, 0.5)
})

test_that("run_pls fits PLS-DA with pre-scaled data", {
  skip_on_cran()
  skip_if_not_installed("ropls")
  m <- droplevels(meta_iris[51:150, ])
  pls <- suppressWarnings(run_pls(scale(iris[51:150, 1:4]), m, method = "plsda",
                                  group = "Species", verbose = FALSE))
  expect_s3_class(pls, "run_pls")
  expect_equal(pls$n_samples, 100)
  expect_true(all(c("Petal.Length", "Petal.Width") %in% names(pls$vip_scores)))
})
