# File export (always inside a temporary directory)

test_that("save_excel writes one sheet per list element and handles an empty list", {
  d <- withr::local_tempdir()
  save_excel(list(Cars = mtcars, Iris = iris), filename = "two", folder = d)
  f <- list.files(d, "\\.xlsx$", full.names = TRUE, recursive = TRUE)
  expect_length(f, 1)
  expect_setequal(openxlsx::getSheetNames(f), c("Cars", "Iris"))
  back <- openxlsx::read.xlsx(f, sheet = "Cars")
  expect_equal(nrow(back), nrow(mtcars))

  expect_warning(save_excel(list(), filename = "empty", folder = d), "empty list")
  expect_length(list.files(d, "empty.*\\.xlsx$", recursive = TRUE), 1)
})

test_that("save_excel exports pondeR result tables", {
  d <- withr::local_tempdir()
  save_excel(run_correl(mtcars[, 1:3], method = "pearson"), filename = "cor", folder = d)
  expect_length(list.files(d, "\\.xlsx$", recursive = TRUE), 1)
})

test_that("save_plots writes ggplots and base-R plot functions", {
  d <- withr::local_tempdir()
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  save_plots(list(sc = p, h = function() graphics::hist(mtcars$mpg)), folder = d,
             format = "png", width = 3, height = 3, dpi = 50)
  expect_setequal(basename(list.files(d, "\\.png$", recursive = TRUE)), c("sc.png", "h.png"))
})
