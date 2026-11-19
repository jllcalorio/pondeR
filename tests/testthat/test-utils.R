# Checkers, selectors and data-wrangling helpers

test_that("is_tabular / is_list / has_names classify objects", {
  expect_true(is_tabular(mtcars))
  expect_true(is_tabular(matrix(1:4, 2)))
  expect_false(is_tabular(1:3, warn = FALSE))
  expect_false(is_tabular(NULL, warn = FALSE))
  expect_warning(is_tabular(list(a = 1)))

  expect_true(is_list(list(a = 1)))
  expect_false(is_list(mtcars, warn = FALSE))

  expect_true(has_names(mtcars))
  expect_false(has_names(list(1, 2), warn = FALSE))
  expect_false(has_names(list(a = 1, 2), warn = FALSE))
})

test_that("has_na / has_zero / has_negative / has_inf detect values", {
  expect_true(has_na(c(1, NA)))
  expect_false(has_na(1:3))
  expect_true(has_zero(data.frame(x = c(1, 0))))
  expect_false(has_zero(c(NA, 1), warn = FALSE))
  expect_true(has_negative(c(0, -0.1)))
  expect_false(has_negative(1:3))
  expect_true(has_inf(matrix(c(1, -Inf), 1)))
  expect_false(has_inf(c(NaN, NA, 1)))
})

test_that("select_numvars / select_catvars split columns by type", {
  expect_named(select_numvars(iris), names(iris)[1:4])
  expect_named(select_catvars(iris), "Species")
})

test_that("run_vector2df converts, filters and sorts a named vector", {
  v <- c(a = 0.5, b = 2, c = 1.2)
  out <- run_vector2df(v, name = "Feature", value = "VIP", min = 1,
                       sort_by = "value_desc")
  expect_named(out, c("Feature", "VIP"))
  expect_equal(out$Feature, c("b", "c"))
})

test_that("run_move_names moves row names into a column", {
  out <- run_move_names(mtcars[1:3, 1:2])
  expect_equal(ncol(out), 3)
  expect_equal(out[[1]], rownames(mtcars)[1:3])
})

test_that("run_relabel recodes and keeps unmapped levels", {
  out <- run_relabel(data.frame(s = c("A", "B", "D")), "s", c(A = "Active"))
  expect_s3_class(out$s, "factor")
  expect_true(all(c("Active", "B", "D") %in% levels(out$s)))
})

test_that("run_randomize samples the requested size reproducibly", {
  a <- run_randomize(mtcars, n = 5, seed = 1)
  b <- run_randomize(mtcars, n = 5, seed = 1)
  expect_equal(nrow(a), 5)
  expect_identical(a, b)
  expect_equal(nrow(run_randomize(mtcars, p = 0.25, seed = 1)), 8)
  g <- run_randomize(iris, per_group = list(list("Species", n = 3)), seed = 1)
  expect_equal(as.vector(table(g$Species)), c(3, 3, 3))
})

test_that("run_randomize returns intact rows by default", {
  r <- run_randomize(mtcars, n = 5, seed = 1)
  expect_equal(nrow(merge(r, mtcars)), 5)          # every row is a real row
  expect_identical(vapply(r, class, ""), vapply(mtcars, class, ""))
  s <- run_randomize(iris, per = "col", seed = 1)
  expect_identical(vapply(s, class, ""), vapply(iris, class, ""))
  expect_equal(sort(s$Sepal.Length), sort(iris$Sepal.Length))
})

test_that("run_reduce keeps common or all columns", {
  d1 <- data.frame(a = 1, b = 2, c = 3)
  d2 <- data.frame(a = 4, b = 5, d = 6)
  common <- run_reduce(d1, d2)
  expect_named(common$data, c("a", "b"))
  expect_equal(nrow(common$data), 2)
  union <- run_reduce(d1, d2, return_common = FALSE)
  expect_setequal(names(union$data), c("a", "b", "c", "d"))
  expect_true(is.na(union$data$d[1]))
})
