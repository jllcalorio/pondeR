test_that("balance = 'all' gives equal counts in every match x group cell", {
  set.seed(7)
  d <- data.frame(
    `Age Group` = sample(c("2-5", "6-10"), 300, TRUE, prob = c(.6, .4)),
    Sex = sample(c("Male", "Female"), 300, TRUE),
    Status = sample(c("Stunted", "Non-stunted"), 300, TRUE, prob = c(.3, .7)),
    check.names = FALSE
  )
  r <- run_subset(d, match = c("Age Group", "Sex"), group = "Status",
                  independent = FALSE, balance = "all")
  tb <- with(r$chosen, table(`Age Group`, Sex, Status))
  expect_length(unique(as.vector(tb)), 1)
  expect_gt(sum(tb), 0)

  r2 <- run_subset(d, n = 40, match = c("Age Group", "Sex"), group = "Status",
                   independent = FALSE, balance = "all")
  tb2 <- with(r2$chosen, table(`Age Group`, Sex, Status))
  expect_length(unique(as.vector(tb2)), 1)
  expect_lte(nrow(r2$chosen), 40)
})
