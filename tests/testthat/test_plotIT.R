context("plotIT")
library(modi)

test_that("Plot of infection times of Epidemic Algorithm", {
  it <- c(rep(NA, 3), rep(1:7, times=c(1, 4, 10, 8, 5, 3, 2)))
  wt <- rep(c(1,2,5), times=12)
  pout <- plotIT(it, wt, 6)
  expect_null(pout)
})
