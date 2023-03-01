context("plotMD")
library(modi)

test_that("Plot of Mahalanobis distance", {
  data(bushfirem, bushfire.weights)
  det.res <- TRC(bushfirem, weights = bushfire.weights)
  PlotMD(det.res$dist, ncol(bushfirem))
})
