# Author: Quentin Grimonprez

context("Predict function")

skip_on_cran()

## common dataset
set.seed(42)

n <- 25
Tmax <- 1
K <- 2
m <- 6
d <- generate_2State(n)
dT <- cut_data(d, Tmax)
row.names(dT) = NULL
##

test_that("predict throws error", {
  fmca2 <- list()
  class(fmca2) = "fmca"
  
  expect_error(predict(dT, dT, nCores = 1, verbose = TRUE))
  
  expect_error(predict(fmca2, data.frame(), nCores = 1, verbose = TRUE), regexp = "Missing columns in data: id, time, state.")
  
  expect_error(predict(fmca2, dT, nCores = 0, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  expect_error(predict(fmca2, dT, nCores = 2.5, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  expect_error(predict(fmca2, dT, nCores = NA, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  expect_error(predict(fmca2, dT, nCores = NaN, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  
  expect_error(predict(fmca2, dT, nCores = 1, verbose = 2), regexp = "verbose must be either TRUE or FALSE.")
})


test_that("predict works", {
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(dT, b, computeCI = FALSE, nCores = 1, verbose = FALSE)

  expect_silent(out <- predict(fmca, dT, nCores = 1, verbose = FALSE))
  
  expect_equivalent(out, fmca$pc)
})

