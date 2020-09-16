# Author: Quentin Grimonprez

context("Bootstrap functions")

test_that("getSignReference works", {
  alpha <- list(matrix(1:9, nrow = 3), matrix(-1 * (9:1), nrow = 3), matrix(c(1:5, 4:1), nrow = 3))
  expectedOut <- list(position = c(9, 1, 5), isNegative = c(FALSE, TRUE, FALSE))
  
  out <- getSignReference(alpha)
  
  expect_equal(out, expectedOut)
})


test_that("unifySign works when all are present", {
  ref <- list(position = c(9, 1, 5), isNegative = c(FALSE, TRUE, FALSE))
  
  encod <- list(list(alpha = list(matrix(1:9, nrow = 3), matrix(1:9, nrow = 3), matrix(1:9, nrow = 3)), 
                     pc = matrix(1:9, nrow = 3), 
                     invF05vec = matrix(1:9, nrow = 3)),
                list(alpha = list(matrix(-1 * (9:1), nrow = 3), matrix(-1 * (9:1), nrow = 3), matrix(-1 * (9:1), nrow = 3)), 
                     pc = matrix(-1 * (9:1), nrow = 3), 
                     invF05vec = matrix(-1 * (9:1), nrow = 3)))
  
  expectedOut <- list(list(alpha = list(matrix(1:9, nrow = 3), matrix(-1*(1:9), nrow = 3), matrix(1:9, nrow = 3)), 
                           pc = matrix(c(1:3, -4, -5, -6, 7:9), nrow = 3), 
                           invF05vec = matrix(c(1:3, -4, -5, -6, 7:9), nrow = 3)),
                      list(alpha = list(matrix((9:1), nrow = 3), matrix(-1 * (9:1), nrow = 3), matrix((9:1), nrow = 3)), 
                           pc = matrix(c(9:7, -6, -5, -4, 3:1), nrow = 3), 
                           invF05vec = matrix(c(9:7, -6, -5, -4, 3:1), nrow = 3)))
  
  out <- unifySign(encod, ref)
  
  expect_equal(out, expectedOut)
})


test_that("unifySign works when there are some NULL elements", {
  ref <- list(position = c(9, 1, 5), isNegative = c(FALSE, TRUE, FALSE))
  
  encod <- list(list(alpha = list(matrix(1:9, nrow = 3), matrix(1:9, nrow = 3), matrix(1:9, nrow = 3)), 
                     pc = matrix(1:9, nrow = 3), 
                     invF05vec = matrix(1:9, nrow = 3)),
                NULL,
                list(alpha = list(matrix(-1 * (9:1), nrow = 3), matrix(-1 * (9:1), nrow = 3), matrix(-1 * (9:1), nrow = 3)), 
                     pc = matrix(-1 * (9:1), nrow = 3), 
                     invF05vec = matrix(-1 * (9:1), nrow = 3)))
  
  expectedOut <- list(list(alpha = list(matrix(1:9, nrow = 3), matrix(-1*(1:9), nrow = 3), matrix(1:9, nrow = 3)), 
                           pc = matrix(c(1:3, -4, -5, -6, 7:9), nrow = 3), 
                           invF05vec = matrix(c(1:3, -4, -5, -6, 7:9), nrow = 3)),
                      NULL,
                      list(alpha = list(matrix((9:1), nrow = 3), matrix(-1 * (9:1), nrow = 3), matrix((9:1), nrow = 3)), 
                           pc = matrix(c(9:7, -6, -5, -4, 3:1), nrow = 3), 
                           invF05vec = matrix(c(9:7, -6, -5, -4, 3:1), nrow = 3)))
  
  out <- unifySign(encod, ref)
  
  expect_equal(out, expectedOut)
})



test_that("compute_optimal_encoding throws error", {
  set.seed(42)
  K <- 2
  d_JK <- generate_2State(n = 10)
  d_JK2 <- cut_data(d_JK, 1)
  
  # create basis object
  m <- 10
  b <- create.bspline.basis(c(0, 1), nbasis = m, norder = 4)
  
  expect_error(compute_optimal_encoding(d_JK2, b, nCores = 1, computeCI = 2), regexp = "computeCI must be either TRUE or FALSE.")
  
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 0, propBootstrap = 0.5), regexp = "nBootstrap must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 10.5, propBootstrap = 0.5), regexp = "nBootstrap must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = NA, propBootstrap = 0.5), regexp = "nBootstrap must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = c(), propBootstrap = 0.5), regexp = "nBootstrap must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = c(20, 50), propBootstrap = 0.5), regexp = "nBootstrap must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = NaN, propBootstrap = 0.5), regexp = "nBootstrap must be an integer > 0.")
  
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = 0), regexp = "propBootstrap must be a real between 0 and 1.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = 1.5), regexp = "propBootstrap must be a real between 0 and 1.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = c(0.5, 0.8)), regexp = "propBootstrap must be a real between 0 and 1.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = NA), regexp = "propBootstrap must be a real between 0 and 1.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = NaN), regexp = "propBootstrap must be a real between 0 and 1.")
  expect_error(compute_optimal_encoding(d_JK2, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = c()), regexp = "propBootstrap must be a real between 0 and 1.")
})




test_that("compute_optimal_encoding works with computeCI = TRUE", {
  set.seed(42)
  n <- 200
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- cut_data(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  expect_silent(fmca <- compute_optimal_encoding(dT, b, computeCI = TRUE, nBootstrap = 50, propBootstrap = 0.5, nCores = 1, verbose = FALSE))
  
  expect_type(fmca, "list")
  expect_named(fmca, c("eigenvalues", "alpha", "pc", "F", "G", "invF05vec", "V", "basisobj", "bootstrap"))
  
  ## bootstrap
  expect_length(fmca$bootstrap, 50)
  expect_named(fmca$bootstrap[[1]], c("eigenvalues", "alpha", "pc", "F", "G", "invF05vec"))
  
  # eigenvalues
  expect_length(fmca$bootstrap[[1]]$eigenvalues, K*m)
  trueEigVal <- 1/((1:m) * (2:(m+1)))
  expect_lte(max(abs(fmca$eigenvalues[1:m] - trueEigVal)), 0.01)
  
  # alpha
  expect_type(fmca$bootstrap[[1]]$alpha, "list")
  expect_length(fmca$bootstrap[[1]]$alpha, m * K)
  expect_equal(dim(fmca$bootstrap[[1]]$alpha[[1]]), c(m, K))
  
  # pc
  expect_equal(dim(fmca$bootstrap[[1]]$pc), c(100, m * K))
  
  # F
  expect_equal(dim(fmca$bootstrap[[1]]$F), c(2*m, 2*m))
  
  # G
  expect_equal(dim(fmca$bootstrap[[1]]$G), c(2*m, 2*m))
  
  # invF05vec
  expect_equal(dim(fmca$bootstrap[[1]]$invF05vec), c(m * K, m * K))
})
