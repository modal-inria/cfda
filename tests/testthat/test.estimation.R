# Author: Quentin Grimonprez

context("Markov estimation")

set.seed(42)

test_that("completeStatetable does not change a square statetable", {
  
  aux <- matrix(c(0, 2, 3, 4, 1, 0, 3, 4, 1, 2, 0, 4, 1, 2, 3, 0), nrow = 4, dimnames = list(1:4, 1:4))
  
  out <- completeStatetable(aux)
  
  expect_equal(out, aux)
})


test_that("completeStatetable completes one missing row", {
  
  # middle
  aux <- matrix(c(0, 3, 4, 1, 3, 4, 1, 0, 4, 1, 3, 0), nrow = 3, dimnames = list(c(1, 3, 4), 1:4))
  
  out <- completeStatetable(aux)
  
  expectedOut <- rbind(aux[1,], rep(0, 4), aux[2:3,])
  rownames(expectedOut) = 1:4
  
  expect_equal(out, expectedOut)
  
  
  # first
  aux <- matrix(c(2, 3, 4, 0, 3, 4, 2, 0, 4, 2, 3, 0), nrow = 3, dimnames = list(c(2, 3, 4), 1:4))
  
  out <- completeStatetable(aux)
  
  expectedOut <- rbind(rep(0, 4), aux)
  rownames(expectedOut) = 1:4
  
  expect_equal(out, expectedOut)
  
  
  # last
  aux <- matrix(c(0, 2, 3, 1, 0, 3, 1, 2, 0, 1, 2, 3), nrow = 3, dimnames = list(c(1, 2, 3), 1:4))
  
  out <- completeStatetable(aux)
  
  expectedOut <- rbind(aux, rep(0, 4))
  rownames(expectedOut) = 1:4
  
  expect_equal(out, expectedOut)
  
})


test_that("completeStatetable completes several missing rows", {
  
  aux <- matrix(c(0, 3, 1, 3, 1, 0, 1, 3), nrow = 2, dimnames = list(c(1, 3), 1:4))
  
  out <- completeStatetable(aux)
  
  expectedOut <- rbind(aux[1,], rep(0, 4), aux[2,], rep(0, 4))
  rownames(expectedOut) = 1:4
  
  expect_equal(out, expectedOut)
  
})


test_that("estimate_Markov estimates well", {
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK = generate_Markov_cfd(n = 500, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 30)
  
  mark <- estimate_Markov(d_JK)
  
  expect_lte(sqrt(mean((mark$lambda-lambda_QJK)^2)), 0.06)
  expect_lte(sqrt(mean((mark$Q-QJK)^2)), 0.02)
})


test_that("plot_Markov does not produce warnings", {
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K, dimnames = list(1:4, 1:4)) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)

  dat <- list(Q = QJK, lambda = lambda_QJK)
  class(dat) = "Markov"
  
  expect_warning(plot(dat), regexp = NA)
})
