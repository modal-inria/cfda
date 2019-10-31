# Author: Quentin Grimonprez

context("Statistics on data_msm object")

test_that("compute_Time_Spent_intern works", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))
  
  out <- compute_Time_Spent_intern(dat, K = 3)
  expectedOut <- c(1.5 + 1, 1, 2.5)
  
  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})


test_that("compute_Time_Spent_intern works with a higher K", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))
  
  out <- compute_Time_Spent_intern(dat, K = 4)
  expectedOut <- c(1.5 + 1, 1, 2.5, 0)
  
  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})

test_that("compute_Time_Spent works", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  
  out <- compute_Time_Spent(dat, K = 3)
  expectedOut <- rbind(c(1.5 + 1, 1, 2.5),
                       c(3, 3, 0))
  colnames(expectedOut) = 1:3
  rownames(expectedOut) = 1:2
  class(expectedOut) = "timeSpent"
  
  expect_equal(out, expectedOut)
})

test_that("boxplot.timeSpent does not produce warnings", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  
  out <- compute_Time_Spent(dat, K = 3)

  expect_warning(boxplot(out), regexp = NA)
  
})


test_that("id_get_state returns the right state", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))

  out <- id_get_state(dat, 3) 
    
  expect_equal(out, 3)
})


test_that("get_state returns the right state", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  
  out <- get_state(dat, 3) 
  
  expect_equivalent(out, c(3, 2))
})


test_that("estimate_pt works with same t", {
  dat <- data.frame(id = rep(1:2, each = 6), time = rep(0:5, 2), state = c(1, 3, 2, 1, 1, 1, 
                                                                           2, 3, 1, 2, 3, 1))
  out <- estimate_pt(dat) 
  
  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, 0:5)
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1/2, 1/2, 0, 0, 0, 1, 1/2, 1/2, 0, 1/2, 1/2, 0, 1/2, 0, 1/2, 1, 0, 0), nrow = 3, dimnames = list(1:3, 0:5)))
})


test_that("estimate_pt works with different t", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1, 3, 2, 1, 1, 1, 
                                                                                         2, 3, 1, 2, 2))
  out <- estimate_pt(dat) 
  
  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, c(0, 1, 1.5, 2, 3, 3.5, 4, 5))
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1/2, 1/2, 0, 0, 1/2, 1/2, 0, 0, 1, 1/2, 1/2, 0, 1, 0, 0, 1/2, 1/2, 0, 1/2, 1/2, 0, 1/2, 1/2, 0), nrow = 3, dimnames = list(1:3, out$t)))
})



test_that("plot_pt_classic does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  pt <- estimate_pt(d_JK2)
  
  expect_warning(plot_pt_classic(pt), regexp = NA)
})


test_that("plot_pt_ribbon does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  pt <- estimate_pt(d_JK2)
  
  expect_warning(plot_pt_ribbon(pt), regexp = NA)
})

test_that("plot_pt does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)

  d_JK2 <- msm2msmTmax(d_JK, 10)

  pt <- estimate_pt(d_JK2)

  expect_warning(plot(pt, ribbon = FALSE), regexp = NA)
  expect_warning(plot(pt, ribbon = TRUE), regexp = NA)
})


test_that("compute_number_jumps works", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1:6, 1:5))
  out <- compute_number_jumps(dat) 
  expectedOut <- c(5, 4)
  class(expectedOut) = "njump"
  
  expect_equivalent(out, expectedOut)
})

