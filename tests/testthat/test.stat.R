# Author: Quentin Grimonprez

context("Statistics on data object")

test_that("compute_time_spent_intern works", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))
  
  out <- compute_time_spent_intern(dat, 1:3)
  expectedOut <- c(1.5 + 1, 1, 2.5)
  
  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})


test_that("compute_time_spent_intern works with more labels", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))
  
  out <- compute_time_spent_intern(dat, 1:4)
  expectedOut <- c(1.5 + 1, 1, 2.5, 0)
  
  expect_equal(sum(out), max(dat$time))
  expect_equal(out, expectedOut)
})

test_that("compute_time_spent works", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  
  out <- compute_time_spent(dat)
  expectedOut <- rbind(c(1.5 + 1, 1, 2.5),
                       c(3, 3, 0))
  colnames(expectedOut) = 1:3
  rownames(expectedOut) = 1:2
  class(expectedOut) = "timeSpent"
  
  expect_equal(out, expectedOut)
})

test_that("boxplot.timeSpent does not produce warnings", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 6, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  
  out <- compute_time_spent(dat)

  expect_warning(boxplot(out), regexp = NA)
  expect_warning(boxplot(out, col = c("red", "blue", "green"), outlier.colour = "black"), regexp = NA)
  
})


test_that("compute_duration works", {
  dat <- data.frame(id = rep(1:2, c(5, 3)), time = c(0, 1.5, 4, 5, 7, 0, 3, 6), state = c(1, 3, 2, 1, 1, 1, 2, 3))
  
  out <- compute_duration(dat)
  expectedOut <- c("1" = 7, "2" = 6)
  class(expectedOut) = "duration"
  
  expect_equivalent(out, expectedOut)
})


test_that("hist.duration does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  duration <- compute_duration(d_JK)
  
  expect_warning(hist(duration), regexp = NA)
  expect_warning(hist(duration, color = "red"), regexp = NA)
})


test_that("id_get_state returns the right state", {
  dat <- data.frame(id = rep(1, 5), time = c(0, 1.5, 4, 5, 6), state = c(1, 3, 2, 1, 1))

  out <- id_get_state(dat, 3, NAafterTmax = FALSE) 
  expect_equal(out, 3)
  
  out <- id_get_state(dat, 7, NAafterTmax = TRUE) 
  expect_equal(out, NA)
  
  out <- id_get_state(dat, 7, NAafterTmax = FALSE) 
  expect_equal(out, 1)
  
  out <- id_get_state(dat, 6, NAafterTmax = TRUE) 
  expect_equal(out, 1)
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
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 6), state = c(1, 3, 2, 1, 1, 1, 
                                                                                        2, 3, 1, 2, 2))
  out <- estimate_pt(dat) 
  
  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, c(0, 1, 1.5, 2, 3, 3.5, 4, 5, 6))
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1/2, 1/2, 0, 0, 1/2, 1/2, 0, 0, 1, 1/2, 1/2, 0, 1, 0, 0, 1/2, 1/2, 0, 1/2, 1/2, 0, 1/2, 1/2, 0, 1/2, 1/2, 0), nrow = 3, dimnames = list(1:3, out$t)))
  
  
  out <- estimate_pt(dat, NAafterTmax = TRUE) 
  
  expect_length(out, 2)
  expect_equal(names(out), c("pt", "t"))
  expect_equal(out$t, c(0, 1, 1.5, 2, 3, 3.5, 4, 5, 6))
  expect_equivalent(colSums(out$pt), rep(1, ncol(out$pt)))
  expect_equal(out$pt, matrix(c(1/2, 1/2, 0, 0, 1/2, 1/2, 0, 0, 1, 1/2, 1/2, 0, 1, 0, 0, 1/2, 1/2, 0, 1/2, 1/2, 0, 1/2, 1/2, 0, 0, 1, 0), nrow = 3, dimnames = list(1:3, out$t)))
})



test_that("plot_pt_classic does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  d_JK2 <- cut_data(d_JK, 10)
  
  pt <- estimate_pt(d_JK2)
  
  expect_warning(plot_pt_classic(pt), regexp = NA)
})


test_that("plot_pt_ribbon does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  d_JK2 <- cut_data(d_JK, 10)
  
  pt <- estimate_pt(d_JK2)
  
  expect_warning(plot_pt_ribbon(pt), regexp = NA)
})

test_that("plot_pt does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)

  d_JK2 <- cut_data(d_JK, 10)

  pt <- estimate_pt(d_JK2)

  expect_warning(plot(pt, ribbon = FALSE), regexp = NA)
  expect_warning(plot(pt, ribbon = FALSE, col = c("red", "blue", "green", "black")), regexp = NA)
  expect_warning(plot(pt, ribbon = TRUE, addBorder = TRUE), regexp = NA)
  expect_warning(plot(pt, ribbon = TRUE, addBorder = FALSE), regexp = NA)
  expect_warning(plot(pt, ribbon = TRUE, addBorder = FALSE, col = c("red", "blue", "green", "black")), regexp = NA)
})


test_that("compute_number_jumps works", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1:6, 1:5))
  out <- compute_number_jumps(dat) 
  expectedOut <- c(5, 4)
  class(expectedOut) = "njump"
  
  expect_equivalent(out, expectedOut)
})

test_that("compute_number_jumpsIntern works with countDuplicated = TRUE", {
  dat <- data.frame(id = 1:20, time = 1:20, state = rep(1:10, each = 2))
  out <- compute_number_jumpsIntern(dat, countDuplicated = TRUE) 
  expectedOut <- 19

  expect_equivalent(out, expectedOut)
  
  
  dat <- data.frame(id = 1:20, time = 1:20, state = rep(letters[1:10], each = 2))
  out <- compute_number_jumpsIntern(dat, countDuplicated = TRUE) 
  expectedOut <- 19
  
  expect_equivalent(out, expectedOut)
})

test_that("compute_number_jumpsIntern works with countDuplicated = FALSE", {
  # without duplicate state
  dat <- data.frame(id = 1:20, time = 1:20, state = 1:20)
  out <- compute_number_jumpsIntern(dat, countDuplicated = FALSE) 
  expectedOut <- 19
  
  expect_equivalent(out, expectedOut)
  
  # with ordered time
  dat <- data.frame(id = 1:20, time = 1:20, state = rep(1:10, each = 2))
  out <- compute_number_jumpsIntern(dat, countDuplicated = FALSE) 
  expectedOut <- 9

  expect_equivalent(out, expectedOut)
  
  # with unordered time
  dat <- data.frame(id = 1:20, time = c(11:20, 1:10), state = as.factor(rep(letters[1:5], each = 4)))
  out <- compute_number_jumpsIntern(dat, countDuplicated = FALSE) 
  expectedOut <- 5
  
  expect_equivalent(out, expectedOut)
})



test_that("hist.njump does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  njump <- compute_number_jumps(d_JK)
  
  expect_warning(hist(njump), regexp = NA)
  expect_warning(hist(njump, color = "red"), regexp = NA)
})


test_that("statetable works", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1:6, 1:5))
  out <- statetable(dat) 
  
  expectedOut <- matrix(0, nrow = 5, ncol = 6)
  expectedOut[1, 2] = expectedOut[2, 3] = expectedOut[3, 4] = expectedOut[4, 5] = 2
  expectedOut[5, 6] = 1

  expect_equivalent(out, as.table(expectedOut))
})


test_that("rep_large_ind works", {
  dat <- data.frame(id = rep(1:2, c(6, 5)), time = c(0:5, 0, 1.5, 2, 3.5, 5), state = c(1:6, 1:5))
  out <- rep_large_ind(dat)
  expectedOut <- data.frame(id = rep(1:2, c(5, 4)), 
                            t_start = c(0:4, 0, 1.5, 2, 3.5), 
                            t_end = c(1:5, 1.5, 2, 3.5, 5),
                            state = c(1:5, 1:4))
  
  expect_equivalent(out, expectedOut)
  
  
  
  out <- rep_large_ind(dat[dat$id == 1, ])
  
  expect_equivalent(out, expectedOut[expectedOut$id == 1, ])
})



test_that("plotData does not produce warnings", {
  # simulate the Jukes Cantor models of nucleotides replacement.
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  d_JKT <- cut_data(d_JK, Tmax = 10)

  expect_warning(plotData(d_JK, addId = TRUE, addBorder = TRUE), regexp = NA)
  expect_warning(plotData(d_JK, addId = FALSE, addBorder = FALSE, col = c("red", "blue", "green", "yellow")), regexp = NA)
  
})
