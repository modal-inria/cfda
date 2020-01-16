# Author: Quentin Grimonprez

context("Data generation")

test_that("generate_Markov returns error", {
  n <- 10
  K <- 4
  Tmax <- 10
  lambda_QJK <- c(1, 1, 1, 1)
  pi0 = c(1, rep(0, K-1))
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  
  expect_error(generate_Markov(n = NA, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "n must be a positive integer.")
  expect_error(generate_Markov(n = NaN, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "n must be a positive integer.")
  expect_error(generate_Markov(n = 3.5, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "n must be a positive integer.")
  expect_error(generate_Markov(n = c(3, 2), K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "n must be a positive integer.")
  
  expect_error(generate_Markov(n = n, K = NA, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "K must be an integer > 1.")
  expect_error(generate_Markov(n = n, K = NaN, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "K must be an integer > 1.")
  expect_error(generate_Markov(n = n, K = 3.5, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "K must be an integer > 1.")
  expect_error(generate_Markov(n = n, K = 0, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "K must be an integer > 1.")
  expect_error(generate_Markov(n = n, K = c(3, 2), Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax), "K must be an integer > 1.")
  
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = NA), "Tmax must be a positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = NaN), "Tmax must be a positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = -0.5), "Tmax must be a positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = c(5, 5.5)), "Tmax must be a positive real.")
  
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = rep(1, 3), pi0 = pi0, Tmax = Tmax), "lambda must be a vector of length K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = rep(-0.5, 4), pi0 = pi0, Tmax = Tmax), "lambda must be a vector of length K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = c(1, NA, 1, 1), pi0 = pi0, Tmax = Tmax), "lambda must be a vector of length K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = c(1, NaN, 1, 1), pi0 = pi0, Tmax = Tmax), "lambda must be a vector of length K of positive real.")
  
  expect_error(generate_Markov(n = n, K = K, Q = matrix(0, nrow = 3, ncol = 4), lambda = lambda, pi0 = pi0, Tmax = Tmax), "Q must be a matrix of size K x K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = matrix(-0.5, nrow = 4, ncol = 4), lambda = lambda, pi0 = pi0, Tmax = Tmax), "Q must be a matrix of size K x K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = matrix(NA, nrow = 4, ncol = 4), lambda = lambda, pi0 = pi0, Tmax = Tmax), "Q must be a matrix of size K x K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = matrix(NaN, nrow = 4, ncol = 4), lambda = lambda, pi0 = pi0, Tmax = Tmax), "Q must be a matrix of size K x K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = rep(4, 4), lambda = lambda, pi0 = pi0, Tmax = Tmax), "Q must be a matrix of size K x K of positive real.")
  
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = rep(2, 3), Tmax = Tmax), "pi0 must be a vector of length K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = rep(-0.5, 4), Tmax = Tmax), "pi0 must be a vector of length K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = rep(NA, 4), Tmax = Tmax), "pi0 must be a vector of length K of positive real.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = rep(NaN, 4), Tmax = Tmax), "pi0 must be a vector of length K of positive real.")
  
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax, labels = letters[1:3]), "labels must be NULL or a vector of length K.")
  expect_error(generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, pi0 = pi0, Tmax = Tmax, labels = rep("a", 4)), "labels must be NULL or a vector of length K.")
})


test_that("generate_Markov output has the right format", {
  n <- 10
  K <- 4
  Tmax <- 10
  lambda_QJK <- c(1, 1, 1, 1)
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  
  d_JK <- generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  
  ## format
  expect_s3_class(d_JK, "data.frame")
  expect_equal(ncol(d_JK), 3)
  expect_equal(colnames(d_JK), c("id", "time", "state"))
  
  
  ## content
  # id are between 1 and n
  expect_equal(sort(unique(d_JK$id)), 1:n)
  
  # state are betwwen 1 and K
  expect_equal(sort(unique(d_JK$state)), 1:K)
  
  # all time are between 0 and Tmax
  expect_true(all(d_JK$time <= Tmax))
  expect_true(all(d_JK$time >= 0))
  
  # each trajectory must start by a 0, so there must have n 0
  expect_equal(sum(d_JK$time == 0), n)
  
  # check time values are ordered per trajectory
  expect_true(all(tapply(d_JK$time, d_JK$id, function(x) {all(order(x) == seq_along(x))})))
})

test_that("generate_Markov output has the right format with labels", {
  n <- 10
  K <- 4
  Tmax <- 10
  lambda_QJK <- c(1, 1, 1, 1)
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  
  d_JK <- generate_Markov(n = n, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax, labels = LETTERS[1:4])
  
  ## format
  expect_s3_class(d_JK, "data.frame")
  expect_equal(ncol(d_JK), 3)
  expect_equal(colnames(d_JK), c("id", "time", "state"))
  
  
  ## content
  # id are between 1 and n
  expect_equal(sort(unique(d_JK$id)), 1:n)
  
  # state are betwwen 1 and K
  expect_equal(sort(unique(d_JK$state)), LETTERS[1:4])
  
  # all time are between 0 and Tmax
  expect_true(all(d_JK$time <= Tmax))
  expect_true(all(d_JK$time >= 0))
  
  # each trajectory must start by a 0, so there must have n 0
  expect_equal(sum(d_JK$time == 0), n)
  
  # check time values are ordered per trajectory
  expect_true(all(tapply(d_JK$time, d_JK$id, function(x) {all(order(x) == seq_along(x))})))
})

test_that("generate_2State returns error", {
  expect_error(generate_2State(3.5), "n must be a positive integer.")
  expect_error(generate_2State(c(3, 2)), "n must be a positive integer.")
  expect_error(generate_2State(NA), "n must be a positive integer.")
  expect_error(generate_2State(NaN), "n must be a positive integer.")
})


test_that("generate_2State output has the right format", {
  n <- 10

  d <- generate_2State(n)
  
  ## format
  expect_s3_class(d, "data.frame")
  expect_equal(ncol(d), 3)
  expect_equal(colnames(d), c("id", "time", "state"))
  
  
  ## content
  # id are between 1 and n
  expect_equal(sort(unique(d$id)), 1:n)
  
  # state are betwwen 1 and 2
  expect_equal(sort(unique(d$state)), 1:2)
  
  # all time are between 0 and 1
  expect_true(all(d$time <= 1))
  expect_true(all(d$time >= 0))
  
  # each trajectory must start by a 0, so there must have n 0
  expect_equal(sum(d$time == 0), n)
  
  # check time values are ordered per trajectory
  expect_true(all(tapply(d$time, d$id, function(x) {all(order(x) == seq_along(x))})))
  
  # check state values are ordered per trajectory
  expect_true(all(tapply(d$state, d$id, function(x) {all(order(x) == seq_along(x))})))
  
  # each individuals must have records
  expect_equivalent(as.numeric(table(d$id)), rep(2, n))
})
