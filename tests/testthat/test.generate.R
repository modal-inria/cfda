# Author: Quentin Grimonprez

context("Data generation")


test_that("generate_Markov_cfd output has the right format", {
  n <- 10
  K <- 4
  Tmax <- 10
  lambda_QJK <- c(1, 1, 1, 1)
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  
  d_JK <- generate_Markov_cfd(n = n, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  
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
