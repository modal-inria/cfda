context("Data generation")


test_that("check generate_Markov_cfd output has the right format",{
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
  
  # all time are betwwen 0 and Tmax
  expect_true(all(d_JK$time <= Tmax))
  expect_true(all(d_JK$time >= 0))
  
  # each trajectory must start by a 0, so there must have n 0
  expect_equal(sum(d_JK$time == 0), n)
  
  # check time values are ordered per trajectory
  expect_true(all(tapply(d_JK$time, d_JK$id, function(x) {all(order(x) == seq_along(x))})))
})
