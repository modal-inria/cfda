context("Data generation")


test_that("run cluster/predict R object",{
  n <- 10
  K <- 4
  lambda_QJK <- c(1, 1, 1, 1)
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  
  d_JK = generate_Markov_cfd (n = n, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
  
  expect_equal(ncol(d_JK), 3)
  expect_equal(colnames(d_JK), c("id", "time", "state"))
  expect_equal(sort(unique(d_JK$id)), 1:10)
  expect_equal(sort(unique(d_JK$state)), 1:4)
  expect_true(all(d_JK$time <= 10))
  expect_true(all(d_JK$time >= 0))
  expect_equal(sum(d_JK$time == 0), n)
})
