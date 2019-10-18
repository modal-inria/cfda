context("Markov estimation")


test_that("estimation_Markov estimates well",{
  K <- 4
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK = generate_Markov_cfd(n = 500, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 30)
  
  mark <- estimation_Markov(d_JK)

  expect_lte(sqrt(mean((mark$lambda-lambda_QJK)^2)), 0.05)
  expect_lte(sqrt(mean((mark$Q-QJK)^2)), 0.02)
})
