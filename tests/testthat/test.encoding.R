context("Encoding functions")


test_that("compute_Uxij works with a simple basis of 1 function", {
  set.seed(42)
  
  K <- 4
  Tmax <- 10
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  m <- 1
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1)# base d'une seule fonction avec fonction constante = 1 entre 0 et Tmax
  I <- diag(rep(1, m))
  phi <- fd(I, b) # fonction constante = 1 entre 0 et Tmax
  
  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Uxij(x, phi, K)
  expectedOut <- rep(0, K)
  for(i in 1:K)
  {
    idx <- which(x$state == i)
    expectedOut[i] = sum(x[idx+1,"time"]-x[idx,"time"], na.rm = TRUE)
  }

  expect_length(out, K*m*m)
  expect_equal(out, expectedOut)
})


test_that("compute_Uxij works with a simple basis of 2 functions", {
  set.seed(42)
  
  K <- 4
  Tmax <- 10
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  m <- 2
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1)# base de deux fonctions:  constante = 1 entre 0 et Tmax/2 puis 0 et réciproquement
  I <- diag(rep(1, m))
  phi <- fd(I, b) 
  
  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Uxij(x, phi, K)
  expectedOut <- rep(0, K*m*m)
  for(i in 1:K)
  {
    idx <- which(x$state == i)
    idx1 <- idx[x[idx,"time"] <= 5]
    expectedOut[1 + (i-1)*m*m] = sum(pmin(x[idx1+1,"time"], 5) - x[idx1,"time"], na.rm = TRUE)
    
    idx2 <- idx[x[idx+1,"time"] > 5]
    expectedOut[4 + (i-1)*m*m] = sum(x[idx2+1,"time"] - pmax(x[idx2,"time"], 5), na.rm = TRUE)
  }
  
  expect_length(out, K*m*m)
  expect_lte(max(abs(out - expectedOut)), 1e-5)
})


test_that("compute_Vxi works with a simple basis of 1 function", {
  set.seed(42)
  
  K <- 4
  Tmax <- 10
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  m <- 1
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1)# base d'une seule fonction avec fonction constante = 1 entre 0 et Tmax
  I <- diag(rep(1, m))
  phi <- fd(I, b) # fonction constante = 1 entre 0 et Tmax
  
  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Vxi(x, phi, K)
  expectedOut <- rep(0, K)
  for(i in 1:K)
  {
    idx <- which(x$state == i)
    expectedOut[i] = sum(x[idx+1,"time"]-x[idx,"time"], na.rm = TRUE)
  }
  
  
  expect_length(out, K*m)
  expect_equal(out, expectedOut)
})



test_that("compute_Vxi works with a simple basis of 2 functions", {
  set.seed(42)
  
  K <- 4
  Tmax <- 10
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  m <- 2
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1)# base de deux fonctions:  constante = 1 entre 0 et Tmax/2 puis 0 et réciproquement
  I <- diag(rep(1, m))
  phi <- fd(I, b) 
  
  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Vxi(x, phi, K)
  expectedOut <- rep(0, K*m)
  for(i in 1:K)
  {
    idx <- which(x$state == i)
    idx1 <- idx[x[idx,"time"] <= 5]
    expectedOut[1 + (i-1)*m] = sum(pmin(x[idx1+1,"time"], 5) - x[idx1,"time"], na.rm = TRUE)
    
    idx2 <- idx[x[idx+1,"time"] > 5]
    expectedOut[2 + (i-1)*m] = sum(x[idx2+1,"time"] - pmax(x[idx2,"time"], 5), na.rm = TRUE)
  }
  
  expect_length(out, K*m)
  expect_lte(max(abs(out - expectedOut)), 1e-5)
})


test_that("compute_optimal_encoding works", {
  set.seed(42)
  n <- 200
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- msm2msmTmax(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(dT, b)
  
  expect_type(fmca, "list")
  expect_named(fmca, c("vp", "alpha", "pc", "F", "G", "V"))
  
  # eigenvaleus
  expect_length(fmca$vp, 2*m)
  trueEigVal <- 1/((1:m) * (2:(m+1)))
  expect_lte(max(abs(fmca$vp[1:m] - trueEigVal)), 0.01)
  
})
