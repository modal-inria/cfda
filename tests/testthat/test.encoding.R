# Author: Quentin Grimonprez

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
    expectedOut[i] = sum(x$time[idx+1]-x$time[idx], na.rm = TRUE)
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
    idx1 <- idx[x$time[idx] <= 5]
    expectedOut[1 + (i-1)*m*m] = sum(pmin(x$time[idx1+1], 5) - x$time[idx1], na.rm = TRUE)
    
    idx2 <- idx[x$time[idx+1] > 5]
    expectedOut[4 + (i-1)*m*m] = sum(x$time[idx2+1] - pmax(x$time[idx2], 5), na.rm = TRUE)
  }
  
  expect_length(out, K*m*m)
  expect_lte(max(abs(out - expectedOut)), 1e-5)
})

oldcompute_Uxij <- function(x, phi, K)
{
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)
  
  for(state in 1:K) 
  {
    idx <- which(x$state == state)
    for(u in idx)
    {
      for(i in 1:nBasis) 
      {
        for(j in 1:nBasis)  
        {
          if(u < nrow(x))
          {
            aux[(state-1)*nBasis*nBasis + (i-1)*nBasis + j] = aux[(state-1)*nBasis*nBasis + (i-1)*nBasis + j] + 
              integrate(function(t) {
                eval.fd(t, phi[i]) * eval.fd(t, phi[j])
              }, lower = x$time[u], upper = x$time[u+1],
              stop.on.error = FALSE)$value
          }
        }
      }
    }
    
  }
  
  return(aux)     
}

test_that("refactor of compute_Uxij keeps the same results", {
  set.seed(42)
  
  K <- 4
  Tmax <- 10
  QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
  lambda_QJK <- c(1, 1, 1, 1)
  d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
  d_JK2 <- msm2msmTmax(d_JK, 10)
  
  
  m <- 10
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  I <- diag(rep(1, m))
  phi <- fd(I, b) 
  
  expectedOut <- by(d_JK2, d_JK2$id, function(x){oldcompute_Uxij(x, phi, K)})
  out <- by(d_JK2, d_JK2$id, function(x){compute_Uxij(x, phi, K)})
  
  expect_equal(out[[1]], expectedOut[[1]])
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
    expectedOut[i] = sum(x$time[idx+1]-x$time[idx], na.rm = TRUE)
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
    idx1 <- idx[x$time[idx] <= 5]
    expectedOut[1 + (i-1)*m] = sum(pmin(x$time[idx1+1], 5) - x$time[idx1], na.rm = TRUE)
    
    idx2 <- idx[x$time[idx+1] > 5]
    expectedOut[2 + (i-1)*m] = sum(x$time[idx2+1] - pmax(x$time[idx2], 5), na.rm = TRUE)
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
  expect_silent(fmca <- compute_optimal_encoding(dT, b, nCores = 1, verbose = FALSE))
  
  expect_type(fmca, "list")
  expect_named(fmca, c("eigenvalues", "alpha", "pc", "F", "G", "V", "basisobj"))
  
  # eigenvalues
  expect_length(fmca$eigenvalues, K*m)
  trueEigVal <- 1/((1:m) * (2:(m+1)))
  expect_lte(max(abs(fmca$eigenvalues[1:m] - trueEigVal)), 0.01)
  
  # alpha
  expect_type(fmca$alpha, "list")
  expect_length(fmca$alpha, m * K)
  expect_equal(dim(fmca$alpha[[1]]), c(m, K))
  
  # pc
  expect_equal(dim(fmca$pc), c(n, m * K))
  
  # F
  expect_equal(dim(fmca$F), c(2*m, 2*m))
  
  # G
  expect_equal(dim(fmca$G), c(2*m, 2*m))
  
  # V
  expect_equal(dim(fmca$V), c(n, 2*m))
  
  # basisobj
  expect_equal(fmca$basisobj, b)
})

test_that("compute_optimal_encoding works verbose", {
  set.seed(42)
  K <- 2
  d_JK <- generate_2State(n = 10)
  d_JK2 <- msm2msmTmax(d_JK, 1)
  
  # create basis object
  m <- 10
  b <- create.bspline.basis(c(0, 1), nbasis = m, norder = 4)
  
  # compute encoding
  expect_output(encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1, verbose = TRUE))
  
})


test_that("getEncoding works", {
  n <- 50
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- msm2msmTmax(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(dT, b, nCores = 1)
  
  out <- getEncoding(fmca, fdObject = TRUE)
  expect_s3_class(out, "fd")
  
  out <- getEncoding(fmca, harm = 3, fdObject = TRUE)
  expect_s3_class(out, "fd")
  
  out <- getEncoding(fmca, fdObject = FALSE)
  expect_named(out, c("x", "y"))
  expect_equal(dim(out$y), c(501, 2))
  expect_length(out$x, 501)
  
  out <- getEncoding(fmca, harm = 3, fdObject = FALSE)
  expect_named(out, c("x", "y"))
  expect_equal(dim(out$y), c(501, 2))
  expect_length(out$x, 501)
  
  out <- getEncoding(fmca, fdObject = FALSE, nx = 100)
  expect_named(out, c("x", "y"))
  expect_equal(dim(out$y), c(100, 2))
  expect_length(out$x, 100)
  
})

test_that("compute_optimal_encoding throws an error when the basis is not well suited", {
  
  data_msm <- data.frame(id = rep(1:2, each = 3), time = c(0, 3, 5, 0, 4, 5), state = c(1, 2, 2, 1, 2, 2))
  b <- create.bspline.basis(c(0, 5), nbasis = 3, norder = 2)
  
  expect_error({fmca <- compute_optimal_encoding(data_msm, b, nCores = 1)}, 
               regexp = "In the support of each basis function, each state must be present at least once (p(x_t) != 0 for t in the support).", 
               fixed = TRUE)
})


test_that("plot.fmca does not produce warnings", {
  n <- 50
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- msm2msmTmax(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(dT, b, nCores = 1)
  
  expect_warning(plot(fmca), regexp = NA)
  expect_warning(plot(fmca, harm = 3, col = c("red", "blue")), regexp = NA)
})


test_that("plotComponent does not produce warnings", {
  n <- 50
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- msm2msmTmax(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(dT, b, nCores = 1)
  
  expect_warning(plotComponent(fmca, addNames = TRUE), regexp = NA)
  expect_warning(plotComponent(fmca, addNames = FALSE), regexp = NA)
})


test_that("plotEigenvalues does not produce warnings", {
  n <- 50
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- msm2msmTmax(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  fmca <- compute_optimal_encoding(dT, b, nCores = 1)
  
  expect_warning(plotEigenvalues(fmca, cumulative = FALSE, normalize = FALSE), regexp = NA)
  expect_warning(plotEigenvalues(fmca, cumulative = FALSE, normalize = TRUE), regexp = NA)
  expect_warning(plotEigenvalues(fmca, cumulative = TRUE, normalize = FALSE), regexp = NA)
  expect_warning(plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE), regexp = NA)
})
