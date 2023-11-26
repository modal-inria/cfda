# Author: Quentin Grimonprez

context("Multivariate Encoding")


## common dataset
set.seed(42)

K <- 4
Tmax <- 10
PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
lambda_PJK <- c(1, 1, 1, 1)
d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
d_JK2 <- cut_data(d_JK, 10)
##

computeUmean_old <- function(data, phi, K, stateColumns, verbose, ...) {
  nId <- length(unique(data$id))

  if (verbose) {
    cat("\n")
    cat("---- Compute U matrix\n")
    pb <- timerProgressBar(min = 0, max = nId, width = 50)
    on.exit(close(pb))
    jj <- 1
  }
  list_of_Uval <- list()

  for (i in unique(data$id)) {
    list_of_Uval[[i]] <- compute_Uxij_multi(data[data$id == i, ], phi, K, stateColumns, verbose = FALSE)
    if (verbose) {
      setTimerProgressBar(pb, jj)
      jj <- jj + 1
    }
  }

  U_list_matrix <- compute_U_list_matrix(list_of_Uval, K)
  U_mean <- U_list_matrix
  for (i in seq_along(K)) {
    for (j in seq_along(K)) {
      U_mean[[i]][[j]] <- colMeans(U_mean[[i]][[j]])
    }
  }

  return(U_mean)
}

computeVlist_old <- function(data, phi, K, stateColumns, verbose, ...) {
  nId <- length(unique(data$id))

  # Compute V
  if (verbose) {
    cat("\n---- Compute V matrix\n")
    pb <- timerProgressBar(min = 0, max = nId, width = 50)
    on.exit(close(pb))
    jj <- 1
  }
  V_multi <- list()

  for (i in unique(data$id)) {
    V_multi[[i]] <- compute_Vxi_multi(data[data$id == i, ], phi, K, stateColumns = stateColumns, ...)
    if (verbose) {
      setTimerProgressBar(pb, jj)
      jj <- jj + 1
    }
  }
  V_multi <- do.call(rbind, V_multi)

  return(V_multi)
}


test_that("computeUmean keeps the same result than computeUmean_old", {
  skip_on_cran()
  set.seed(42)
  Tmax <- 2
  x1 <- generate_Markov(n = 50, K = 2)
  x1 <- cut_data(x1, Tmax = Tmax)
  x2 <- generate_Markov(n = 50, K = 2)
  x2 <- cut_data(x2, Tmax = Tmax)

  x <- list(x1, x2)
  x <- convert2mvcfd(x)
  basisobj <- create.bspline.basis(c(0, Tmax), nbasis = 10, norder = 4)
  nBasis <- basisobj$nbasis
  phi <- fd(diag(nBasis), basisobj)

  oldRes <- computeUmean_old(x, phi, K = c(2, 2), stateColumns = c("state1", "state2"), verbose = FALSE)
  newRes <- computeUmean(x, phi, K = c(2, 2), stateColumns = c("state1", "state2"), verbose = FALSE, nCores = 2)

  expect_equivalent(newRes, oldRes)
})


test_that("computeVlist keeps the same result than computeVlist_old", {
  skip_on_cran()
  set.seed(42)
  Tmax <- 2
  x1 <- generate_Markov(n = 50, K = 2)
  x1 <- cut_data(x1, Tmax = Tmax)
  x2 <- generate_Markov(n = 50, K = 2)
  x2 <- cut_data(x2, Tmax = Tmax)

  x <- list(x1, x2)
  x <- convert2mvcfd(x)
  basisobj <- create.bspline.basis(c(0, Tmax), nbasis = 10, norder = 4)
  nBasis <- basisobj$nbasis
  phi <- fd(diag(nBasis), basisobj)

  oldRes <- computeVlist_old(x, phi, K = c(2, 2), stateColumns = c("state1", "state2"), verbose = FALSE)
  newRes <- computeVlist(x, phi, K = c(2, 2), stateColumns = c("state1", "state2"), verbose = FALSE, nCores = 2)

  expect_equivalent(newRes, oldRes)
})

test_that("compute_optimal_encoding_multivariate throws error", {
  # create basis object
  m <- 10
  b <- create.bspline.basis(c(0, 1), nbasis = m, norder = 4)

  expect_error(
    compute_optimal_encoding_multivariate(m, basisobj = b, nCores = 1, verbose = TRUE),
    regexp = "data must be a data.frame."
  )

  expect_error(
    compute_optimal_encoding_multivariate(d_JK2, 3, nCores = 1, verbose = TRUE),
    regexp = "basisobj is not a basis object."
  )

  expect_error(
    compute_optimal_encoding_multivariate(d_JK2, b, nCores = 0, verbose = TRUE),
    regexp = "nCores must be an integer > 0."
  )
  expect_error(
    compute_optimal_encoding_multivariate(d_JK2, b, nCores = 2.5, verbose = TRUE),
    regexp = "nCores must be an integer > 0."
  )
  expect_error(
    compute_optimal_encoding_multivariate(d_JK2, b, nCores = NA, verbose = TRUE),
    regexp = "nCores must be an integer > 0."
  )
  expect_error(
    compute_optimal_encoding_multivariate(d_JK2, b, nCores = NaN, verbose = TRUE),
    regexp = "nCores must be an integer > 0."
  )
  expect_error(
    compute_optimal_encoding_multivariate(d_JK2, b, nCores = 1, verbose = 2),
    regexp = "verbose must be either TRUE or FALSE."
  )
})
