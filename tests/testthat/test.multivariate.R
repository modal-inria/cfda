# Author: Quentin Grimonprez

context("Multivariate Encoding")

test_that("convert2mvcfd works with 1 indiv and different time", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 0.9), state = c(1, 2, 3))
  x <- list(x1, x2)

  expectedOut <- data.frame(
    id = rep(1, 5), time = c(0, 0.25, 0.5, 0.9, 1), state1 = c(1, 1, 2, 2, 1), state2 = c(1, 2, 2, 3, 3)
  )
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd works with several ind", {
  x1 <- data.frame(id = c(1, 1, 1, 2, 2), time = c(0, 0.5, 1, 0, 1.5), state = c(1, 2, 1, 1, 2))
  x2 <- data.frame(id = c(1, 1, 1, 2, 2), time = c(0, 0.25, 0.9, 0, 2), state = c(1, 2, 3, 1, 2))
  x <- list(x1, x2)

  expectedOut <- data.frame(
    id = rep(c(1, 2), c(5, 3)),
    time = c(0, 0.25, 0.5, 0.9, 1, 0, 1.5, 2),
    state1 = c(1, 1, 2, 2, 1, 1, 2, 2),
    state2 = c(1, 2, 2, 3, 3, 1, 1, 2)
  )
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd works with 1 indiv and same time", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 1), state = c(1, 2, 3))
  x <- list(x1, x2)

  expectedOut <- data.frame(id = c(1, 1, 1, 1), time = c(0, 0.25, 0.5, 1), state1 = c(1, 1, 2, 1), state2 = c(1, 2, 2, 3))
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})

test_that("convert2mvcfd works with more than 2 dataframes", {
  x1 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 2, 1))
  x2 <- data.frame(id = c(1, 1, 1), time = c(0, 0.25, 1), state = c(1, 2, 3))
  x3 <- data.frame(id = c(1, 1, 1), time = c(0, 0.5, 1), state = c(1, 3, 1))
  x <- list(x1, x2, x3)

  expectedOut <- data.frame(
    id = c(1, 1, 1, 1), time = c(0, 0.25, 0.5, 1), state1 = c(1, 1, 2, 1), state2 = c(1, 2, 2, 3), state3 = c(1, 1, 3, 1)
  )
  out <- convert2mvcfd(x)
  expect_equal(out, expectedOut)
})


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
