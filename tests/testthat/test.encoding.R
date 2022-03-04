# Author: Quentin Grimonprez

context("Encoding functions")

## common dataset
set.seed(42)

K <- 4
Tmax <- 10
PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
lambda_PJK <- c(1, 1, 1, 1)
d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
d_JK2 <- cut_data(d_JK, 10)
##

test_that("compute_Uxij works with a simple basis of 1 function", {
  set.seed(42)

  m <- 1
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1) # base d'une seule fonction avec fonction constante = 1 entre 0 et Tmax
  I <- diag(rep(1, m))
  phi <- fd(I, b) # fonction constante = 1 entre 0 et Tmax

  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Uxij(x, phi, K)
  expectedOut <- rep(0, K)
  for (i in 1:K)
  {
    idx <- which(x$state == i)
    expectedOut[i] <- sum(x$time[idx + 1] - x$time[idx], na.rm = TRUE)
  }

  expect_length(out, K * m * m)
  expect_equal(out, expectedOut)
})


test_that("compute_Uxij works with a simple basis of 2 functions", {
  skip_on_cran()
  m <- 2
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1) # base de deux fonctions:  constante = 1 entre 0 et Tmax/2 puis 0 et réciproquement
  I <- diag(rep(1, m))
  phi <- fd(I, b)

  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Uxij(x, phi, K)
  expectedOut <- rep(0, K * m * m)
  for (i in 1:K)
  {
    idx <- which(x$state == i)
    idx1 <- idx[x$time[idx] <= 5]
    expectedOut[1 + (i - 1) * m * m] <- sum(pmin(x$time[idx1 + 1], 5) - x$time[idx1], na.rm = TRUE)

    idx2 <- idx[x$time[idx + 1] > 5]
    expectedOut[4 + (i - 1) * m * m] <- sum(x$time[idx2 + 1] - pmax(x$time[idx2], 5), na.rm = TRUE)
  }

  expect_length(out, K * m * m)
  expect_lte(max(abs(out - expectedOut)), 1e-5)
})

oldcompute_Uxij <- function(x, phi, K) {
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)

  for (state in 1:K)
  {
    idx <- which(x$state == state)
    for (u in idx)
    {
      for (i in 1:nBasis)
      {
        for (j in 1:nBasis)
        {
          if (u < nrow(x)) {
            aux[(state - 1) * nBasis * nBasis + (i - 1) * nBasis + j] <- aux[(state - 1) * nBasis * nBasis + (i - 1) * nBasis + j] +
              integrate(function(t) {
                eval.fd(t, phi[i]) * eval.fd(t, phi[j])
              },
              lower = x$time[u], upper = x$time[u + 1],
              stop.on.error = FALSE
              )$value
          }
        }
      }
    }
  }

  return(aux)
}

test_that("refactor of compute_Uxij keeps the same results", {
  skip_on_cran()
  skip_on_ci()
  skip_on_covr()

  m <- 10
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  I <- diag(rep(1, m))
  phi <- fd(I, b)

  expectedOut <- by(d_JK2, d_JK2$id, function(x) {
    oldcompute_Uxij(x, phi, K)
  })
  out <- by(d_JK2, d_JK2$id, function(x) {
    compute_Uxij(x, phi, K)
  })

  expect_equal(out[[1]], expectedOut[[1]])
})


test_that("compute_Vxi works with a simple basis of 1 function", {
  m <- 1
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1) # base d'une seule fonction avec fonction constante = 1 entre 0 et Tmax
  I <- diag(rep(1, m))
  phi <- fd(I, b) # fonction constante = 1 entre 0 et Tmax

  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Vxi(x, phi, K)
  expectedOut <- rep(0, K)
  for (i in 1:K)
  {
    idx <- which(x$state == i)
    expectedOut[i] <- sum(x$time[idx + 1] - x$time[idx], na.rm = TRUE)
  }


  expect_length(out, K * m)
  expect_equal(out, expectedOut)
})



test_that("compute_Vxi works with a simple basis of 2 functions", {
  skip_on_cran()
  m <- 2
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 1) # base de deux fonctions:  constante = 1 entre 0 et Tmax/2 puis 0 et réciproquement
  I <- diag(rep(1, m))
  phi <- fd(I, b)

  x <- d_JK2[d_JK2$id == 1, ]
  out <- compute_Vxi(x, phi, K)
  expectedOut <- rep(0, K * m)
  for (i in 1:K)
  {
    idx <- which(x$state == i)
    idx1 <- idx[x$time[idx] <= 5]
    expectedOut[1 + (i - 1) * m] <- sum(pmin(x$time[idx1 + 1], 5) - x$time[idx1], na.rm = TRUE)

    idx2 <- idx[x$time[idx + 1] > 5]
    expectedOut[2 + (i - 1) * m] <- sum(x$time[idx2 + 1] - pmax(x$time[idx2], 5), na.rm = TRUE)
  }

  expect_length(out, K * m)
  expect_lte(max(abs(out - expectedOut)), 1e-5)
})

test_that("computeVmatrix keeps the id order", {
  all_dat_P1S1 <- data.frame(state = c(1, 2, 1), time = c(0, 0.5, 1), id = rep("P1S1", 3))
  all_dat_P1S2 <- data.frame(state = c(1, 2, 1), time = c(0, 0.4, 1), id = rep("P1S2", 3))
  all_dat_P1S3 <- data.frame(state = c(1, 2, 1), time = c(0, 0.6, 1), id = rep("P1S3", 3))
  all_dat_P2S1 <- data.frame(state = c(1, 3, 1), time = c(0, 0.5, 1), id = rep("P2S1", 3))
  all_dat_P2S2 <- data.frame(state = c(1, 3, 1), time = c(0, 0.4, 1), id = rep("P2S2", 3))
  all_dat_P2S3 <- data.frame(state = c(1, 3, 1), time = c(0, 0.6, 1), id = rep("P2S3", 3))
  all_dat_P3S1 <- data.frame(state = c(3, 1, 1), time = c(0, 0.5, 1), id = rep("P3S1", 3))
  all_dat_P3S2 <- data.frame(state = c(3, 1, 1), time = c(0, 0.4, 1), id = rep("P3S2", 3))
  all_dat_P3S3 <- data.frame(state = c(3, 1, 1), time = c(0, 0.6, 1), id = rep("P3S3", 3))

  # Original dataset
  all_dat_a <- rbind(all_dat_P1S1, all_dat_P1S2, all_dat_P1S3,
                     all_dat_P2S1, all_dat_P2S2, all_dat_P2S3,
                     all_dat_P3S1, all_dat_P3S2, all_dat_P3S3)

  all_dat_b <- rbind(all_dat_P2S1, all_dat_P2S2, all_dat_P2S3,
                     all_dat_P3S1, all_dat_P3S2, all_dat_P3S3,
                     all_dat_P1S1, all_dat_P1S2, all_dat_P1S3)
  uniqueIdA <- unique(all_dat_a$id)
  uniqueIdB <- unique(all_dat_b$id)

  basisobj <- create.bspline.basis(c(0, 1), nbasis = 6, norder = 1)

  Va <- computeVmatrix(all_dat_a, basisobj, K = 3, uniqueId = uniqueIdA, nCores = 1, verbose = FALSE)
  Vb <- computeVmatrix(all_dat_b, basisobj, K = 3, uniqueId = uniqueIdB, nCores = 1, verbose = FALSE)

  expect_equal(Va, Vb[c(7, 8, 9, 1, 2, 3, 4, 5, 6), ])
})

test_that("computeUmatrix keeps the id order", {
  all_dat_P1S1 <- data.frame(state = c(1, 2, 1), time = c(0, 0.5, 1), id = rep("P1S1", 3))
  all_dat_P1S2 <- data.frame(state = c(1, 2, 1), time = c(0, 0.4, 1), id = rep("P1S2", 3))
  all_dat_P1S3 <- data.frame(state = c(1, 2, 1), time = c(0, 0.6, 1), id = rep("P1S3", 3))
  all_dat_P2S1 <- data.frame(state = c(1, 3, 1), time = c(0, 0.5, 1), id = rep("P2S1", 3))
  all_dat_P2S2 <- data.frame(state = c(1, 3, 1), time = c(0, 0.4, 1), id = rep("P2S2", 3))
  all_dat_P2S3 <- data.frame(state = c(1, 3, 1), time = c(0, 0.6, 1), id = rep("P2S3", 3))
  all_dat_P3S1 <- data.frame(state = c(3, 1, 1), time = c(0, 0.5, 1), id = rep("P3S1", 3))
  all_dat_P3S2 <- data.frame(state = c(3, 1, 1), time = c(0, 0.4, 1), id = rep("P3S2", 3))
  all_dat_P3S3 <- data.frame(state = c(3, 1, 1), time = c(0, 0.6, 1), id = rep("P3S3", 3))

  # Original dataset
  all_dat_a <- rbind(all_dat_P1S1, all_dat_P1S2, all_dat_P1S3,
                     all_dat_P2S1, all_dat_P2S2, all_dat_P2S3,
                     all_dat_P3S1, all_dat_P3S2, all_dat_P3S3)

  all_dat_b <- rbind(all_dat_P2S1, all_dat_P2S2, all_dat_P2S3,
                     all_dat_P3S1, all_dat_P3S2, all_dat_P3S3,
                     all_dat_P1S1, all_dat_P1S2, all_dat_P1S3)
  uniqueIdA <- unique(all_dat_a$id)
  uniqueIdB <- unique(all_dat_b$id)

  basisobj <- create.bspline.basis(c(0, 1), nbasis = 6, norder = 1)

  Ua <- computeUmatrix(all_dat_a, basisobj, K = 3, uniqueId = uniqueIdA, nCores = 1, verbose = FALSE)
  Ub <- computeUmatrix(all_dat_b, basisobj, K = 3, uniqueId = uniqueIdB, nCores = 1, verbose = FALSE)

  expect_equal(Ua, Ub[c(7, 8, 9, 1, 2, 3, 4, 5, 6), ])
})

test_that("compute_optimal_encoding throws error", {
  set.seed(42)

  # create basis object
  m <- 10
  b <- create.bspline.basis(c(0, 1), nbasis = m, norder = 4)

  expect_error(compute_optimal_encoding(d_JK2, 3, nCores = 1, verbose = TRUE), regexp = "basisobj is not a basis object.")

  expect_error(compute_optimal_encoding(d_JK2, b, nCores = 0, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, nCores = 2.5, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, nCores = NA, verbose = TRUE), regexp = "nCores must be an integer > 0.")
  expect_error(compute_optimal_encoding(d_JK2, b, nCores = NaN, verbose = TRUE), regexp = "nCores must be an integer > 0.")

  expect_error(compute_optimal_encoding(d_JK2, b, nCores = 1, verbose = 2), regexp = "verbose must be either TRUE or FALSE.")
})


test_that("compute_optimal_encoding works", {
  set.seed(42)
  n <- 200
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- cut_data(d, Tmax)
  row.names(dT) <- NULL

  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  expect_silent(fmca <- compute_optimal_encoding(dT, b, computeCI = FALSE, nCores = 1, verbose = FALSE))

  expect_type(fmca, "list")
  expect_named(fmca, c("eigenvalues", "alpha", "pc", "F", "G", "V", "basisobj", "label", "pt", "runTime"))

  # eigenvalues
  expect_length(fmca$eigenvalues, K * m)
  trueEigVal <- 1 / ((1:m) * (2:(m + 1)))
  expect_lte(max(abs(fmca$eigenvalues[1:m] - trueEigVal)), 0.01)

  # alpha
  expect_type(fmca$alpha, "list")
  expect_length(fmca$alpha, m * K)
  expect_equal(dim(fmca$alpha[[1]]), c(m, K))

  # pc
  expect_equal(dim(fmca$pc), c(n, m * K))

  # F
  expect_equal(dim(fmca$F), c(m * K, m * K))

  # G
  expect_equal(dim(fmca$G), c(m * K, m * K))

  # V
  expect_equal(dim(fmca$V), c(n, m * K))

  # basisobj
  expect_equal(fmca$basisobj, b)

  # label
  expect_equal(fmca$label, data.frame(label = 0:1, code = 1:2))
})

skip_on_cran()

## data and results used in the next tests
set.seed(42)
n <- 10
Tmax <- 1
K <- 2
m <- 5
d <- generate_2State(n)
dT <- cut_data(d, Tmax)
row.names(dT) <- NULL

b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
fmca <- compute_optimal_encoding(dT, b, nCores = 1, computeCI = FALSE, verbose = FALSE)
##

test_that("summary.cfda does not produce warnings/errors", {
  expect_warning(summary(fmca), regexp = NA)
  expect_error(summary(fmca), regexp = NA)
})

test_that("print.cfda does not produce warnings/errors", {
  expect_warning(print(fmca), regexp = NA)
  expect_error(print(fmca), regexp = NA)
})

test_that("compute_optimal_encoding works verbose", {
  skip_on_cran()
  expect_output(encoding <- compute_optimal_encoding(dT[dT$id <= 10, ], b, computeCI = FALSE, nCores = 1, verbose = TRUE))
})

test_that("compute_optimal_encoding throws a warning when the basis is not well suited", {
  skip_on_cran()
  data_msm <- data.frame(id = rep(1:2, each = 3), time = c(0, 3, 5, 0, 4, 5), state = c(1, 2, 2, 1, 2, 2))
  b <- create.bspline.basis(c(0, 5), nbasis = 3, norder = 2)

  expect_warning(
    {
      fmca <- compute_optimal_encoding(data_msm, b, computeCI = FALSE, nCores = 1)
    },
    regexp = "The F matrix contains at least one column of 0s. At least one state is not present in the support of one basis function. Corresponding coefficients in the alpha output will have a 0 value.",
    fixed = TRUE
  )
})



test_that("get_encoding throws error", {
  fmca <- list(alpha = rep(1, 5))
  class(fmca) <- "fmca"

  expect_error(get_encoding(3), regexp = "x must be a fmca object.")

  expect_error(get_encoding(fmca, harm = c(1, 2)), regexp = "harm must be an integer between 1 and the number of components.")
  expect_error(get_encoding(fmca, harm = 2.5), regexp = "harm must be an integer between 1 and the number of components.")
  expect_error(get_encoding(fmca, harm = 0), regexp = "harm must be an integer between 1 and the number of components.")
  expect_error(get_encoding(fmca, harm = NA), regexp = "harm must be an integer between 1 and the number of components.")
  expect_error(get_encoding(fmca, harm = NaN), regexp = "harm must be an integer between 1 and the number of components.")
  expect_error(get_encoding(fmca, harm = 10), regexp = "harm must be an integer between 1 and the number of components.")


  expect_error(get_encoding(fmca, harm = 1, nx = 0), regexp = "nx must be a positive integer.")
  expect_error(get_encoding(fmca, harm = 1, nx = c(3, 2)), regexp = "nx must be a positive integer.")
  expect_error(get_encoding(fmca, harm = 1, nx = NA), regexp = "nx must be a positive integer.")
  expect_error(get_encoding(fmca, harm = 1, nx = NaN), regexp = "nx must be a positive integer.")
})


test_that("get_encoding works", {
  out <- get_encoding(fmca, fdObject = TRUE)
  expect_s3_class(out, "fd")

  out <- get_encoding(fmca, harm = 3, fdObject = TRUE)
  expect_s3_class(out, "fd")

  out <- get_encoding(fmca, fdObject = FALSE)
  expect_named(out, c("x", "y"))
  expect_equal(dim(out$y), c(128, 2))
  expect_length(out$x, 128)

  out <- get_encoding(fmca, harm = 3, fdObject = FALSE)
  expect_named(out, c("x", "y"))
  expect_equal(dim(out$y), c(128, 2))
  expect_length(out$x, 128)

  out <- get_encoding(fmca, fdObject = FALSE, nx = 100)
  expect_named(out, c("x", "y"))
  expect_equal(dim(out$y), c(100, 2))
  expect_length(out$x, 100)
})

test_that("plot.fmca does not produce warnings", {
  expect_warning(plot(fmca), regexp = NA)
  expect_warning(plot(fmca, harm = 3, col = c("red", "blue")), regexp = NA)
  expect_warning(plot(fmca, addCI = TRUE), regexp = NA)
})


test_that("plotComponent throws error", {
  fmca <- list(pc = matrix(nrow = 3, ncol = 5))
  class(fmca) <- "fmca"

  expect_error(plotComponent(3), regexp = "x must be a fmca object.")

  expect_error(plotComponent(fmca, comp = 1), regexp = "comp must be a vector of positive integers of length 2.")
  expect_error(plotComponent(fmca, comp = c(1, NA)), regexp = "comp must be a vector of positive integers of length 2.")
  expect_error(plotComponent(fmca, comp = c(1, NaN)), regexp = "comp must be a vector of positive integers of length 2.")
  expect_error(plotComponent(fmca, comp = c(1, 1.5)), regexp = "comp must be a vector of positive integers of length 2.")
  expect_error(plotComponent(fmca, comp = c(1, 6)), regexp = "comp must be a vector of positive integers of length 2.")
})


test_that("plotComponent does not produce warnings", {
  expect_warning(plotComponent(fmca, addNames = TRUE), regexp = NA)
  expect_warning(plotComponent(fmca, comp = c(2, 3), addNames = FALSE), regexp = NA)
  expect_warning(plotComponent(fmca, shape = 23), regexp = NA)
})

test_that("plotEigenvalues throws error", {
  expect_error(plotEigenvalues(2), regexp = "x must be a fmca object.")
})

test_that("plotEigenvalues does not produce warnings", {
  expect_warning(plotEigenvalues(fmca, cumulative = FALSE, normalize = FALSE), regexp = NA)
  expect_warning(plotEigenvalues(fmca, cumulative = FALSE, normalize = TRUE), regexp = NA)
  expect_warning(plotEigenvalues(fmca, cumulative = TRUE, normalize = FALSE), regexp = NA)
  expect_warning(plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE, shape = 23), regexp = NA)
})
