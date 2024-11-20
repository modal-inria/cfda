#' Compute the optimal encoding for each state
#'
#' Compute the optimal encoding for categorical functional data using an extension of the multiple correspondence analysis
#' to a stochastic process.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state. All individuals must begin at the same time T0 and end at the same time Tmax
#' (use \code{\link{cut_data}}).
#' @param basisobj basis created using the \code{fda} package (cf. \code{\link[fda]{create.basis}}).
#' @param computeCI if TRUE, perform a bootstrap to estimate the variance of encoding functions coefficients
#' @param nBootstrap number of bootstrap samples
#' @param propBootstrap size of bootstrap samples relative to the number of individuals: propBootstrap * number of individuals
#' @param method computation method: "parallel" or "precompute": precompute all integrals
#' (efficient when the number of unique time values is low)
#' @param verbose if TRUE print some information
#' @param nCores number of cores used for parallelization (only if method == "parallel"). Default is half the cores.
#' @param ... parameters for \code{\link{integrate}} function (see details).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{eigenvalues} eigenvalues
#'   \item \code{alpha} optimal encoding coefficients associated with each eigenvectors
#'   \item \code{pc} principal components
#'   \item \code{F} matrix containing the \eqn{F_{(x,i)(y,j)}}
#'   \item \code{V} matrix containing the \eqn{V_{(x,i)}}
#'   \item \code{G} covariance matrix of \code{V}
#'   \item \code{basisobj} \code{basisobj} input parameter
#'   \item \code{pt} output of \code{\link{estimate_pt}} function
#'   \item \code{bootstrap} Only if \code{computeCI = TRUE}. Output of every bootstrap run
#'   \item \code{varAlpha} Only if \code{computeCI = TRUE}. Variance of alpha parameters
#'   \item \code{runTime} Total elapsed time
#' }
#'
#' @details
#' See the vignette for the mathematical background: \code{RShowDoc("cfda", package = "cfda")}
#'
#' Extra parameters (\emph{...}) for the \code{\link{integrate}} function can be:
#' \itemize{
#'  \item \emph{subdivisions} the maximum number of subintervals.
#'  \item \emph{rel.tol} relative accuracy requested.
#'  \item \emph{abs.tol} absolute accuracy requested.
#' }
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 5
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(
#'   n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax,
#'   labels = c("A", "C", "G", "T")
#' )
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 5
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#'
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#' summary(encoding)
#'
#' # plot the optimal encoding
#' plot(encoding)
#'
#' # plot the two first components
#' plotComponent(encoding, comp = c(1, 2))
#'
#' # extract the optimal encoding
#' get_encoding(encoding, harm = 1)
#' @seealso \code{link{plot.fmca}} \code{link{print.fmca}} \code{link{summary.fmca}} \code{link{plotComponent}} \code{link{get_encoding}}
#'
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @references
#' \itemize{
#'   \item Deville J.C. (1982) Analyse de données chronologiques qualitatives : comment analyser des calendriers ?,
#' Annales de l'INSEE, No 45, p. 45-104.
#'   \item Deville J.C. et  Saporta G. (1980) Analyse harmonique qualitative, DIDAY et al. (editors), Data Analysis and
#' Informatics, North Holland, p. 375-389.
#'   \item Saporta G. (1981) Méthodes exploratoires d'analyse de données temporelles, Cahiers du B.U.R.O, Université
#' Pierre et Marie Curie, 37-38, Paris.
#'   \item Preda C, Grimonprez Q, Vandewalle V. Categorical Functional Data Analysis. The cfda R Package.
#' Mathematics. 2021; 9(23):3074. https://doi.org/10.3390/math9233074
#' }
#'
#' @family encoding functions
#'
#' @export
compute_optimal_encoding <- function(
    data, basisobj, computeCI = TRUE, nBootstrap = 50, propBootstrap = 1, method = c("precompute", "parallel"),
    verbose = TRUE, nCores = max(1, ceiling(detectCores() / 2)), ...) {
  t1 <- proc.time()

  ## check parameters
  check_compute_optimal_encoding_parameters(data, basisobj, nCores, verbose, computeCI, nBootstrap, propBootstrap)
  method <- match.arg(method)
  nCores <- min(max(1, nCores), detectCores() - 1)
  ## end check

  # used to determine the moments where the probability is 0 in plot.fmca
  pt <- estimate_pt(data)

  if (verbose) {
    cat("######### Compute encoding #########\n")
  }


  # change state as integer
  out <- stateToInteger(data$state)
  data$state <- out$state
  label <- out$label
  rm(out)

  uniqueId <- as.character(unique(data$id))
  nId <- length(uniqueId)
  K <- length(label$label)
  nBasis <- basisobj$nbasis

  if (verbose) {
    cat(paste0("Number of individuals: ", nId, "\n"))
    cat(paste0("Number of states: ", K, "\n"))
    cat(paste0("Basis type: ", basisobj$type, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
    cat(paste0("Number of cores: ", nCores, "\n"))
    cat(paste0("Method: ", method, "\n"))
  }

  if (method == "precompute") {
    uniqueTime <- sort(unique(data$time))

    V <- computeVmatrix2(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose, ...)

    Uval <- computeUmatrix2(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose, ...)
  } else {
    V <- computeVmatrix(data, basisobj, K, uniqueId, nCores, verbose, ...)

    Uval <- computeUmatrix(data, basisobj, K, uniqueId, nCores, verbose, ...)
  }

  fullEncoding <- computeEncoding(Uval, V, K, nBasis, uniqueId, label, verbose, manage0 = TRUE)


  if (computeCI) {
    signRef <- getSignReference(fullEncoding$alpha)

    bootEncoding <- computeBootStrapEncoding(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, signRef, verbose)
    if (length(bootEncoding) == 0) {
      warning("All bootstrap samples return an error. Try to change the basis.")
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt))
    } else {
      varAlpha <- computeVarianceAlpha(bootEncoding, K, nBasis)
      out <- c(fullEncoding, list(
        V = V, basisobj = basisobj, label = label, pt = pt,
        bootstrap = bootEncoding, varAlpha = varAlpha
      ))
    }
  } else {
    out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt))
  }

  if (is.complex(out$eigenvalues)) {
    warning("Eigenvalues contain complex values. Only the real part is returned.")
    out$eigenvalues <- Re(out$eigenvalues)
    out$pc <- Re(out$pc)
  }

  class(out) <- "fmca"
  t2 <- proc.time()

  out$runTime <- as.numeric((t2 - t1)[3])
  if (verbose) {
    cat(paste0("Run Time: ", round(out$runTime, 2), "s\n"))
  }

  return(out)
}

check_compute_optimal_encoding_parameters <- function(data, basisobj, nCores, verbose, computeCI, nBootstrap, propBootstrap) {
  checkData(data)
  checkDataBeginTime(data)
  checkDataEndTmax(data)
  checkDataNoDuplicatedTimes(data)

  if (!is.basis(basisobj)) {
    stop("basisobj is not a basis object.")
  }
  if (any(is.na(nCores)) || !is.whole.number(nCores) || (nCores < 1)) {
    stop("nCores must be an integer > 0.")
  }
  checkLogical(verbose, "verbose")
  checkLogical(computeCI, "computeCI")
  if (computeCI) {
    checkInteger(nBootstrap, minValue = 0, minEqual = FALSE, paramName = "nBootstrap")

    if (any(is.na(propBootstrap)) || !is.numeric(propBootstrap) || (length(propBootstrap) != 1) || (propBootstrap > 1) || (
      propBootstrap <= 0)) {
      stop("propBootstrap must be a real between 0 and 1.")
    }
  }
}

# return a matrix with nId rows and nBasis * nState columns
computeVmatrix <- function(data, basisobj, K, uniqueId, nCores, verbose, ...) {
  nBasis <- basisobj$nbasis

  phi <- fd(diag(nBasis), basisobj) # basis function as functional data

  # declare parallelization
  if (nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }


  if (verbose) {
    cat("---- Compute V matrix:\n")
    pbo <- pboptions(type = "timer", char = "=")
  } else {
    pboptions(type = "none")
  }
  t2 <- proc.time()


  # we compute V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  V <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Vxi, phi = phi, K = K, ...)[uniqueId])
  rownames(V) <- NULL

  t3 <- proc.time()

  if (verbose) {
    cat(paste0("\nDONE in ", round((t3 - t2)[3], 2), "s\n"))
  }

  # stop parallelization
  if (nCores > 1) {
    stopCluster(cl)
  }

  return(V)
}



# compute_Vxi  (Vxi = \int_0^Tphi_i(t)X_t=x dt)
#
# @param x one individual (id, time, state)
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
# @param ... parameters for integrate function
#
# @return vector of size K*nBasis: V[(x=1,i=1)],
# V[(x=1,i=2)],..., V[(x=1,i=nBasis)], V[(x=2,i=1)], V[(x=2,i=2)], ...
#
# @examples
# K <- 4
# Tmax <- 6
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
# d_JK2 <- cut_data(d_JK, Tmax)
#
# m <- 6
# b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# I <- diag(rep(1, m))
# phi <- fd(I, b)
# compute_Vxi(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda
compute_Vxi <- function(x, phi, K, ...) {
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis) # V11, V12,...V1m, V21, V22, ..., V2m, ... etc VK1... VKm

  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]

    for (j in seq_len(nBasis)) { # j = la base
      ind <- (state - 1) * nBasis + j
      aux[ind] <- aux[ind] +
        integrate(
          function(t) {
            eval.fd(t, phi[j])
          },
          lower = x$time[u], upper = x$time[u + 1],
          stop.on.error = FALSE, ...
        )$value
    }
  }

  return(aux)
}


# return a matrix with nId rows and nBasis * nState columns
computeVmatrix2 <- function(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose, ...) {
  nBasis <- basisobj$nbasis

  phi <- fd(diag(nBasis), basisobj) # basis function as functional data

  t2 <- proc.time()

  index <- data.frame(seq_along(uniqueTime), row.names = uniqueTime)

  # V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  integrals <- compute_integral_V(phi, uniqueTime, verbose)
  # V <- do.call(
  #   rbind,
  #   pblapply(cl = cl, split(data, data$id), fill_V, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  # )
  V <- do.call(
    rbind,
    lapply(split(data, data$id), fill_V, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  )
  rownames(V) <- NULL

  # V <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Vxi, phi = phi, K = K, ...)[uniqueId])
  # rownames(V) <- NULL

  t3 <- proc.time()

  if (verbose) {
    cat(paste0("\nDONE in ", round((t3 - t2)[3], 2), "s\n"))
  }

  return(V)
}


compute_integral_V <- function(phi, uniqueTime, verbose, ...) {
  nBasis <- phi$basis$nbasis
  if (verbose) {
    pb <- timerProgressBar(min = 0, max = nBasis, width = 50)
    on.exit(close(pb))
    jj <- 1
  }
  integrals <- list()
  for (i in seq_len(nBasis)) { # TODO parallel
    integrals[[i]] <- rep(0., length(uniqueTime))
    for (ii in seq_len(length(uniqueTime) - 1)) {
      integrals[[i]][ii] <- integrate(
        function(t) {
          eval.fd(t, phi[i])
        },
        lower = uniqueTime[ii], upper = uniqueTime[ii + 1],
        stop.on.error = FALSE, ...
      )$value
    }
    if (verbose) {
      setTimerProgressBar(pb, jj)
      jj <- jj + 1
    }
  }

  return(integrals)
}


fill_V <- function(x, integrals, index, K, nBasis) {
  aux <- rep(0, K * nBasis)
  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    s <- as.character(x$time[u])
    e <- as.character(x$time[u + 1])
    for (i in seq_len(nBasis)) {
      integral <- sum(integrals[[i]][index[s, ]:(index[e, ] - 1)])
      ind <- (state - 1) * nBasis + i
      aux[ind] <- aux[ind] + integral
    }
  }
  return(aux)
}


# return a matrix with nId rows and nBasis * nState columns
computeUmatrix <- function(data, basisobj, K, uniqueId, nCores, verbose, ...) {
  nBasis <- basisobj$nbasis

  phi <- fd(diag(nBasis), basisobj) # basis functions as functional data

  # declare parallelization
  if (nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }

  t3 <- proc.time()
  if (verbose) {
    cat(paste0("---- Compute U matrix:\n"))
    pbo <- pboptions(type = "timer", char = "=")
  } else {
    pbo <- pboptions(type = "none")
  }

  Uval <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Uxij, phi = phi, K = K, ...)[uniqueId])

  # stop parallelization
  if (nCores > 1) {
    stopCluster(cl)
  }


  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }


  return(Uval)
}



# compute Uxij
#
# compute int_0^T phi_i(t) phi_j(t) 1_{X_{t}=x}
#
# @param x one individual (id, time, state)
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
# @param ... parameters for integrate function
#
# @return vector of size K*nBasis*nBasis: U[(x=1,i=1),(x=1,j=1)],
# U[(x=1,i=1), (x=1,j=2)],..., U[(x=1,i=1), (x=1,j=nBasis)], U[(x=2,i=1), (x=2,j=1)], U[(x=2,i=1), (x=2,j=2)], ...
#
# @examples
# K <- 4
# Tmax <- 6
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
# d_JK2 <- cut_data(d_JK, Tmax)
#
# m <- 6
# b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# I <- diag(rep(1, m))
# phi <- fd(I, b)
# compute_Uxij(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda, Quentin Grimonprez
compute_Uxij <- function(x, phi, K, ...) {
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)

  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    for (i in seq_len(nBasis)) {
      for (j in i:nBasis) { # symmetry between i and j
        integral <- integrate(
          function(t) {
            eval.fd(t, phi[i]) * eval.fd(t, phi[j])
          },
          lower = x$time[u], upper = x$time[u + 1],
          stop.on.error = FALSE, ...
        )$value

        ind <- (state - 1) * nBasis * nBasis + (i - 1) * nBasis + j
        aux[ind] <- aux[ind] + integral

        # when i == j, we are on the diagonal of the matrix, no symmetry to apply
        if (i != j) {
          ind <- (state - 1) * nBasis * nBasis + (j - 1) * nBasis + i
          aux[ind] <- aux[ind] + integral
        }
      }
    }
  }

  return(aux)
}


# return a matrix with nId rows and nBasis * nState columns
computeUmatrix2 <- function(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose, ...) {
  nBasis <- basisobj$nbasis

  phi <- fd(diag(nBasis), basisobj) # basis function as functional data
  # TODO https://r-coder.com/progress-bar-r/

  t3 <- proc.time()
  if (verbose) {
    cat(paste0("---- Compute U matrix:\n"))
  }
  index <- data.frame(seq_along(uniqueTime), row.names = uniqueTime)

  integrals <- compute_integral_U(phi, uniqueTime, verbose)

  # Uval <- do.call(
  #   rbind,
  #   pblapply(cl = cl, split(data, data$id), fill_U, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  # )

  cl <- NULL

  if (verbose) {
    pbo <- pboptions(type = "timer", char = "=")
  } else {
    pbo <- pboptions(type = "none")
  }


  Uval <- do.call(
    rbind,
    pblapply(cl = cl, split(data, data$id), fill_U, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  )

  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }

  return(Uval)
}


compute_integral_U <- function(phi, uniqueTime, verbose, ...) {
  nBasis <- phi$basis$nbasis
  if (verbose) {
    nIter <- nBasis * (1 + nBasis) / 2
    pb <- timerProgressBar(min = 0, max = nIter, width = 50)
    on.exit(close(pb))
    jj <- 1
  }

  integrals <- list()
  for (i in seq_len(nBasis)) { # TODO parallel
    integrals[[i]] <- list()
    for (j in i:nBasis) {
      integrals[[i]][[j]] <- rep(0, length(uniqueTime))
      for (ii in seq_len(length(uniqueTime) - 1)) {
        integrals[[i]][[j]][ii] <- integrate(
          function(t) {
            eval.fd(t, phi[i]) * eval.fd(t, phi[j])
          },
          lower = uniqueTime[ii], upper = uniqueTime[ii + 1],
          stop.on.error = FALSE, ...
        )$value
      }
      if (verbose) {
        setTimerProgressBar(pb, jj)
        jj <- jj + 1
      }
    }
  }
  return(integrals)
}

fill_U <- function(x, integrals, index, K, nBasis) {
  aux <- rep(0, K * nBasis * nBasis)
  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    s <- as.character(x$time[u])
    e <- as.character(x$time[u + 1])
    for (i in seq_len(nBasis)) {
      for (j in i:nBasis) { # symmetry between i and j
        integral <- sum(integrals[[i]][[j]][index[s, ]:(index[e, ] - 1)])
        ind <- (state - 1) * nBasis * nBasis + (i - 1) * nBasis + j
        aux[ind] <- aux[ind] + integral

        # when i == j, we are on the diagonal of the matrix, no symmetry to apply
        if (i != j) {
          ind <- (state - 1) * nBasis * nBasis + (j - 1) * nBasis + i
          aux[ind] <- aux[ind] + integral
        }
      }
    }
  }
  return(aux)
}

# @author Cristian Preda, Quentin Grimonprez
computeEncoding <- function(Uval, V, K, nBasis, uniqueId, label, verbose, manage0 = FALSE) {
  t4 <- proc.time()
  if (verbose) {
    cat(paste0("---- Compute encoding: "))
  }

  G <- cov(V) * (nrow(V) - 1) / nrow(V)

  # create F matrix
  Fval <- colMeans(Uval)
  Fmat <- matrix(0, ncol = K * nBasis, nrow = K * nBasis) # diagonal-block matrix with K blocks of size nBasis*nBasis
  for (i in seq_len(K)) {
    Fmat[((i - 1) * nBasis + 1):(i * nBasis), ((i - 1) * nBasis + 1):(i * nBasis)] <-
      matrix(Fval[((i - 1) * nBasis * nBasis + 1):(i * nBasis * nBasis)], ncol = nBasis, byrow = TRUE)
  }

  # save matrices before modifying them
  outMat <- list(F = Fmat, G = G)


  # manage the case where there is a column full of 0 (non invertible matrix)
  # if TRUE, we remove the 0-columns (and rows) and throw a warning
  # otherwise we throw an error and the process stops
  if (manage0) {
    # column full of 0
    ind0 <- (colSums(Fmat == 0) == nrow(Fmat))
    if (sum(ind0) > 0) {
      warning(paste(
        "The F matrix contains at least one column of 0s.",
        "At least one state is not present in the support of one basis function.",
        "Corresponding coefficients in the alpha output will have a 0 value."
      ))
    }

    F05 <- t(mroot(Fmat[!ind0, !ind0])) # F  = t(F05)%*%F05
    G <- G[!ind0, !ind0]
    V <- V[, !ind0]
  } else {
    ind0 <- (colSums(Fmat == 0) == nrow(Fmat))

    # res = eigen(solve(F)%*%G)
    F05 <- t(mroot(Fmat)) # F  = t(F05)%*%F05

    if (any(dim(F05) != rep(K * nBasis, 2))) {
      cat("\n")
      if (any(colSums(Fmat) == 0)) {
        stop(paste(
          "F matrix is not invertible. In the support of each basis function,",
          "each state must be present at least once (p(x_t) != 0 for t in the support).",
          "You can try to change the basis."
        ))
      }

      stop("F matrix is not invertible. You can try to change the basis.")
    }
  }


  invF05 <- solve(F05)
  # res <- eigen(F05 %*% solve(F) %*% G %*% solve(F05))
  res <- eigen(t(invF05) %*% G %*% invF05)

  # eigenvectors (they give the coefficients of the m=nBasis encoding, for each eigenvalue) as a list of matrices of size m x K
  # The first matrix is the first eigenvector. This vector (1rst column in res$vectors) contains the coefficients of the
  # encoding for the state 1 in the first m (=nBasis) positions, then for the state 2 in the m next positions and so on
  # until the k-th state.
  # I put this first column as and matrix and the coefficients are column of length m

  # transform the matrix of eigenvectors as a list

  # aux1 = split(res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))
  invF05vec <- invF05 %*% res$vectors
  aux1 <- split(invF05vec, rep(seq_len(ncol(res$vectors)), each = nrow(res$vectors)))

  # we create the matrices m x K for each eigenvalue: 1rst column = coefficients for state 1,
  # 2nd col = coefficients for state 2, etc

  alpha <- lapply(aux1, function(w) {
    wb <- rep(NA, nBasis * K)
    wb[!ind0] <- w

    return(matrix(wb, ncol = K, dimnames = list(NULL, label$label)))
  })

  pc <- V %*% invF05vec
  rownames(pc) <- uniqueId

  t5 <- proc.time()

  if (verbose) {
    cat(paste0("\nDONE in ", round((t5 - t4)[3], 2), "s\n"))
  }

  return(list(eigenvalues = res$values, alpha = alpha, pc = pc, F = outMat$F, G = outMat$G))
}
