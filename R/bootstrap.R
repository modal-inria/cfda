
# compute bootstrap estimates of encoding
#
# @param Uval output of computeUmatrix function
# @param V output of computeVmatrix function
# @param K number of states
# @param nBasis number of basis
# @param label label element from stateToInteger function
# @param nId number of individuals
# @param propBootstrap proportion of individuals used for bootstrap sample
# @param nBootstrap number of bootstrap samples
# @param output of getSignReference function the full encoding
# @param verbose if TRUE print some information
#
# @return a list of computeEncoding function output
#
# @author Quentin Grimonprez
computeBootStrapEncoding <- function(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, signReference, verbose) {
  outEnc <- list()
  t3 <- proc.time()
  if (verbose) {
    cat("---- Compute Bootstrap Encoding:\n")
  }

  for (i in seq_len(nBootstrap)) {
    if (verbose) {
      cat("*")
    }
    idToKeep <- sample(nId, floor(propBootstrap * nId), replace = TRUE)

    try({
      outEnc[[i]] <- computeEncoding(Uval[idToKeep, ], V[idToKeep, ], K, nBasis, idToKeep,
                                     label, verbose = FALSE, manage0 = TRUE)
    })

    # outEnc[[i]] = c(outEnc[[i]] , list(basisobj = basisobj))
    # class(outEnc[[i]]) = "fmca"
  }

  # reorder alpha, pc such that representation have the same sign for each bootstrap sample
  outEnc <- unifySign(outEnc, signReference)


  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }

  return(outEnc)
}

# get the sign of the alpha of the full encoding to later ensure that bootstrap samples have the same sign
# @author Quentin Grimonprez
getSignReference <- function(alpha) {
  pos <- rep(0, length(alpha))
  isNeg <- rep(FALSE, length(alpha))
  for (i in seq_along(alpha)) {
    pos[i] <- which.max(abs(alpha[[i]]))
    isNeg[i] <- Re(alpha[[i]][pos[i]]) < 0
  }

  return(list(position = pos, isNegative = isNeg, allNegative = lapply(alpha, function(x) Re(x) < 0)))
}


# eigenvectors are equivalent up to the sign
# we try to have the same sign for each bootstrap sample
#
# @param out a list of computeEncoding function output
#
# @author Quentin Grimonprez
unifySign <- function(out, signReference) {
  for (i in seq_along(out)) {
    if (!is.null(out[[i]])) { # an output can be NULL due to inversion problem
      # we look if the element at the given position is negative
      for (j in seq_along(out[[i]]$alpha)) {
        signNeg <- Re(out[[i]]$alpha[[j]][signReference$position[j]]) < 0

        # if there is a NA at he given position, we try an other position
        if (!is.na(signNeg)) {
          if (signNeg != signReference$isNegative[j]) {
            out[[i]]$alpha[[j]] <- out[[i]]$alpha[[j]] * -1
            out[[i]]$pc[, j] <- out[[i]]$pc[, j] * -1
          }
        } else {
          newPos <- which.max(out[[i]]$alpha[[j]])
          signNeg <- Re(out[[i]]$alpha[[j]][newPos]) < 0

          if (signNeg != signReference$allNegative[[j]][newPos]) {
            out[[i]]$alpha[[j]] <- out[[i]]$alpha[[j]] * -1
            out[[i]]$pc[, j] <- out[[i]]$pc[, j] * -1
          }
        }
      }
    }
  }


  return(out)
}


# Compute the variance of alpha
#
# @param bootEncoding output of computeBootStrapEncoding function
# @param nState number of states
# @param nBasis number of basis functions
#
# @return a list (length number of harmonics) of list (length number of states) of variance matrix
#
# @author Quentin Grimonprez
computeVarianceAlpha <- function(bootEncoding, nState, nBasis) {
  nHarm <- nState * nBasis

  varAlpha <- list()
  for (harm in seq_len(nHarm)) {
    varAlpha[[harm]] <- list()
    for (iState in seq_len(nState)) {
      tryCatch(
        {
          varAlpha[[harm]][[iState]] <- var(do.call(
            rbind,
            lapply(bootEncoding, function(x) {
              if (!is.null(x$alpha[[harm]][, iState])) {
                return(Re(x$alpha[[harm]][, iState]))
              }
            })
          ),
          use = "pairwise.complete.obs"
          )
        },
        error = function(e) e
      )
    }
  }

  return(varAlpha)
}


# compute the variance of encoding
#
# a_x(t) = sum_i alpha_ix * phi_i(t)
# Var(a_x(t)) = sum_i var(alpha_ix) * phi_i^2(t) + sum_{i<j}  2 * phi_i(t) * phi_j(t) * cov(alpha_ix, alpha_jx)
#
# @author Quentin Grimonprez
computeVarianceEncoding <- function(varAlpha, basisobj, harm = 1, nx = 200) {
  nBasis <- basisobj$nbasis
  phi <- fd(diag(nBasis), basisobj)
  nState <- length(varAlpha[[harm]])

  timeVal <- seq(basisobj$rangeval[1], basisobj$rangeval[2], length = nx)

  Phi <- matrix(nrow = length(timeVal), ncol = nBasis)
  for (i in seq_len(nBasis)) {
    Phi[, i] <- eval.fd(timeVal, phi[i])
  }

  funcVar <- list()
  for (iState in seq_len(nState)) {
    funcVar[[iState]] <- rep(NA, nx)
    for (j in seq_along(timeVal)) {
      varAlpha[[harm]][[iState]][is.na(varAlpha[[harm]][[iState]])] <- 0
      funcVar[[iState]][j] <- Phi[j, , drop = FALSE] %*% varAlpha[[harm]][[iState]] %*% t(Phi[j, , drop = FALSE])
    }
  }

  return(funcVar)
}
