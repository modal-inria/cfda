compute_Vxi_multi <- function(x, phi, K, stateColumns, ...) {
  aux <- list()
  for (j in seq_along(K)) {
    aux[[j]] <- compute_Vxi(x, phi, K[j], stateColumns[j], ...)
  }
  return(do.call(c, aux))
}

computeVlist <- function(data, phi, K, stateColumns, uniqueId, verbose, nCores, ...) {
  # declare parallelization
  if (nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }

  t3 <- proc.time()
  if (verbose) {
    cat("---- Compute V matrix:\n")
    pbo <- pboptions(type = "timer", char = "=")
  } else {
    pboptions(type = "none")
  }

  V_multi <- do.call(
    rbind,
    pblapply(cl = cl, split(data, data$id), compute_Vxi_multi, phi = phi, K = K, stateColumns = stateColumns, ...)[uniqueId]
  )
  rownames(V_multi) <- NULL

  # stop parallelization
  if (nCores > 1) {
    stopCluster(cl)
  }

  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }

  return(V_multi)
}

compute_Uxij_multi <- function(x, phi, K, stateColumns, verbose = FALSE, ...) {
  m <- phi$basis$nbasis

  integrals <- list()
  for (i in seq_along(K)) {
    integrals[[i]] <- list()
    for (j in seq_along(K)) {
      integrals[[i]][[j]] <- rep(0, m * m * K[i] * K[j])
    }
  }

  for (u in seq_len(nrow(x) - 1)) {
    states <- as.numeric(x[u, stateColumns])
    if (verbose) {
      cat(paste0("------ row ", u, ": state \n"))
      print(states)
    }
    for (l1 in seq_len(m)) {
      for (l2 in seq_len(m)) {
        integral <- integrate(
          function(t) {
            eval.fd(t, phi[l1]) * eval.fd(t, phi[l2])
          },
          lower = x$time[u], upper = x$time[u + 1],
          stop.on.error = FALSE, ...
        )$value
        for (i in seq_along(K)) {
          for (j in seq_along(K)) {
            state1 <- states[i]
            state2 <- states[j]
            num_row <- (state1 - 1) * m + l1
            num_col <- (state2 - 1) * m + l2

            ind2 <- (num_col - 1) * (K[i] * m) + num_row
            integrals[[i]][[j]][ind2] <- integrals[[i]][[j]][ind2] + integral
          }
        }
      }
    }
  }
  return(integrals)
}


compute_U_list_matrix <- function(list_of_Uval, K) {
  U_list <- list_of_Uval[[1]]
  for (i_sample in seq_along(list_of_Uval[-1])) {
    for (i in seq_along(K)) {
      for (j in seq_along(K)) {
        U_list[[i]][[j]] <- rbind(U_list[[i]][[j]], list_of_Uval[[1 + i_sample]][[i]][[j]])
      }
    }
  }
  return(U_list)
}

computeUmean <- function(data, phi, K, stateColumns, uniqueId, verbose, nCores, ...) {
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

  list_of_Uval <- pblapply(cl = cl, split(data, data$id), compute_Uxij_multi, phi = phi, K = K, stateColumns = stateColumns, ...)[uniqueId]

  U_list_matrix <- compute_U_list_matrix(list_of_Uval, K)
  U_mean <- U_list_matrix
  for (i in seq_along(K)) {
    for (j in seq_along(K)) {
      U_mean[[i]][[j]] <- colMeans(U_mean[[i]][[j]])
    }
  }

  # stop parallelization
  if (nCores > 1) {
    stopCluster(cl)
  }

  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }

  return(U_mean)
}

#' Compute the multivariate optimal encoding
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' and several columns for the associated state for each dimension.
#' All individuals must begin at the same time T0 and end at the same time Tmax.
#' @param basisobj basis created using the \code{fda} package (cf. \code{\link{create.basis}}).
#' The same basis is used for every dimension.
#' @param epsilon epsilon added to the diagonal of F in order to invert it. If NULL, an epsilon is computed with regards to F.
#' It can be a vector to test several values.
#' @param stateColumns column names for multivariate states. By default, "state1", "state2", ...
#' @param verbose if TRUE print some information
#' @param nCores number of cores used for parallelization. Default is all cores except 1.
#'
#' @return a fmca object
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' Tmax <- 2
#' x1 <- generate_Markov(n = 500, K = 2)
#' x1 <- cut_data(x1, Tmax = Tmax)
#' x2 <- generate_Markov(n = 500, K = 2)
#' x2 <- cut_data(x2, Tmax = Tmax)
#'
#' x <- list(x1, x2)
#'
#' mvcfd <- convert2mvcfd(x)
#' basisobj <- create.bspline.basis(c(0, Tmax), nbasis = 10, norder = 4)
#'
#' multEnc <- compute_optimal_encoding_multivariate(mvcfd, basisobj, verbose = FALSE, nCores = 1)
#'
#' plot(multEnc)
#' }
#' @export
compute_optimal_encoding_multivariate <- function(
    data, basisobj, epsilon = NULL, stateColumns = NULL, verbose = TRUE, nCores = max(1, detectCores() - 1)) {
  t1 <- proc.time()
  if (verbose) {
    cat("######### Compute multivariate encoding #########\n")
  }

  nBasis <- basisobj$nbasis
  phi <- fd(diag(nBasis), basisobj)

  uniqueId <- as.character(unique(data$id))
  if (is.null(stateColumns)) {
    stateColumns <- sort(colnames(data)[grep("state", colnames(data))])
  }
  K <- sapply(data[stateColumns], FUN = function(x) length(unique(x)))

  nId <- length(uniqueId)

  if (verbose) {
    cat(paste0("Number of individuals: ", nId, "\n"))
    cat(paste0("Number of categorical variables: ", length(K), "\n"))
    cat(paste0("Number of states: ", paste(K, collapse = ", "), "\n"))
    cat(paste0("Basis type: ", basisobj$type, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
  }


  V_multi <- computeVlist(data, phi, K, stateColumns, uniqueId, verbose, nCores)

  U_mean <- computeUmean(data, phi, K, stateColumns, uniqueId, verbose, nCores)

  if (verbose) {
    cat(paste0("\n---- Compute encoding:\n"))
  }

  Fmat_normal <- matrix(0, nrow = nBasis * sum(K), ncol = nBasis * sum(K))

  for (i in seq_along(K)) {
    for (j in seq_along(K)) {
      ind_start_dim1 <- (sum(K) - sum(K[i:length(K)])) * nBasis + 1
      ind_end_dim1 <- sum(K[seq_len(i)]) * nBasis
      ind_start_dim2 <- (sum(K) - sum(K[j:length(K)])) * nBasis + 1
      ind_end_dim2 <- sum(K[seq_len(j)]) * nBasis
      Fmat_normal[ind_start_dim1:ind_end_dim1, ind_start_dim2:ind_end_dim2] <- U_mean[[i]][[j]]
    }
  }

  G <- cov(V_multi)
  V <- V_multi

  if (is.null(epsilon)) {
    epsilon <- min(abs(diag(Fmat_normal)[diag(Fmat_normal) != 0])) * c(1e-12, 1e-9, 1e-6, 1e-4, 1e-2, 1e-1)
    if (verbose) {
      cat(paste0("You did not provide a value for epsilon. Several values will be tested: ", paste(epsilon, collapse = ", ")))
    }
  }

  epsilon <- sort(epsilon)
  isInverted <- FALSE

  for (eps in epsilon) {
    tryCatch(
      {
        Fmat <- Fmat_normal + eps * diag(ncol(Fmat_normal))
        ind0 <- (colSums(Fmat == 0) == nrow(Fmat))
        F05 <- t(mroot(Fmat[!ind0, !ind0])) # F  = t(F05)%*%F05
        invF05 <- solve(F05)
        isInverted <- TRUE
        break
      },
      error = function(e) {
        print(paste0("Failed for epsilon=", eps))
      }
    )
  }

  if (!isInverted) {
    stop("solve.default(F05) : 'F05' must be squared. The epsilon is perharps too big.")
  }

  G <- G[!ind0, !ind0]
  V <- V[, !ind0]

  res <- eigen(t(invF05) %*% G %*% invF05)

  invF05vec <- invF05 %*% res$vectors
  aux1 <- split(invF05vec, rep(seq_len(ncol(res$vectors)), each = nrow(res$vectors)))

  alpha <- lapply(aux1, function(w) {
    wb <- rep(NA, nBasis * sum(K))
    wb[!ind0] <- Re(w)
    name <- paste0(paste0("Dim ", rep(seq_along(K), K), " K"), do.call(c, lapply(K, seq_len)))
    return(matrix(wb, ncol = sum(K), dimnames = list(NULL, name)))
  })

  pc <- V %*% invF05vec
  rownames(pc) <- uniqueId

  # en attendant mieux, besoin pour le plot
  uniqueTime <- sort(unique(data$time))
  pt <- list(
    pt = matrix(
      1 / ncol(alpha[[1]]),
      nrow = ncol(alpha[[1]]),
      ncol = length(uniqueTime),
      dimnames = list(colnames(alpha[[1]]), uniqueTime)
    ),
    t = uniqueTime
  )
  class(pt) <- "pt"

  mult_enc <- list(
    "pc" = pc,
    "alpha" = alpha,
    "eigenvalues" = Re(res$values),
    "F" = Fmat,
    "G" = G,
    "V" = V,
    "basisobj" = basisobj,
    "pt" = pt,
    "epsilon" = eps
  )
  class(mult_enc) <- "fmca"

  t2 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t2 - t1)[3], 2), "s\n"))
  }

  mult_enc$runTime <- as.numeric((t2 - t1)[3])
  if (verbose) {
    cat(paste0("Run Time: ", round(mult_enc$runTime, 2), "s\n"))
  }

  return(mult_enc)
}


#' Convert a list of cfd into a multivariate cfd
#'
#' @param x list of cfd
#' @param state_columns column names for multivariate states. By default, "state1", "state2", ...
#'
#' @return a multivariate cfd
#'
#' @examples
#' set.seed(42)
#' x1 <- generate_Markov(n = 10, K = 2)
#' x1 <- cut_data(x1, Tmax = 1)
#' x2 <- generate_Markov(n = 10, K = 2)
#' x2 <- cut_data(x2, Tmax = 1)
#'
#' x <- list(x1, x2)
#'
#' mvcfd <- convert2mvcfd(x)
#'
#' @export
convert2mvcfd <- function(x, state_columns = NULL) {
  nDim <- length(x)
  if (is.null(state_columns)) {
    state_columns <- paste0("state", seq_len(nDim))
  }

  # rename state columns to each data frame
  for (i in seq_len(nDim)) {
    x[[i]] <- x[[i]] %>% rename(!!state_columns[i] := "state")
  }

  # merge data frames
  x <- Reduce(function(x1, x2) merge(x1, x2, by = c("id", "time"), all = TRUE), x)

  # order by id and time
  x <- arrange(x, .data[["id"]], .data[["time"]])

  # fill missing values with the state before
  x <- as.data.frame(x %>% group_by(.data[["id"]]) %>% fill(all_of(state_columns), .direction = "downup"))

  return(distinct(x))
}
