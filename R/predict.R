#' @title Predict the principal components for new trajectories
#'
#' @description Predict the principal components for new trajectories
#'
#' @param object output of \link{compute_optimal_encoding} function.
#' @param newdata data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state. All individuals must begin at the same time T0 and end at the same time Tmax
#' (use \code{\link{cut_data}}).
#' @param method computation method: "parallel" or "precompute": precompute all integrals
#' (efficient when the number of unique time values is low)
#' @param verbose if TRUE print some information
#' @param nCores number of cores used for parallelization (only if method == "parallel"). Default is half the cores.
#' @param ... parameters for \code{\link{integrate}} function (see details).
#'
#' @return principal components for the individuals
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 6
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(
#'   n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax,
#'   labels = c("A", "C", "G", "T")
#' )
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' \donttest{
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#'
#' # predict principal components
#' d_JK_predict <- generate_Markov(
#'   n = 5, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax,
#'   labels = c("A", "C", "G", "T")
#' )
#' d_JK_predict2 <- cut_data(d_JK, Tmax)
#'
#' pc <- predict(encoding, d_JK_predict2, nCores = 1)
#' }
#'
#' @method predict fmca
#'
#' @family encoding functions
#' @author Quentin Grimonprez
#' @export
predict.fmca <- function(
  object, newdata = NULL, method = c("precompute", "parallel"), verbose = TRUE,
  nCores = max(1, ceiling(detectCores() / 2)), ...
  ) {
  if (is.null(newdata)) {
    return(object$pc)
  }

  ## check parameters
  checkData(newdata)
  checkLogical(verbose, "verbose")
  if (any(is.na(nCores)) || !is.whole.number(nCores) || (nCores < 1)) {
    stop("nCores must be an integer > 0.")
  }
  method <- match.arg(method)
  ##

  if (verbose) {
    cat("######### Predict Principal Components #########\n")
  }

  nCores <- min(max(1, nCores), detectCores() - 1)

  # change state as integer
  newdata$state <- refactorCategorical(newdata$state, object$label$label, object$label$code)

  uniqueId <- as.character(unique(newdata$id))

  K <- length(object$label$label)

  if (method == "precompute") {
    uniqueTime <- sort(unique(newdata$time))
    V <- computeVmatrix2(newdata, object$basisobj, K, uniqueId, uniqueTime, nCores, verbose, ...)
  } else {
    V <- computeVmatrix(newdata, object$basisobj, K, uniqueId, nCores, verbose, ...)
  }

  invF05vec <- sapply(object$alpha, as.vector)
  invF05vec[is.na(invF05vec)] <- 0

  pc <- V %*% invF05vec

  return(pc)
}
