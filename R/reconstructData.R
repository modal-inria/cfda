#' @title Reconstruct the indicators using encoding
#'
#' @description
#' The reconstruction formula is:
#' \deqn{1^{x}(t) = p^x(t) ( 1 + \sum_{i\geq 1} z_i*a_i^x(t)})
#'
#' with \eqn{z_i}, the i-th principal component,
#' encoding \eqn{a_i^x = \sum_j \alpha_{(x, j)} * \phi_j(t)}
#' and \eqn{p^x(t) =  1 / (\sum_{i \geq 1} a_i^x(t)^2)}
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param nComp number of components to use for the reconstruction. By default, all are used.
#' @param timeValues vector containing time values at which compute the indicators. If NULL, the time values from the data
#' @param propMinEigenvalues Only if nComp = NULL. Minimal proportion used to estimate the number of non-null eigenvalues
#'
#' @return a data.frame with columns: time, id, state1, ..., stateK, state.
#' state1 contains the estimated indicator values for the first state.
#' state contains the state with the maximum values of all indicators
#'
#' @examples
#' set.seed(42)
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 3
#' Tmax <- 1
#' d_JK <- generate_Markov(n = 100, K = K, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 20
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' \donttest{
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#'
#' indicators <- reconstructIndicators(encoding)
#'
#' # we plot the first path and its reconstructed indicators
#' iInd <- 3
#' plotData(d_JK2[d_JK2$id == iInd, ])
#'
#' plotIndicatorsReconstruction(indicators, id = iInd)
#'
#' # the column state contains the state associated with the greatest indicator.
#' # So, the output can be used with plotData function
#' plotData(remove_duplicated_states(indicators[indicators$id == iInd, ]))
#' }
#' @author Quentin Grimonprez
#' @seealso \code{\link{plotIndicatorsReconstruction}}
#' @export
reconstructIndicators <- function(x, nComp = NULL, timeValues = NULL, propMinEigenvalues = 1e-4) {
  checkFmca(x)
  if (
    any(is.na(propMinEigenvalues)) || !is.numeric(propMinEigenvalues) || (length(propMinEigenvalues) != 1) || (
      propMinEigenvalues > 10) || (propMinEigenvalues <= 0)
  ) {
    stop("propMinEigenvalues must be a positive real number.")
  }
  if (is.null(nComp)) {
    nComp <- sum(x$eigenvalues / sum(x$eigenvalues) > propMinEigenvalues)
    print(paste0("Reconstruct data using ", nComp, " components (out of ", length(x$eigenvalues), ")"))
  }
  checkInteger(
    nComp,
    minValue = 1, maxValue = length(x$alpha), minEqual = TRUE, maxEqual = TRUE,
    customMessage = "nComp must be an integer between 1 and the number of components."
  )

  if (is.null(timeValues)) {
    timeValues <- x$pt$t
  }

  # a_i^x : encoding = \sum_j \alpha_{(x, j)} * phi_j(t)
  a_i_t_x <- list()
  for (j in seq_len(nComp)) {
    a <- fda::fd(x$alpha[[j]], x$basisobj)
    a_i_t_x[[j]] <- fda::eval.fd(timeValues, a)
  }

  # Compute p^x(t) =  1 / (1 + sum_{i \geq 1} a_i^x(t)^2 * lambda_i)
  denominator <- 1
  for (j in seq_len(nComp)) {
    denominator <- denominator + a_i_t_x[[j]]^2 * x$eigenvalues[j]
  }
  pt <- 1 / denominator

  allDf <- list()
  # 1^{x}(t) = p^x(t) * (1 + \sum_{i\geq 1} z_i*a_i^x(t))
  z <- x$pc[, seq_len(nComp), drop = FALSE]
  for (iInd in seq_len(nrow(z))) {
    y <- 1
    for (j in seq_len(nComp)) {
      y <- y + z[iInd, j] * a_i_t_x[[j]]
    }

    yy <- y * pt

    df <- data.frame("time" = timeValues, "id" = rownames(z)[iInd])
    for (i in seq_along(x$label$label)) {
      df[[paste0("state", x$label$label[i])]] <- yy[, i]
    }
    allDf[[iInd]] <- df
    allDf[[iInd]]$state <- x$label$label[apply(allDf[[iInd]][paste0("state", x$label$label)], 1, which.max)]
  }

  return(do.call(rbind, allDf))
  # return(list(reconstruction = do.call(rbind, allDf), pt = pt, t = timeValues))
}

#' @title Plot reconstructed indicators
#'
#' @param reconstruction output of \code{\link{reconstructIndicators}}
#' @param id id of the individual to plot. \code{id} must be in \code{reconstruction$id}
#'
#' @return ggplot
#'
#' @inherit reconstructIndicators examples
#'
#' @author Quentin Grimonprez
#' @seealso \code{\link{reconstructIndicators}}
#' @export
plotIndicatorsReconstruction <- function(reconstruction, id) {
  x <- reconstruction[reconstruction$id == id, ]
  x$id <- NULL
  x$state <- NULL

  p <- ggplot(
    data = pivot_longer(x, cols = -1, values_to = "probability", names_to = "state"),
    aes(x = .data$time, y = .data$probability, color = .data$state)
  ) +
    geom_line()
  return(p)
}


# functional object (p^x)-1
invProba <- function(x, nComp = NULL) {
  if (is.null(nComp)) {
    nComp <- sum(x$eigenvalues / sum(x$eigenvalues) > 1e-6)
  }
  # a_i^x : encoding = \sum_j \alpha_{(x, j)} * phi_j(t)
  a <- list()
  for (j in seq_len(nComp)) {
    a[[j]] <- fda::fd(x$alpha[[j]], x$basisobj)
  }

  # Compute p^x(t) =  1 / (1 + sum_{i \geq 1} a_i^x(t)^2 * lambda_i)
  denominator <- 1
  for (j in seq_len(nComp)) {
    denominator <- denominator + a[[j]]^2 * x$eigenvalues[j]
  }

  return(denominator)
}
