
#' Plot the optimal encoding
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param harm harmonic to use for the encoding
#' @param states states to plot (default = NULL, it plots all states)
#' @param addCI if TRUE, plot confidence interval (only when \code{computeCI = TRUE} in \link{compute_optimal_encoding})
#' @param coeff the confidence interval is computed with +- coeff * the standard deviation
#' @param col a vector containing color for each state
#' @param nx number of time points used to plot
#' @param ... not used
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @details
#' The encoding for the harmonic \code{h} is \eqn{a_{x}^{(h)} \approx \sum_{i=1}^m \alpha_{x,i}^{(h)}\phi_i}.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 6
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' \donttest{
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#'
#' # plot the encoding produced by the first harmonic
#' plot(encoding)
#'
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plot(encoding, harm = 2, col = c("red", "blue", "darkgreen", "yellow")) +
#'   labs(title = "Optimal encoding")
#' }
#'
#' @family encoding functions
#' @author Quentin Grimonprez
#'
#' @export
plot.fmca <- function(x, harm = 1, states = NULL, addCI = FALSE, coeff = 1.96, col = NULL, nx = 128, ...) {
  checkLogical(addCI, "addCI")

  if (any(is.na(coeff)) || (length(coeff) > 1) || !is.numeric(coeff) || (coeff < 0)) {
    stop("coeff must be a positive real.")
  }

  fdmat <- get_encoding(x, harm = harm, fdObject = FALSE, nx = nx) # harm and nx are checked in this function

  if (is.null(states)) {
    states <- colnames(fdmat$y)
  }

  states <- intersect(states, colnames(fdmat$y))
  if (length(states) == 0) {
    stop("No correct states given.")
  }

  if (addCI && ("bootstrap" %in% names(x))) {
    variance <- computeVarianceEncoding(x$varAlpha, x$basisobj, harm = harm, nx = nx)
    p <- plotEncodingCI(fdmat, variance, coeff, states, harm, col)
  } else {
    p <- plotEncoding(fdmat, states, harm, col)
  }

  return(p)
}

# plot the encoding and the associated confidence interval for each state
# @author Quentin Grimonprez
plotEncodingCI <- function(fdmat, variance, coeff = 2, states = NULL, harm = 1, col = NULL) {
  p <- ggplot()
  for (i in match(states, colnames(fdmat$y))) {
    df <- data.frame(
      time = fdmat$x,
      ymin = fdmat$y[, i] - sqrt(variance[[i]]) * coeff,
      ymax = fdmat$y[, i] + sqrt(variance[[i]]) * coeff,
      State = factor(rep(colnames(fdmat$y)[i], each = nrow(fdmat$y)), levels = colnames(fdmat$y))
    )
    p <- p + geom_ribbon(
      data = df, aes_string(ymin = "ymin", ymax = "ymax", x = "time", fill = "State"),
      colour = NA, alpha = 0.8
    )
  }

  df <- data.frame(
    x = rep(fdmat$x, ncol(fdmat$y)),
    y = as.vector(fdmat$y),
    State = factor(rep(colnames(fdmat$y), each = nrow(fdmat$y)), levels = colnames(fdmat$y))
  )
  df <- df[df$State %in% states, ]

  p <- p +
    geom_line(data = df, mapping = aes_string(x = "x", y = "y", group = "State", colour = "State"), alpha = 1) +
    scale_colour_hue(l = 30, drop = FALSE)

  p <- p +
    labs(x = "Time", y = expression(paste("a"["x"], "(t)")), title = paste0("Encoding function for harmonic number ", harm))

  if (!is.null(col)) {
    p <- p + scale_fill_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData

  return(p)
}

# plot the encoding for each state
# @author Quentin Grimonprez
plotEncoding <- function(fdmat, states = NULL, harm = 1, col = NULL) {
  df <- data.frame(
    x = rep(fdmat$x, ncol(fdmat$y)),
    y = as.vector(fdmat$y),
    State = factor(rep(colnames(fdmat$y), each = nrow(fdmat$y)), levels = colnames(fdmat$y))
  )

  df <- df[df$State %in% states, ]
  p <- ggplot(df, aes_string(x = "x", y = "y", group = "State", colour = "State")) +
    geom_line()

  p <- p +
    labs(x = "Time", y = expression(paste("a"["x"], "(t)")), title = paste0("Encoding function for harmonic number ", harm))

  if (!is.null(col)) {
    p <- p + scale_colour_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_colour_hue(drop = FALSE)
  } # keep the same color order as plotData

  return(p)
}



#' Extract the computed encoding
#'
#' Extract the encoding as an \code{fd} object or as a matrix
#'
#' @param x Output of \code{\link{compute_optimal_encoding}}
#' @param harm harmonic to use for the encoding
#' @param fdObject If TRUE returns a \code{fd} object else a matrix
#' @param nx (Only if \code{fdObject = TRUE}) Number of points to evaluate the encoding
#'
#' @return a \code{fd} object or a list of two elements \code{y}, a matrix with \code{nx} rows containing
#' the encoding of the state and \code{x}, the vector with time values.
#'
#' @details
#' The encoding is \eqn{a_{x} \approx \sum_{i=1}^m \alpha_{x,i}\phi_i}.
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 6
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' \donttest{
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#'
#' # extract the encoding using 1 harmonic
#' encodFd <- get_encoding(encoding, fdObject = TRUE)
#' encodMat <- get_encoding(encoding, nx = 200)
#' }
#'
#' @author Cristian Preda
#' @family encoding functions
#' @export
get_encoding <- function(x, harm = 1, fdObject = FALSE, nx = NULL) {
  ## check parameters
  if (!inherits(x, "fmca")) {
    stop("x must be a fmca object.")
  }
  checkLogical(fdObject, "fdObject")
  if (!is.null(nx)) {
    if (any(is.na(nx)) || (length(nx) > 1) || !is.whole.number(nx) || (nx <= 0)) {
      stop("nx must be a positive integer.")
    }
  }
  if (any(is.na(harm)) || (length(harm) > 1) || !is.whole.number(harm) || (harm < 1) || (harm > length(x$alpha))) {
    stop("harm must be an integer between 1 and the number of components.")
  }
  ##


  alpha <- x$alpha[[harm]]

  if (fdObject) {
    fdObj <- fd(x$alpha[[harm]], x$basisobj)

    return(fdObj)
  } else {
    alpha[is.na(alpha)] <- 0
    fdObj <- fd(alpha, x$basisobj)

    rangex <- fdObj$basis$rangeval
    nBasis <- fdObj$basis$nbasis

    if (is.null(nx)) {
      nx <- max(length(x$pt$t), 128, 10 * nBasis + 1)
    }

    timeVal <- seq(rangex[1], rangex[2], length = nx)

    fdmat <- eval.fd(timeVal, fdObj)

    fdmat <- removeTimeAssociatedWithNACoeff(fdmat, timeVal, x$pt)

    return(list(x = timeVal, y = fdmat))
  }
}


# when the probability at a given time is 0, the encoding at this time is returned as NA
removeTimeAssociatedWithNACoeff <- function(fdmat, timeVal, pt) {
  p <- t(sapply(timeVal, get_proba, pt = pt))
  p <- p[, match(colnames(fdmat), colnames(p))]
  p[p != 0] <- 1
  p[p == 0] <- NA

  return(p * fdmat)
}


#' Plot Components
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param comp a vector of two elements indicating the components to plot
#' @param addNames if TRUE, add the id labels on the plot
#' @param nudge_x,nudge_y horizontal and vertical adjustment to nudge labels by
#' @param size size of labels
#' @param ... \code{geom_point} parameters
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 6
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' \donttest{
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#'
#' plotComponent(encoding, comp = c(1, 2))
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plotComponent(encoding, comp = c(1, 2), shape = 23) +
#'   labs(title = "Two first components")
#' }
#'
#'
#' @author Quentin Grimonprez
#' @family encoding functions
#' @export
plotComponent <- function(x, comp = c(1, 2), addNames = TRUE, nudge_x = 0.1, nudge_y = 0.1, size = 4, ...) {
  ## check parameters
  if (!inherits(x, "fmca")) {
    stop("x must be a fmca object.")
  }
  checkLogical(addNames, "addNames")
  if (length(comp) != 2) {
    stop("comp must be a vector of positive integers of length 2.")
  }
  if (any(is.na(comp)) || any(!is.whole.number(comp)) || any(comp < 0) || any(comp > ncol(x$pc))) {
    stop("comp must be a vector of positive integers of length 2.")
  }
  ##

  df <- as.data.frame(Re(x$pc))
  df$name <- rownames(x$pc)

  p <- ggplot(df, aes_string(x = paste0("V", comp[1]), y = paste0("V", comp[2]))) +
    geom_point(...) +
    labs(x = paste0("Comp ", comp[1]), y = paste0("Comp ", comp[2]))

  if (addNames) {
    p <- p + geom_text(aes_string(label = "name"), nudge_x = nudge_x, nudge_y = nudge_y, size = size)
  }

  p
}


#' Plot Eigenvalues
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param cumulative if TRUE, plot the cumulative eigenvalues
#' @param normalize if TRUE eigenvalues are normalized for summing to 1
#' @param ... \code{geom_point} parameters
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 6
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' \donttest{
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#'
#' # plot eigenvalues
#' plotEigenvalues(encoding, cumulative = TRUE, normalize = TRUE)
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plotEigenvalues(encoding, shape = 23) +
#'   labs(caption = "Jukes-Cantor model of nucleotide replacement")
#' }
#'
#'
#' @author Quentin Grimonprez
#' @family encoding functions
#' @export
plotEigenvalues <- function(x, cumulative = FALSE, normalize = FALSE, ...) {
  ## check parameters
  if (!inherits(x, "fmca")) {
    stop("x must be a fmca object.")
  }
  checkLogical(cumulative, "cumulative")
  checkLogical(normalize, "normalize")
  ##

  if (normalize) {
    eigenv <- x$eigenvalues / sum(x$eigenvalues)
  } else {
    eigenv <- x$eigenvalues
  }

  comp <- seq_along(eigenv)
  if (cumulative) {
    eigenv <- c(0, cumsum(eigenv))
    comp <- seq_along(eigenv) - 1
  }


  df <- data.frame(eigenvalues = eigenv, component = comp)

  p <- ggplot(df, aes_string(x = "component", y = "eigenvalues")) +
    geom_point(...) +
    geom_step() +
    labs(
      title = ifelse(cumulative, "Cumulative eigenvalues", "Eigenvalues"),
      x = ifelse(cumulative, "Number of components", "Components"),
      y = ifelse(cumulative, "Cumulative eigenvalues", "Eigenvalues")
    )

  p
}
