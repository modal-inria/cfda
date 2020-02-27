

#' Plot the optimal encoding
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param harm harmonic to use for the encoding
#' @param col a vector containing color for each state.
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
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding 
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # plot the encoding produced by the first harmonic
#' plot(encoding)
#' 
#' 
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plot(encoding, harm = 2, col = c("red" , "blue", "darkgreen", "yellow")) +
#'    labs(title = "Optimal encoding")
#' 
#' @seealso \link{plotComponent} \link{plotEigenvalues}
#' 
#' @author Quentin Grimonprez
#' 
#' @export
plot.fmca <- function(x, harm = 1, col = NULL, ...)
{
  fdmat <- get_encoding(x, harm = harm, fdObject = FALSE)
  df <- data.frame(x = rep(fdmat$x, ncol(fdmat$y)), y = as.vector(fdmat$y), State = factor(rep(colnames(fdmat$y), each = nrow(fdmat$y)), levels = colnames(fdmat$y)))
  
  p <- ggplot(df, aes_string(x = "x", y = "y", group = "State", colour = "State")) +
    geom_line() +
    labs(x = "Time", y = expression(paste("a"["x"], "(t)")), title = paste0("Encoding with harmonic number ", harm))
  
  if(!is.null(col))
    p = p + scale_colour_manual(values = col)
  
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
#' @return a \code{fd} object or a list of two elements \code{y}, a matrix with \code{nx} rows containing the encoding of the state and \code{x}, the vector with time values.
#' 
#' @details 
#' The encoding is \eqn{a_{x} \approx \sum_{i=1}^m \alpha_{x,i}\phi_i}.
#' 
#' 
#' @examples 
#' # Simulate the Jukes-Cantor model of nucleotide replacement 
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # extract the encoding using 1 harmonic
#' encodFd <- get_encoding(encoding, fdObject = TRUE)
#' encodMat <- get_encoding(encoding, nx = 200)
#' 
#' @author Cristian Preda
#'
#' @export
get_encoding <- function(x, harm = 1, fdObject = FALSE, nx = NULL)
{
  ## check parameters
  if(class(x) != "fmca")
    stop("x must be a fmca object.")
  checkLogical(fdObject, "fdObject")
  if(!is.null(nx))
  {
    if(any(is.na(nx)) || (length(nx) > 1) || !is.whole.number(nx) || (nx <= 0))
      stop("nx must be a positive integer.")
  }
  if(any(is.na(harm)) || (length(harm) > 1) || !is.whole.number(harm) || (harm < 1) || (harm > length(x$alpha)))
    stop("harm must be an integer between 1 and the number of components.")
  ##
  
  fdObj <- fd(x$alpha[[harm]], x$basisobj)
  
  
  if(fdObject)
  {
    return(fdObj)
  }else{
    rangex <- fdObj$basis$rangeval
    nBasis <- fdObj$basis$nbasis
    
    if(is.null(nx))
      nx = max(501, 10 * nBasis + 1)
    
    timeVal <- seq(rangex[1], rangex[2], length = nx)
    
    fdmat <- eval.fd(timeVal, fdObj)
    
    return(list(x = timeVal, y = fdmat))
  }
  
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
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' plotComponent(encoding, comp = c(1, 2))
#' 
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plotComponent(encoding, comp = c(1, 2), shape = 23) +
#'    labs(title = "Two first components")
#' 
#' @seealso \link{plot.fmca} \link{plotEigenvalues}
#' 
#' @author Quentin Grimonprez
#' 
#' @export
plotComponent <- function(x, comp = c(1, 2), addNames = TRUE, nudge_x = 0.1, nudge_y = 0.1, size = 4, ...)
{
  ## check parameters
  if(class(x) != "fmca")
    stop("x must be a fmca object.")
  checkLogical(addNames, "addNames")
  if(length(comp) != 2)
    stop("comp must be a vector of positive integers of length 2.")
  if(any(is.na(comp)) || any(!is.whole.number(comp)) || any(comp < 0) || any(comp > ncol(x$pc)))
    stop("comp must be a vector of positive integers of length 2.")
  ##
  
  df <- as.data.frame(x$pc)
  df$name = rownames(x$pc)
  
  p <- ggplot(df, aes_string(x = paste0("V", comp[1]), y = paste0("V", comp[2]))) +
    geom_point(...) + 
    labs(x = paste0("Comp ", comp[1]), y = paste0("Comp ", comp[2]))
  
  if(addNames)
    p = p + geom_text(aes_string(label = "name"), nudge_x = nudge_x, nudge_y = nudge_y, size = size)
  
  p
}


#' Plot Eigenvalues
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param cumulative if TRUE, plot the cumualtive eigenvalues
#' @param normalize if TRUE eigenvalues are normalized for summing to 1
#' @param ... \code{geom_point} parameters
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package. 
#' 
#' @examples 
#' # Simulate the Jukes-Cantor model of nucleotide replacement  
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # plot eigenvalues
#' plotEigenvalues(encoding, cumulative = TRUE, normalize = TRUE)
#' 
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plotEigenvalues(encoding, shape = 23) +
#'    labs(caption = "Jukes-Cantor model of nucleotide replacement")
#' 
#' @seealso \link{plot.fmca} \link{plotComponent}
#' 
#' @author Quentin Grimonprez
#' 
#' @export
plotEigenvalues <- function(x, cumulative = FALSE, normalize = FALSE, ...)
{
  ## check parameters
  if(class(x) != "fmca")
    stop("x must be a fmca object.")
  checkLogical(cumulative, "cumulative")
  checkLogical(normalize, "normalize")
  ##
  
  if(normalize)
    eigenv <- x$eigenvalues/sum(x$eigenvalues)
  else
    eigenv <- x$eigenvalues
  
  comp <- seq_along(eigenv)
  if(cumulative)
  {
    eigenv = c(0, cumsum(eigenv))
    comp = seq_along(eigenv)-1
  }
  
  
  df <- data.frame(eigenvalues = eigenv, component = comp)
  
  p <- ggplot(df, aes_string(x = "component", y = "eigenvalues")) +
    geom_point(...) + geom_step() +
    labs(title = ifelse(cumulative, "Cumulative eigenvalues", "Eigenvalues"), 
         x = ifelse(cumulative, "Number of components", "Components"),
         y = ifelse(cumulative, "Cumulative eigenvalues", "Eigenvalues"))
  
  p
}
