#' @import diagram doSNOW fda ggplot2 msm mgcv
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster stopCluster
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom stats cov integrate rexp runif
#' @importFrom graphics boxplot
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @title Categorical Functional Data Analysis
#' @docType package
#' @aliases cfda-package
#' @name cfda-package
#' @description  
#' Categorical Functional Data Analysis
#' 
#' @details 
#' See the vignette for a detailled example and mathematical background:
#' \code{RShowDoc("cfda", package = "cfda")}
#' 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, Tmax)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # plot the optimal encoding
#' plot(encoding)
#' 
#' # plot the two first components
#' plotComponent(encoding, comp = c(1, 2))
#' 
#' # extract the optimal encoding
#' getEncoding(encoding)
#' 
#' 
#' @seealso \link{compute_optimal_encoding}
#' 
#' @keywords package
NULL
