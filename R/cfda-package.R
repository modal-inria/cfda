#' @import diagram fda ggplot2 msm mgcv pbapply
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom stats cov integrate rexp runif var
#' @importFrom graphics boxplot
#' @importFrom utils head
#'
#' @title Categorical Functional Data Analysis
#' @docType package
#' @aliases cfda-package
#' @name cfda-package
#' @description
#' cfda provides functions for the analysis of categorical functional data.
#'
#' The main contribution is the computation of an optimal encoding (real functional variable) of each state of the categorical functional data.
#' This can be done using the \code{\link{compute_optimal_encoding}} function that takes in arguments the data in a specific format and
#' a basis of functions created using the \code{fda} package (cf. \code{\link{create.basis}}). The output can be analysed with \code{\link{summary.fmca}},
#' \code{\link{plot.fmca}}, \code{\link{get_encoding}}, \code{\link{plotEigenvalues}} and \code{\link{plotComponent}}.
#'
#' Moreover, \code{cfda} contains functions to visualize and compute some statistics about categorical functional data.
#' A summary of the dataset is available with \code{\link{summary_cfd}}.
#' \code{\link{plotData}} shows a graphical representation of the dataset.
#' Basic statistics can be computed: the number of jumps (\code{\link{compute_number_jumps}}), the duration (\code{\link{compute_duration}}),
#' the time spent in each state (\code{\link{compute_time_spent}}), the probability to be in each state at any given time (\code{\link{estimate_pt}}),
#' the transition table (\code{\link{statetable}}).
#'
#' The parameters of a Markov process can be estimated using \code{\link{estimate_Markov}} function.
#'
#' In order to test the different functions, a real dataset is provided (\code{\link{biofam2}}) as well as two functions for generating data:
#' (\code{\link{generate_Markov}} and \code{\link{generate_2State}}).
#'
#' @details
#' See the vignette for a detailed example and mathematical background:
#' \code{RShowDoc("cfda", package = "cfda")}
#'
#'
#' @references
#' \itemize{
#'   \item Deville J.C. (1982) Analyse de données chronologiques qualitatives : comment analyser des calendriers ?, Annales de l'INSEE, No 45, p. 45-104.
#'   \item Deville J.C. et  Saporta G. (1980) Analyse harmonique qualitative, DIDAY et al. (editors), Data Analysis and Informatics, North Holland, p. 375-389.
#'   \item Saporta G. (1981) Méthodes exploratoires d'analyse de données temporelles, Cahiers du B.U.R.O, Université Pierre et Marie Curie, 37-38, Paris.
#'   \item Preda C, Grimonprez Q, Vandewalle V. Categorical Functional Data Analysis. The cfda R Package. Mathematics. 2021; 9(23):3074. https://doi.org/10.3390/math9233074
#' }
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' Tmax <- 5
#' PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
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
#' # plot eigenvalues
#' plotEigenvalues(encoding, cumulative = TRUE, normalize = TRUE)
#'
#' # plot the two first components
#' plotComponent(encoding, comp = c(1, 2))
#'
#' # plot the encoding using the first harmonic
#' plot(encoding)
#'
#' # extract the encoding using the first harmonic
#' encod <- get_encoding(encoding)
#'
#'
#' @seealso \link{compute_optimal_encoding}
#'
#' @keywords package
NULL
