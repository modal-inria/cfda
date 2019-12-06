#' Generate Markov Trajectories
#' 
#' Simulate individuals from a Markov process defined by a transition matrix, time spent in each time and initial probabilities.
#' 
#' @param n number of trajectories to generate
#' @param K number of states
#' @param Q matrix which indicates the allowed transitions in the continuous-time Markov chain, and optionally also the initial values of those transitions
#' @param lambda time spent in each state
#' @param pi0 initial distribution of states
#' @param Tmax maximal duration of trajectories
#' @param labels state names. If \code{NULL}, integers are used
#' 
#' @return 
#' a data.frame with 3 columns: \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, new state.
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 100, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10, 
#'                             labels = c("A", "C", "G", "T"))
#' 
#' head(d_JK)
#' 
#' 
#' @author Cristian Preda
#' @export
generate_Markov <- function(n = 5, K = 2, Q = 1 - diag(K), lambda = rep(1, K), pi0 = c(1, rep(0, K-1)), Tmax = 1, labels = NULL)
{
  ## check parameters
  if((length(n) != 1) || !is.whole.number(n))
    stop("n must be a positive integer.")
  if((length(K) != 1) || !is.whole.number(K) || (K <= 1))
    stop("K must be an integer > 1.")
  if(!is.numeric(Tmax) || (length(Tmax) != 1))
    stop("Tmax must be a real.")
  if(!is.matrix(Q) || (nrow(Q) != K) || (nrow(Q) != K) || any(Q < 0))
    stop("Q must be a matrix of size K x K of positive real.")
  if(!is.vector(lambda) || (length(lambda) != K) || any(lambda < 0))
     stop("lambda must be a vector of length K of positive real.")
  if(!is.vector(pi0) || (length(pi0) != K) || any(pi0 < 0))
    stop("pi0 must be a vector of length K of positive real.")
  if(!is.null(labels))
  {
    if(!is.vector(labels) || (length(labels) != K))
      stop("labels must be NULL or a vector of length K.")  
  }
  ## end check
  
  
  d <- data.frame(id = numeric(0), time = numeric(0), state = numeric(0))
  
  for(i in 1:n)
  {
    e <- sample(K, 1, prob = pi0) 
    t <- 0
    while(t <= Tmax)
    {
      d = rbind(d, data.frame(id = i, time = t, state = e))
      sej = rexp(1, lambda[e]) 
      t = t + sej
      e = sample(K, 1, prob = Q[e,])
    }
  }
  
  if(!is.null(labels))
    d$state = labels[d$state]
  
  return(d)
}

#' Generate data following a 2 states model
#'
#' @param n number of individuals
#'
#' @return 
#' a data.frame with 3 columns: \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, new state.
#'
#' @details 
#' Let \eqn{\theta\sim \mathcal{U}[0, 1]}
#'
#' The state at time t is defined by \eqn{X_t(w) = 1} if \eqn{t < \theta(w)}, 2 otherwise.
#'
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
generate_2State <- function(n)
{
  ## check parameters
  if((length(n) != 1) || !is.whole.number(n))
    stop("n must be a positive integer.")
  ## end check
  
  temps <- rep(0, 2*n)
  temps[(1:n*2)] = runif(n)
  d <- data.frame(id = rep(1:n, each = 2), time = temps,  state = rep(1:2, n))
  
  return(d)
}