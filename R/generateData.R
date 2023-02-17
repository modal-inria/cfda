#' Generate Markov Trajectories
#'
#' Simulate individuals from a Markov process defined by a transition matrix,
#' time spent in each time and initial probabilities.
#'
#' @param n number of trajectories to generate
#' @param K number of states
#' @param P matrix containing the transition probabilities from one state to another.
#' Each row contains positive reals summing to 1.
#' @param lambda time spent in each state
#' @param pi0 initial distribution of states
#' @param Tmax maximal duration of trajectories
#' @param labels state names. If \code{NULL}, integers are used
#'
#' @return
#' a data.frame with 3 columns: \code{id}, id of the trajectory, \code{time},
#' time at which a change occurs and \code{state}, new state.
#'
#' @details
#' For one individual, assuming the current state is \eqn{s_j} at time \eqn{t_j},
#' the next state and time is simulated as follows:
#' \enumerate{
#' \item generate one sample, \eqn{d}, of an exponential law of parameter \code{lambda[s_j]}
#' \item define the next time values as: \eqn{t_{j+1} = t_j + d}
#' \item generate the new state \eqn{s_{j+1}} using a multinomial law with probabilities \code{Q[s_j,]}
#' }
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(
#'   n = 100, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10,
#'   labels = c("A", "C", "G", "T")
#' )
#'
#' head(d_JK)
#' @author Cristian Preda
#'
#' @export
generate_Markov <- function(n = 5, K = 2, P = (1 - diag(K)) / (K - 1), lambda = rep(1, K), pi0 = c(1, rep(0, K - 1)), Tmax = 1, labels = NULL) {
  ## check parameters
  if (any(is.na(n)) || (length(n) != 1) || !is.whole.number(n)) {
    stop("n must be a positive integer.")
  }
  if (any(is.na(K)) || (length(K) != 1) || !is.whole.number(K) || (K <= 1)) {
    stop("K must be an integer > 1.")
  }
  if (any(is.na(Tmax)) || !is.numeric(Tmax) || (length(Tmax) != 1) || (Tmax <= 0)) {
    stop("Tmax must be a positive real.")
  }
  if (any(is.na(P)) || !is.matrix(P) || (nrow(P) != K) || (nrow(P) != K) || any(P < 0)) {
    stop("P must be a matrix of size K x K of positive real.")
  }
  if (any(is.na(lambda)) || !is.vector(lambda) || (length(lambda) != K) || any(lambda < 0)) {
    stop("lambda must be a vector of length K of positive real.")
  }
  if (any(is.na(pi0)) || !is.vector(pi0) || (length(pi0) != K) || any(pi0 < 0)) {
    stop("pi0 must be a vector of length K of positive real.")
  }
  if (!is.null(labels)) {
    if (!is.vector(labels) || (length(unique(labels)) != K)) {
      stop("labels must be NULL or a vector of length K.")
    }
  }
  ## end check


  d <- data.frame(id = numeric(0), time = numeric(0), state = numeric(0))

  for (i in seq_len(n)) {
    e <- sample(K, 1, prob = pi0)
    t <- 0
    while (t <= Tmax) {
      d <- rbind(d, data.frame(id = i, time = t, state = e))
      sej <- rexp(1, lambda[e])
      t <- t + sej
      e <- sample(K, 1, prob = P[e, ])
    }
  }

  if (!is.null(labels)) {
    d$state <- labels[d$state]
  }

  row.names(d) <- NULL

  return(d)
}

#' Generate data following a 2 states model
#'
#' Generate individuals such that each individual starts at time 0 with state 0 and then an unique change
#' to state 1 occurs at a time \eqn{t} generated using an uniform law between 0 and 1.
#'
#' @param n number of individuals
#'
#' @return
#' a data.frame with 3 columns: \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, new state.
#'
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
generate_2State <- function(n) {
  ## check parameters
  if (any(is.na(n)) || (length(n) != 1) || !is.whole.number(n)) {
    stop("n must be a positive integer.")
  }
  ## end check

  temps <- rep(0, 2 * n)
  temps[seq_len(n) * 2] <- runif(n)
  d <- data.frame(id = rep(seq_len(n), each = 2), time = temps, state = rep(0:1, n))

  row.names(d) <- NULL

  return(d)
}
