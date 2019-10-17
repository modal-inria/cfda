#' Generate Markov Trajectories
#' 
#' @param n number of trajectories to generate
#' @param K number of states
#' @param Q matrix which indicates the allowed transitions in the continuous-time Markov chain, and optionally also the initial values of those transitions
#' @param lambda time spent in each state
#' @param pi_0 initial distribution of states
#' @param Tmax length of trajectories
#' 
#' @return 
#' a data.frame with 3 columns : \code{id}, id of the trajectory, \code{time}, time at chich a change occurs and \code{state}, new state.
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK = generate_Markov_cfd (n = 100, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' head(d_JK)
#' 
#' 
#' @author Cristian Preda
#' @export
generate_Markov_cfd <- function(n = 5, K = 2, Q = 1 - diag(K), lambda = rep(1, K), pi_0 = c(1, rep(0, K-1)), Tmax = 1)
{
  d <- data.frame(id = numeric(0), time = numeric(0), state = numeric(0))
  for(i in 1:n)
  {
    #tirage de l'etat initial
    etat_initial <- sample(1:K, 1, prob = pi_0) 
    e <- etat_initial
    t = 0
    while(t <= Tmax)
    {
      d = rbind(d, data.frame(id = i, time=t, state = e))
      sej = rexp(1, lambda[e]) 
      t = t + sej
      e = sample(1:K, 1, prob = Q[e,])
    }
  }
  
  return(d)
}
