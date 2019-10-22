#' Compute time spent in each state
#'
#' For each individual, compute the time spent in each state
#'
#' @param data_msm data.frame containing \code{id}, \code{time} and \code{state} (see \code{\link{generate_Markov_cfd}})
#' @param K maximal number of state
#'
#' @return a vector of length \code{K} containing the total time spent in each state
#'
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' # cut at Tmax = 8
#' d_JK2 <- msm2msmTmax(d_JK, Tmax = 8)
#'
#' #  compute time spent by each id in each state
#' timeSpent <- compute_Time_Spent(d_JK2, K)
#' 
#' 
#' @author Cristian Preda
#'
#' @export
compute_Time_Spent <- function(data_msm, K)
{
  res <- by(data_msm, data_msm$id, function(x){compute_Time_Spent_intern(x, K)})
  out <- do.call(rbind, res)
  colnames(out) = 1:K
  
  return(out)
}
  

# combien de temps passe un id dans [0,T] dans chaque etat, x vient d'un msmT
compute_Time_Spent_intern <- function(data_msm, K)
{
  aux <- rep(0, K)            
  for(state in 1:K)
  {
    idx <- which(data_msm$state == state)
    for(u in idx)
    {
      if(u < nrow(data_msm))
        aux[state] = aux[state] + data_msm[u+1, "time"] - data_msm[u, "time"]
    }  
  }
  return(aux)
}          

