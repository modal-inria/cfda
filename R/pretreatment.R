#' Cut data to a maximal given time
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state (integer starting at 1).
#' @param Tmax max time considered
#' 
#' @return a data.frame with the same format as \code{data_msm} where each individual has \code{Tmax} as last time entry.
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK = generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' tail(d_JK)
#' 
#' # cut at Tmax = 8
#' d_JK2 <- msm2msmTmax(d_JK, Tmax = 8)
#' tail(d_JK2)
#' 
#' @author Cristian Preda
#' 
#' @export
msm2msmTmax <- function(data_msm, Tmax)
{
  ## check parameters
  checkDataMsm(data_msm)
  if(!is.numeric(Tmax) || (length(Tmax) != 1))
    stop("Tmax must be a real.")
  ## end check
  
  do.call(rbind, by(data_msm, data_msm$id, function(x){cut_cfd(x, Tmax)}))
}


cut_cfd <- function(data_msm, Tmax)
{
  l <- nrow(data_msm)
  currTmax <- max(data_msm$time)
  if(Tmax > currTmax) 
  {
    return(rbind(data_msm, data.frame(id = data_msm[1, "id"], state = data_msm[l, "state"], time = Tmax)))
  }
  else
  {
    if(currTmax == Tmax)
    {
     return(data_msm) 
    }
    else
    {
      k <- max(which(data_msm$time <= Tmax))  
      return(rbind(data_msm[1:k,], data.frame(state = data_msm[k, "state"], time = Tmax, id = data_msm[1, "id"])))  
    }
    
  }
  
}
