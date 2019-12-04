#' Cut data to a maximal given time
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
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
#' d_JK2 <- cut_data(d_JK, Tmax = 8)
#' tail(d_JK2)
#' 
#' @author Cristian Preda
#' 
#' @export
cut_data <- function(data_msm, Tmax)
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
    return(rbind(data_msm, data.frame(id = data_msm$id[1], state = data_msm$state[l], time = Tmax)))
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
      return(rbind(data_msm[1:k,], data.frame(state = data_msm$state[k], time = Tmax, id = data_msm$id[1])))  
    }
    
  }
  
}

# change the labels into integer
#
# @param state vector with labels
# @return a lsit with state containing the new formatted state and label, 
# a data.frame containing the labels and corresponding integers
# 
# @author Quentin Grimonprez
stateToInteger <- function(state)
{
  lab <- data.frame(label = sort(unique(state)), code = 1:length(unique(state)))
  
  newstate <- refactorCategorical(state, lab$label, lab$code)
  
  return(list(state = newstate, label = lab))
}


# Rename a categorical value
#
# @param data matrix/data.frame/vector containing the data
# @param oldCateg vector containing categories to change
# @param newCateg vector containing new categorical values
#
# @return Data with new categorical values
#
# @examples
# dat <- c("single", "married", "married", "divorced", "single")
# refactorCategorical(dat, c("single", "married", "divorced"), 1:3)
#
# @author Quentin Grimonprez
refactorCategorical <- function(data, oldCateg = unique(data), newCateg = 1:length(oldCateg))
{
  ind <- match(data, oldCateg)
  
  if(any(is.na(ind[!is.na(data)])))
    warning("NA produced.")
  
  return(newCateg[ind])
}


