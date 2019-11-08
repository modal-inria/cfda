
#' Estimate transition matrix and spent time
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state (integer starting at 1).
#'
#' @return list of two elements: \code{Q}, the estimated transition matrix, adn \code{lambda}, the estimated time spent in each state
#'
#'
#' @examples
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 100, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' # estimation  
#' mark <- estimation_Markov(d_JK)
#' mark$Q
#' mark$lambda
#'
#' @seealso \link{plot.Markov}
#'
#' @author Cristian Preda
#' 
#' @export
estimation_Markov <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  aux <- statetable.msm(data_msm$state, data_msm$id) # il se peut que la matrice aux ne soit pas carré si au moins un état absorbant existe.

  # identifier les etats d'ou on part jamais dans [0,T] (absorbants ou pas)
  # add row for the missing state
  aux1 <- completeStatetable(aux)

  
  matA <- crudeinits.msm(data_msm$state ~ data_msm$time, subject = data_msm$id, 
                         qmatrix = matrix(as.numeric(aux1 > 0), ncol = length(unique(data_msm$state))))
  
  # estimation of the time spent in each state
  lambda_est <- -diag(matA) 
  
  # estimation of the transition matrix 
  Q_est <- diag(1/lambda_est)%*%(matA + diag(lambda_est))
  colnames(Q_est) = rownames(Q_est) = colnames(aux1)

  out <- list(Q = Q_est, lambda = lambda_est)
  class(out) = "Markov"
  
  return(out)  
}


# identifier les etats d'ou on part jamais dans [0,T] (absorbants ou pas)
# add row (containing 0s) for the missing state
#
# @param aux output of statetable.msm
completeStatetable <- function(aux)
{
  aux1 <- aux
  for(i in which(!as.numeric(colnames(aux))%in%as.numeric(row.names(aux))))
  {
    if(i == 1) 
    {
      aux1 =  rbind(rep(0, ncol(aux)), aux1)
    }
    else if(i == ncol(aux))
    {
      aux1 = rbind(aux1, rep(0, ncol(aux)))
    }
    else 
    {
      aux1 = rbind(aux1[1:(i-1), ], rep(0, ncol(aux)), aux1[i:nrow(aux1), ])  
    }
  }
  
  rownames(aux1) = colnames(aux)
  
  return(aux1)
}


#' Plot the transition graph
#' 
#' Plot the transition graph between the different states. A node corresponds to a state with the mean 
#' time spent in this state. Each arrow represents the probability of transtion between states.
#'
#' @param x output of \code{\link{estimation_Markov}} function
#' @param ... not used
#' 
#' @examples
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 100, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' # estimation  
#' mark <- estimation_Markov(d_JK)
#' 
#' plot(mark)
#'
#' @author Cristian Preda
#' 
#' @export
plot.Markov <- function(x, ...)
{ 
  plotmat(t(round(x$Q, 2)), main = "The transition graph", box.prop = 0.3, box.type = "circle", 
          box.col = "yellow", relsize = 0.9, arr.length = 0.2, lcol = "black", dtext = 0.3, arr.pos = 0.5,
          name = paste(colnames(x$Q), rep("(", ncol(x$Q)), round(1/x$lambda, 2),
                     rep(")", ncol(x$Q)), sep = ""))
}
