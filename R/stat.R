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




#' Extract the state of each individual at a given time
#' 
#' @param data_msm data.frame containing \code{id}, \code{time} and \code{state} (see \code{\link{generate_Markov_cfd}})
#' @param t real
#' 
#' @return a vector containing the state of each individual at time t 
#' 
#' 
#' @examples
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' get_state(d_JK, 6)
#' 
#' @author Cristian Preda
#' 
#' @export
get_state <- function(data_msm, t) 
{
  out <- by(data_msm, data_msm$id, function(x){id_get_state(x, t)})
  out2 <- as.vector(out)
  names(out2) = names(out)
  
  return(out2)
}


# x un individu de type msm et t un temps
id_get_state <- function(x, t) 
{
  aux <- max(which(x[,"time"]  <= t))
  return(x[aux,"state"])
}



#' Estimate probabilities to be in each state
#' 
#' @param data_msm data.frame containing \code{id}, \code{time} and \code{state} (see \code{\link{generate_Markov_cfd}}). All individual must end at the same Tmax (use \code{\link{msm2msmTmax}}).
#' 
#' @return A list of two elements:
#' \itemize{
#'   \item{t: vector of time}
#'   \item{pt: a matrix with K (= number of states) rows and with \code{length(t)} columns containing the probabilities to be in each state at each time.}
#' }
#' 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' d_JK2 <- msm2msmTmax(d_JK, 10)
#' 
#' estimate_pt(d_JK2)
#' 
#' @author Cristian Preda
#' 
#' @export   
estimate_pt <- function(data_msm)
{
  t_jumps <- sort(unique(data_msm$time)) 
  n <- length(unique(data_msm$id))
  states <- unique(data_msm[,"state"])  
  res <- matrix(0, nrow = length(states), ncol = length(t_jumps), dimnames = list(1:length(states), round(t_jumps, 3)))
  
  for(i in seq_along(t_jumps))
  {
    aux <- as.vector(by(data_msm, data_msm$id, function(x){id_get_state(x, t_jumps[i])})) 
    # donne pour chaque temps, le nb d'individus dans chacun des état
    res[,i] = tabulate(aux, nbins = length(states))/n
  }
  return(list(pt = res, t = t_jumps))
}


#' Plot probabilities
#'
#' Plot the probabilities of each state at each given time
#'
#' @param pt output of \code{\link{estimate_pt}}
#' @param ribbon if TRUE, use ribbon to plot probabilities 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' d_JK2 <- msm2msmTmax(d_JK, 10)
#' 
#' pt <- estimate_pt(d_JK2)
#' 
#' plot_pt(pt, ribbon = TRUE)
#' 
#' @author Quentin Grimonprez
#' 
#' @export
plot_pt <- function(pt, ribbon = FALSE)
{
  if(ribbon)
    p <- plot_pt_ribbon(pt)
  else
    p <- plot_pt_classic(pt)
  
  return(p)
}


# plot line
# @author Quentin Grimonprez
plot_pt_classic <- function(pt)
{
  plot_data <- data.frame(state = as.factor(rep(1:nrow(pt$pt), each = ncol(pt$pt))), 
                          proba = as.vector(t(pt$pt)), 
                          time = rep(pt$t, nrow(pt$pt)))
  
  p <- ggplot(plot_data, aes_string(x = "time", y = "proba", group = "state", colour = "state")) +
    geom_point() + geom_line() + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)")
  
  return(p)
}



# plot probabilities using ribbon
# @author Quentin Grimonprez
plot_pt_ribbon <- function(pt)
{
  plot_data <- as.data.frame(t(apply(pt$pt, 2, cumsum)))
  nState <- ncol(plot_data)
  names(plot_data) = paste0("state", names(plot_data))
  plot_data$time = pt$t
  plot_data$state0 = rep(0, nrow(plot_data))
  
  p <- ggplot(plot_data)
  for(i in 1:nState)
    p = p + geom_ribbon(aes_string(ymin = paste0("state", i-1), 
                                   ymax = paste0("state", i), x = "time", fill = factor(i, levels = 1:nState)), colour = "black", 
                         alpha = 0.8)
  
  p = p  + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)", fill = "state")
  
  return(p)
}

