#' Compute time spent in each state
#'
#' For each individual, compute the time spent in each state
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#'
#' @return a matrix with \code{K} columns containing the total time spent in each state for each individuals
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
#' # compute time spent by each id in each state
#' timeSpent <- compute_Time_Spent(d_JK2)
#' 
#' @seealso \link{boxplot.timeSpent}
#' 
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_Time_Spent <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  labels <- sort(unique(data_msm$state))
  
  res <- by(data_msm, data_msm$id, function(x){compute_Time_Spent_intern(x, labels)})
  out <- do.call(rbind, res)
  colnames(out) = labels
  class(out) = "timeSpent"
  
  return(out)
}


# combien de temps passe un id dans [0,T] dans chaque etat, x vient d'un msmT
compute_Time_Spent_intern <- function(data_msm, labels)
{
  aux <- rep(0, length(labels))            
  for(i in seq_along(labels))
  {
    idx <- which(data_msm$state == labels[i])
    for(u in idx)
    {
      if(u < nrow(data_msm))
        aux[i] = aux[i] + data_msm[u+1, "time"] - data_msm[u, "time"]
    }  
  }
  return(aux)
}          

#' Boxplot of time spent in each state
#'
#' @param x output of \code{\link{compute_Time_Spent}} function
#' @param ... not used
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
#' # compute time spent by each id in each state
#' timeSpent <- compute_Time_Spent(d_JK2)
#' 
#' # plot the result
#' boxplot(timeSpent)
#' 
#' @author Quentin Grimonprez
#'
#' @export
boxplot.timeSpent <- function(x, ...)
{
  ggplot(data.frame(timeSpent = as.vector(x), state = rep(colnames(x), each = nrow(x))), 
         aes_string(x = "state", y = "timeSpent", fill = "state")) + 
    geom_boxplot() + labs(x = "State", y = "Time Spent", fill = "State")
}


#' Compute duration of individuals
#'
#' For each individual, compute the duration
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#'
#' @return a vector containing the duration of each trajectories
#'
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#'
#' # compute duration of each individual
#' duration <- compute_Duration(d_JK)
#' 
#' 
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_Duration <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  out <- tapply(data_msm$time, as.factor(data_msm$id), function(x) diff(range(x)))

  return(out)
}

#' Extract the state of each individual at a given time
#' 
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
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
#' # get the state of each individuals at time t = 6
#' get_state(d_JK, 6)
#' 
#' @author Cristian Preda
#' 
#' @export
get_state <- function(data_msm, t) 
{
  ## check parameters
  checkDataMsm(data_msm)
  if(!is.numeric(t) || (length(t) != 1))
    stop("t must be a real.")
  ## end check
  
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
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state. All individual must end at the same time Tmax (use \code{\link{msm2msmTmax}}).
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
#' # estimate probabilities
#' estimate_pt(d_JK2)
#' 
#' @seealso \link{plot.pt}
#' 
#' @author Cristian Preda
#' 
#' @export   
estimate_pt <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  t_jumps <- sort(unique(data_msm$time)) 
  n <- length(unique(data_msm$id))
  states <- sort(unique(data_msm$state))
  res <- matrix(0, nrow = length(states), ncol = length(t_jumps), dimnames = list(states, round(t_jumps, 3)))
  
  for(i in seq_along(t_jumps))
  {
    aux <- as.vector(by(data_msm, data_msm$id, function(x){id_get_state(x, t_jumps[i])})) 
    
    # donne pour chaque temps, le nb d'individus dans chacun des état
    tab <- table(aux)/n
    res[match(names(tab), states), i] = tab
  }
  
  out <- list(pt = res, t = t_jumps)
  class(out) = "pt"
  
  return(out)
}


#' Plot probabilities
#'
#' Plot the probabilities of each state at each given time
#'
#' @param x output of \code{\link{estimate_pt}}
#' @param ribbon if TRUE, use ribbon to plot probabilities 
#' @param ... unused
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
#' plot(pt, ribbon = TRUE)
#' 
#' @author Quentin Grimonprez
#' 
#' @method plot pt
#' 
#' @export
plot.pt <- function(x, ribbon = FALSE, ...)
{
  ## check parameters
  checkLogical(ribbon, "ribbon")
  ## end check
  
  if(ribbon)
    p <- plot_pt_ribbon(x)
  else
    p <- plot_pt_classic(x)
  
  return(p)
}


# plot line
# @author Quentin Grimonprez
plot_pt_classic <- function(pt)
{
  plot_data <- data.frame(state = as.factor(rep(rownames(pt$pt), each = ncol(pt$pt))), 
                          proba = as.vector(t(pt$pt)), 
                          time = rep(pt$t, nrow(pt$pt)))
  
  p <- ggplot(plot_data, aes_string(x = "time", y = "proba", group = "state", colour = "state")) +
    geom_line() + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)")
  
  return(p)
}



# plot probabilities using ribbon
# @author Quentin Grimonprez
plot_pt_ribbon <- function(pt)
{
  plot_data <- as.data.frame(t(apply(pt$pt, 2, cumsum)))
  nState <- ncol(plot_data)
  labels <- paste0("state", names(plot_data))
  shortLabels <- names(plot_data)
  names(plot_data) = labels
  plot_data$time = pt$t
  plot_data$state0 = rep(0, nrow(plot_data))
  labels = c("state0", labels)
  
  p <- ggplot(plot_data)
  for(i in 1:nState)
    p = p + geom_ribbon(aes_string(ymin = labels[i], 
                                   ymax = labels[i+1], x = "time", fill = factor(shortLabels[i], levels = shortLabels)), colour = "black", 
                         alpha = 0.8)
  
  p = p  + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)", fill = "state")
  
  return(p)
}


#' Compute the number of jumps
#' 
#' For each individual, compute the number of jumps performed
#' 
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#' @param countDuplicated if \code{TRUE}, jumps in the same state are counted as jump
#' 
#' @return A vector containing the number of jumps for each individual 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' # compute the number of jumps
#' nJump <- compute_number_jumps(d_JK)
#' 
#' @seealso \link{hist.njump}
#' 
#' @author Cristian Preda, Quentin Grimonprez
#' 
#' @export  
compute_number_jumps <- function(data_msm, countDuplicated = TRUE)
{
  ## check parameters
  checkDataMsm(data_msm)
  checkLogical(countDuplicated, "countDuplicated")
  ## end check
  
  out <- by(data_msm, data_msm$id, function(x){compute_number_jumpsIntern(x, countDuplicated)})
  nom <- names(out)
  out <- as.vector(out)
  names(out) = nom
  class(out) = "njump"
  
  return(out)
}

# @param state vector with state, ordered by time
# @param countDuplicated if TRUE jump in the same state are counted
# @author Quentin Grimonprez
compute_number_jumpsIntern <- function(x, countDuplicated = TRUE)
{
  if(countDuplicated)
  {
    return(length(x$state)-1)
  }else{
    out <- rle(x$state[order(x$time)])$values
    return(length(out)-1)
  }
}


#' Plot the number of jumps
#' 
#' 
#' @param x output of \code{\link{compute_number_jumps}} function
#' @param breaks number of breaks. If not given use the sturges rule
#' @param ... not used
#' 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' nJump <- compute_number_jumps(d_JK)
#' 
#' hist(nJump)
#' 
#' @author Quentin Grimonprez
#' 
#' @export  
hist.njump <- function(x, breaks = NULL, ...)
{
  # choose the number of breaks using sturges rule
  if(is.null(breaks))
    breaks <- floor(1 + log2(length(x)))
  
  ggplot(data.frame(njump = as.vector(x)), aes_string(x = "njump"))+
    geom_histogram(fill = "lightblue", color = "black", bins = breaks) +
    labs(x = "Number of jumps", y = "Frequency")
}


#' Table of transitions
#'
#' Calculates a frequency table counting the number of times each pair of states were observed in successive observation times.
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
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
#' # table of transitions
#' statetable(d_JK)
#' 
#'
#' @export
statetable <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  statetable.msm(data_msm$state, data_msm$id)
}


#' Plot categorical functional data
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#' @param addLabel If TRUE, add id labels
#' @param addBorder If TRUE add black border to each individuals
#' 
#' @return
#' each row represent an individual over [0:Tmax]. The color at a given time gives the state of the individual. 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = 10)
#' 
#' # add a line with time Tmax at the end of each individuals
#' d_JKT <- msm2msmTmax(d_JK, Tmax = 10)
#' 
#' plotData(d_JKT)
#'   
#' @author Cristian Preda, Quentin Grimonprez
#' 
#' @export
plotData <- function(data_msm, addLabel = TRUE, addBorder = TRUE)
{
  ## check parameters
  checkDataMsm(data_msm)
  checkLogical(addLabel, "addLabel")
  checkLogical(addBorder, "addBorder")
  ## end check
  
  d_graph <- rep_large_ind(data_msm)

  d_graph$state = factor(d_graph$state, levels = unique(data_msm$state))
  
  nInd <- length(unique(d_graph$id))
  d_graph$id2 <- unclass(factor(d_graph$id))
  
  p <- ggplot() + 
    scale_x_continuous(name = "time") + 
    geom_rect(data = d_graph, mapping = aes_string(xmin = "t_start", xmax = "t_end", ymin = "id2 - 0.5", ymax = "id2 + 0.5", fill = "state"), 
              color = ifelse(addBorder, "black", NA), alpha = 0.7)
  
  if(addLabel)
  {
    p = p + scale_y_continuous(name = "id", breaks = 1:nInd)
  }else{
    p = p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }
    
  p
}


# transform the data_msm format to a new format with 4 columns: id, t_stast, t_end, state.
# usefull for ggplot
rep_large_ind <- function(data_msm)
{
  out <- by(data_msm, data_msm$id, function(x){
    data.frame(id = x[1:(nrow(x)-1), "id"] , 
               t_start = x[1:(nrow(x)-1), "time"], 
               t_end = x[2:nrow(x), "time"], 
               state = x[1:(nrow(x)-1), "state"])
  })
  
  return(do.call(rbind, out))
}

