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
#' timeSpent <- compute_time_spent(d_JK2)
#' 
#' @seealso \link{boxplot.timeSpent}
#' 
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_time_spent <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  labels <- sort(unique(data_msm$state))
  
  res <- by(data_msm, data_msm$id, function(x){compute_time_spent_intern(x, labels)})
  out <- do.call(rbind, res)
  colnames(out) = labels
  class(out) = "timeSpent"
  
  return(out)
}


# combien de temps passe un id dans [0,T] dans chaque etat, x vient d'un msmT
compute_time_spent_intern <- function(data_msm, labels)
{
  aux <- rep(0, length(labels))            
  for(i in seq_along(labels))
  {
    idx <- which(data_msm$state == labels[i])
    for(u in idx)
    {
      if(u < nrow(data_msm))
        aux[i] = aux[i] + data_msm$time[u+1] - data_msm$time[u]
    }  
  }
  return(aux)
}          

#' Boxplot of time spent in each state
#'
#' @param x output of \code{\link{compute_time_spent}} function
#' @param col a vector containing color for each state
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
#' timeSpent <- compute_time_spent(d_JK2)
#' 
#' # plot the result
#' boxplot(timeSpent, col = c("#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"))
#' 
#' @author Quentin Grimonprez
#'
#' @export
boxplot.timeSpent <- function(x, col = NULL, ...)
{
  p <- ggplot(data.frame(timeSpent = as.vector(x), state = rep(colnames(x), each = nrow(x))), 
              aes_string(x = "state", y = "timeSpent", fill = "state")) + 
    geom_boxplot() + labs(x = "State", y = "Time Spent", fill = "State")
  
  if(!is.null(col))
    p = p + scale_fill_manual(values = col)
  
  return(p)
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
#' duration <- compute_duration(d_JK)
#' 
#' hist(duration)
#' 
#' @seealso \link{hist.duration}
#' 
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_duration <- function(data_msm)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  out <- tapply(data_msm$time, as.factor(data_msm$id), function(x) diff(range(x)))
  class(out) = "duration"
  
  return(out)
}

#' Plot the duration
#' 
#' 
#' @param x output of \code{\link{compute_duration}} function
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
#'
#' # compute duration of each individual
#' duration <- compute_duration(d_JK)
#' 
#' hist(duration)
#' 
#' @author Quentin Grimonprez
#' 
#' @export 
hist.duration <- function(x, breaks = NULL, ...)
{
  # choose the number of breaks using sturges rule
  if(is.null(breaks))
    breaks <- floor(1 + log2(length(x)))
  
  ggplot(data.frame(duration = as.vector(x)), aes_string(x = "duration"))+
    geom_histogram(fill = "lightblue", color = "black", bins = breaks) +
    labs(x = "Duration", y = "Frequency")
}


#' Extract the state of each individual at a given time
#' 
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#' @param t time at which extract the state
#' @param NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax (usefull when individuals has different lengths)
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
#' 
#' # get the state of each individuals at time t = 12 (> Tmax)
#' get_state(d_JK, 12)
#' # if NAafterTmax = TRUE, it will return NA for t > Tmax
#' get_state(d_JK, 12, NAafterTmax = TRUE)
#' 
#' @author Cristian Preda, Quentin Grimonprez
#' 
#' @export
get_state <- function(data_msm, t, NAafterTmax = FALSE) 
{
  ## check parameters
  checkDataMsm(data_msm)
  if(!is.numeric(t) || (length(t) != 1))
    stop("t must be a real.")
  ## end check
  
  out <- by(data_msm, data_msm$id, function(x){id_get_state(x, t, NAafterTmax)})
  out2 <- as.vector(out)
  names(out2) = names(out)
  
  return(out2)
}

# return the state at time t
#
# x un individu de type msm et t un temps
# NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
id_get_state <- function(x, t, NAafterTmax = FALSE) 
{
  if(NAafterTmax && (t > x$time[length(x$time)]))
    return(NA)
  
  aux <- which(x$time <= t)
  return(x$state[aux[length(aux)]])
}



#' Estimate probabilities to be in each state
#' 
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state. All individual must end at the same time Tmax (use \code{\link{msm2msmTmax}}).
#' @param NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax (usefull when individuals has different lengths)
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
#' @author Cristian Preda, Quentin Grimonprez
#' 
#' @export   
estimate_pt <- function(data_msm, NAafterTmax = FALSE)
{
  ## check parameters
  checkDataMsm(data_msm)
  ## end check
  
  t_jumps <- sort(unique(data_msm$time)) 
  uniqueId <- unique(data_msm$id)
  n <- length(uniqueId)
  states <- sort(unique(data_msm$state))
  res <- matrix(0, nrow = length(states), ncol = length(t_jumps), dimnames = list(states, round(t_jumps, 3)))
  
  for(id in uniqueId)
  {
    x <- data_msm[data_msm$id == id,]
    for(i in seq_along(t_jumps))
    {
      aux <- id_get_state(x, t_jumps[i], NAafterTmax) 
      
      res[match(aux, states), i] = res[match(aux, states), i] + 1
    }
  }
  
  res = prop.table(res, margin = 2)
  
  out <- list(pt = res, t = t_jumps)
  class(out) = "pt"
  
  return(out)
}


#' Plot probabilities
#'
#' Plot the probabilities of each state at each given time
#'
#' @param x output of \code{\link{estimate_pt}}
#' @param col a vector containing color for each state
#' @param ribbon if TRUE, use ribbon to plot probabilities 
#' @param ... only if \code{ribbon = TRUE}, parameter \code{addBorder}, if TRUE, add black border to the ribbons.
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
plot.pt <- function(x, col = NULL, ribbon = FALSE, ...)
{
  ## check parameters
  checkLogical(ribbon, "ribbon")
  ## end check
  
  if(ribbon)
    p <- plot_pt_ribbon(x, col, ...)
  else
    p <- plot_pt_classic(x, col)
  
  return(p)
}


# plot line
# @author Quentin Grimonprez
plot_pt_classic <- function(pt, col = NULL)
{
  plot_data <- data.frame(state = as.factor(rep(rownames(pt$pt), each = ncol(pt$pt))), 
                          proba = as.vector(t(pt$pt)), 
                          time = rep(pt$t, nrow(pt$pt)))
  
  p <- ggplot(plot_data, aes_string(x = "time", y = "proba", group = "state", colour = "state")) +
    geom_line() + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)")
  
  if(!is.null(col))
    p = p + scale_colour_manual(values = col)
  
  return(p)
}



# plot probabilities using ribbon
# @author Quentin Grimonprez
plot_pt_ribbon <- function(pt, col = NULL, addBorder = TRUE)
{
  ## check parameters
  checkLogical(addBorder, "addBorder")
  ## end check
  
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
    p = p + geom_ribbon(aes_string(ymin = paste0("`", labels[i], "`"), 
                                   ymax = paste0("`", labels[i+1], "`"), x = "time", 
                                   fill = factor(shortLabels[i], levels = shortLabels)), 
                        colour = ifelse(addBorder, "black", NA), alpha = 0.8)
  
  if(!is.null(col))
    p = p + scale_fill_manual(values = col)
  
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
    out <- rle(as.character(x$state[order(x$time)]))$values #rle does not manage factor, as.character allows it
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
  
  newState = stateToInteger(data_msm$state)
  
  out <- statetable.msm(newState$state, data_msm$id)
  colnames(out) = newState$label$label[match(colnames(out), newState$label$code)]
  rownames(out) = newState$label$label[match(rownames(out), newState$label$code)]
  
  return(out)
}


#' Plot categorical functional data
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#' @param col a vector containing color for each state
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
plotData <- function(data_msm, col = NULL, addLabel = TRUE, addBorder = TRUE)
{
  ## check parameters
  checkDataMsm(data_msm)
  checkLogical(addLabel, "addLabel")
  checkLogical(addBorder, "addBorder")
  ## end check
  
  d_graph <- rep_large_ind(data_msm)

  d_graph$state = factor(d_graph$state, levels = sort(unique(data_msm$state)))
  
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
    
  if(!is.null(col))
    p = p + scale_fill_manual(values = col)
  
  return(p)
}


# transform the data_msm format to a new format with 4 columns: id, t_stast, t_end, state.
# usefull for ggplot
rep_large_ind <- function(data_msm)
{
  out <- by(data_msm, data_msm$id, function(x){
    data.frame(id = x$id[1:(nrow(x)-1)] , 
               t_start = x$time[1:(nrow(x)-1)], 
               t_end = x$time[2:nrow(x)], 
               state = x$state[1:(nrow(x)-1)])
  })
  
  return(do.call(rbind, out))
}

