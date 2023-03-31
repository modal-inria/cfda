#' Compute time spent in each state
#'
#' For each individual, compute the time spent in each state
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#'
#' @return a matrix with \code{K} columns containing the total time spent in each state for each individual
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # cut at Tmax = 8
#' d_JK2 <- cut_data(d_JK, Tmax = 8)
#'
#' # compute time spent by each id in each state
#' timeSpent <- compute_time_spent(d_JK2)
#' @seealso \link{boxplot.timeSpent}
#' @family Descriptive statistics
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_time_spent <- function(data) {
  ## check parameters
  checkData(data)
  ## end check

  if (is.factor(data$state)) {
    labels <- levels(data$state)
  } else {
    labels <- sort(unique(data$state))
  }

  res <- by(data, data$id, function(x) {
    compute_time_spent_intern(x, labels)
  })
  out <- do.call(rbind, res)
  colnames(out) <- labels
  class(out) <- "timeSpent"

  return(out)
}


# How long an individual stays in each state
# @author Cristian Preda
compute_time_spent_intern <- function(data, labels) {
  aux <- rep(0, length(labels))
  for (i in seq_along(labels)) {
    idx <- which(data$state == labels[i])
    for (u in idx) {
      if (u < nrow(data)) {
        aux[i] <- aux[i] + data$time[u + 1] - data$time[u]
      }
    }
  }
  return(aux)
}

#' Boxplot of time spent in each state
#'
#' @param x output of \code{\link{compute_time_spent}} function
#' @param col a vector containing color for each state
#' @param ... extra parameters for \code{geom_boxplot}
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # cut at Tmax = 8
#' d_JK2 <- cut_data(d_JK, Tmax = 8)
#'
#' # compute time spent by each id in each state
#' timeSpent <- compute_time_spent(d_JK2)
#'
#' # plot the result
#' boxplot(timeSpent, col = c("#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"))
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' boxplot(timeSpent, notch = TRUE, outlier.colour = "black") +
#'   coord_flip() +
#'   labs(title = "Time spent in each state")
#' @author Quentin Grimonprez
#' @seealso \link{compute_time_spent}
#' @family Descriptive statistics
#'
#' @export
boxplot.timeSpent <- function(x, col = NULL, ...) {
  df <- data.frame(timeSpent = as.vector(x), state = factor(rep(colnames(x), each = nrow(x)), levels = colnames(x)))
  p <- ggplot(df, aes_string(x = "state", y = "timeSpent", fill = "state")) +
    geom_boxplot(...) +
    labs(x = "State", y = "Time Spent", fill = "State")

  if (!is.null(col)) {
    p <- p + scale_fill_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData

  return(p)
}


#' Compute duration of individuals
#'
#' For each individual, compute the duration
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#'
#' @return a vector containing the duration of each trajectories
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#'
#' # compute duration of each individual
#' duration <- compute_duration(d_JK)
#'
#' hist(duration)
#' @seealso \link{hist.duration}
#' @family Descriptive statistics
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_duration <- function(data) {
  ## check parameters
  checkData(data)
  ## end check

  out <- tapply(data$time, as.factor(data$id), function(x) diff(range(x)))
  class(out) <- "duration"

  return(out)
}

#' Plot the duration
#'
#'
#' @param x output of \code{\link{compute_duration}} function
#' @param breaks number of breaks. If not given, use the Sturges rule
#' @param ... parameters for \code{geom_histogram}
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#'
#' # compute duration of each individual
#' duration <- compute_duration(d_JK)
#'
#' hist(duration)
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' hist(duration) +
#'   labs(title = "Distribution of the duration")
#' @author Quentin Grimonprez
#' @seealso \link{compute_duration}
#' @family Descriptive statistics
#'
#' @export
hist.duration <- function(x, breaks = NULL, ...) {
  # choose the number of breaks using Sturges rule
  if (is.null(breaks)) {
    breaks <- floor(1 + log2(length(x)))
  }

  extraParam <- list(...)
  defaultParam <- list(fill = "lightblue", color = "black", bins = breaks)
  param <- c(extraParam, defaultParam[which(!(names(defaultParam) %in% names(extraParam)))])

  ggplot(data.frame(duration = as.vector(x)), aes_string(x = "duration")) +
    do.call(geom_histogram, param) +
    labs(x = "Duration", y = "Frequency")
}


#' Extract the state of each individual at a given time
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
#' @param t time at which extract the state
#' @param NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
#' (useful when individuals has different lengths)
#'
#' @return a vector containing the state of each individual at time t
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # get the state of each individual at time t = 6
#' get_state(d_JK, 6)
#'
#'
#' # get the state of each individual at time t = 12 (> Tmax)
#' get_state(d_JK, 12)
#' # if NAafterTmax = TRUE, it will return NA for t > Tmax
#' get_state(d_JK, 12, NAafterTmax = TRUE)
#'
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
get_state <- function(data, t, NAafterTmax = FALSE) {
  ## check parameters
  checkData(data)
  if (any(is.na(t)) || !is.numeric(t) || (length(t) != 1)) {
    stop("t must be a real.")
  }
  ## end check

  out <- by(data, data$id, function(x) {
    id_get_state(x, t, NAafterTmax)
  })
  out2 <- as.vector(out)
  names(out2) <- names(out)

  return(out2)
}

# return the state at time t
#
# x cfda dataframe
# t time value
# NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
# @author Cristian Preda, Quentin Grimonprez
id_get_state <- function(x, t, NAafterTmax = FALSE) {
  if (NAafterTmax && (t > x$time[length(x$time)])) {
    return(NA)
  }

  aux <- which(x$time <= t)
  return(x$state[aux[length(aux)]])
}



#' Estimate probabilities to be in each state
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#' @param NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
#' (useful when individuals has different lengths)
#'
#' @return A list of two elements:
#' \itemize{
#'   \item{t: vector of time}
#'   \item{pt: a matrix with K (= number of states) rows and with \code{length(t)} columns containing the
#' probabilities to be in each state at each time.}
#' }
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' d_JK2 <- cut_data(d_JK, 10)
#'
#' # estimate probabilities
#' estimate_pt(d_JK2)
#' @author Cristian Preda, Quentin Grimonprez
#' @seealso \link{plot.pt}
#' @family Descriptive statistics
#'
#' @export
estimate_pt <- function(data, NAafterTmax = FALSE) {
  ## check parameters
  checkData(data)
  ## end check

  t_jumps <- sort(unique(data$time))
  uniqueId <- unique(data$id)

  if (is.factor(data$state)) {
    states <- levels(data$state)
  } else {
    states <- sort(unique(data$state))
  }

  res <- matrix(0, nrow = length(states), ncol = length(t_jumps), dimnames = list(states, round(t_jumps, 3)))

  for (id in uniqueId) {
    x <- data[data$id == id, ]
    for (i in seq_along(t_jumps)) {
      aux <- id_get_state(x, t_jumps[i], NAafterTmax)

      res[match(aux, states), i] <- res[match(aux, states), i] + 1
    }
  }

  res <- prop.table(res, margin = 2)

  out <- list(pt = res, t = t_jumps)
  class(out) <- "pt"

  return(out)
}

# Extract probability to be in each state at a given time
#
# @param pt output of \link{estimate_pt} function
# @param t time value at which the probability is required
#
# @return  probability to be in each state at time t
#
#
# @examples
# # Simulate the Jukes-Cantor model of nucleotide replacement
# K <- 4
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#
# d_JK2 <- cut_data(d_JK, 10)
#
# # estimate probabilities
# pt <- estimate_pt(d_JK2)
#
# get_proba(pt, 1.5)
#
#
# @seealso \link{estimate_pt}
#
# @author Quentin Grimonprez
#
# @export
get_proba <- function(pt, t) {
  # if we do not export this function, there is no need to do theses checks
  # if(!inherits(x, "pt"))
  #   stop("pt must be an object of class pt.")
  # if(any(is.na(t)) || !is.numeric(t) || length(t) != 1)
  #   stop("t must be a real.")

  # find the index containing the first time value greater than the given time t
  i <- sum(t >= pt$t)

  # if i == 0, the given time is lower than any time in pt, we can't estimate probabilities, we return NA
  if (i == 0) {
    p <- rep(NA, nrow(pt$pt))
    names(p) <- rownames(pt$pt)

    return(p)
  }

  return(pt$pt[, i])
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
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' d_JK2 <- cut_data(d_JK, 10)
#'
#' pt <- estimate_pt(d_JK2)
#'
#' plot(pt, ribbon = TRUE)
#' @author Quentin Grimonprez
#' @method plot pt
#' @seealso \link{estimate_pt}
#' @family Descriptive statistics
#'
#' @export
plot.pt <- function(x, col = NULL, ribbon = FALSE, ...) {
  ## check parameters
  checkLogical(ribbon, "ribbon")
  ## end check

  if (ribbon) {
    p <- plot_pt_ribbon(x, col, ...)
  } else {
    p <- plot_pt_classic(x, col)
  }

  return(p)
}


# plot line
# @author Quentin Grimonprez
plot_pt_classic <- function(pt, col = NULL) {
  plot_data <- data.frame(
    State = factor(rep(rownames(pt$pt), each = ncol(pt$pt)), levels = rownames(pt$pt)),
    proba = as.vector(t(pt$pt)),
    time = rep(pt$t, nrow(pt$pt))
  )

  p <- ggplot(plot_data, aes_string(x = "time", y = "proba", group = "State", colour = "State")) +
    geom_line() +
    ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)")

  if (!is.null(col)) {
    p <- p + scale_colour_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData

  return(p)
}



# plot probabilities using ribbon
# @author Quentin Grimonprez
plot_pt_ribbon <- function(pt, col = NULL, addBorder = TRUE) {
  ## check parameters
  checkLogical(addBorder, "addBorder")
  ## end check

  plot_data <- as.data.frame(t(apply(pt$pt, 2, cumsum)))
  nState <- ncol(plot_data)
  labels <- paste0("state", names(plot_data))
  shortLabels <- factor(names(plot_data), levels = names(plot_data))
  names(plot_data) <- labels
  plot_data$time <- pt$t
  plot_data$state0 <- rep(0, nrow(plot_data))
  labels <- c("state0", labels)

  p <- ggplot(plot_data)
  for (i in seq_len(nState)) {
    p <- p + geom_ribbon(aes_string(
      ymin = paste0("`", labels[i], "`"),
      ymax = paste0("`", labels[i + 1], "`"), x = "time",
      fill = shortLabels[i]
    ),
    colour = ifelse(addBorder, "black", NA), alpha = 0.8
    )
  }

  if (!is.null(col)) {
    p <- p + scale_fill_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData

  p <- p + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)", fill = "State")

  return(p)
}


#' Compute the number of jumps
#'
#' For each individual, compute the number of jumps performed
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
#' @param countDuplicated if \code{TRUE}, jumps in the same state are counted as jump
#'
#' @return A vector containing the number of jumps for each individual
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # compute the number of jumps
#' nJump <- compute_number_jumps(d_JK)
#' @seealso \link{hist.njump}
#' @family Descriptive statistics
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_number_jumps <- function(data, countDuplicated = FALSE) {
  ## check parameters
  checkData(data)
  checkLogical(countDuplicated, "countDuplicated")
  ## end check

  out <- by(data, data$id, function(x) {
    compute_number_jumpsIntern(x, countDuplicated)
  })
  nom <- names(out)
  out <- as.vector(out)
  names(out) <- nom
  class(out) <- "njump"

  return(out)
}

# @param state vector with state, ordered by time
# @param countDuplicated if TRUE jump in the same state are counted
# @author Quentin Grimonprez
compute_number_jumpsIntern <- function(x, countDuplicated = TRUE) {
  if (countDuplicated) {
    return(length(x$state) - 1)
  } else {
    out <- rle(as.character(x$state[order(x$time)]))$values # rle does not manage factor, as.character allows it
    return(length(out) - 1)
  }
}


#' Plot the number of jumps
#'
#'
#' @param x output of \code{\link{compute_number_jumps}} function
#' @param breaks number of breaks. If not given, use the Sturges rule
#' @param ... parameters for \code{geom_histogram}
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' nJump <- compute_number_jumps(d_JK)
#'
#' hist(nJump)
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' hist(nJump, fill = "#984EA3") +
#'   labs(title = "Distribution of the number of jumps")
#' @author Quentin Grimonprez
#' @seealso \link{compute_number_jumps}
#' @family Descriptive statistics
#'
#' @export
hist.njump <- function(x, breaks = NULL, ...) {
  # choose the number of breaks using Sturges rule
  if (is.null(breaks)) {
    breaks <- min(floor(1 + log2(length(x))), max(x) + 1)
  }

  extraParam <- list(...)
  defaultParam <- list(fill = "lightblue", color = "black", bins = breaks, center = 0)
  param <- c(extraParam, defaultParam[which(!(names(defaultParam) %in% names(extraParam)))])

  ggplot(data.frame(njump = as.vector(x)), aes_string(x = "njump")) +
    do.call(geom_histogram, param) +
    labs(x = "Number of jumps", y = "Frequency") +
    scale_x_continuous(breaks = function(x) pretty(seq(ceiling(x[1]), floor(x[2]), by = 1)))
}


#' Table of transitions
#'
#' Calculates a frequency table counting the number of times each pair of states were observed in successive observation times.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#' @param removeDiagonal if TRUE, does not count transition from a state i to i
#'
#' @return a matrix of size \code{K*K} containing the number of transition for each pair
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # table of transitions
#' statetable(d_JK)
#' @author Quentin Grimonprez
#' @family Descriptive statistics
#'
#' @export
statetable <- function(data, removeDiagonal = FALSE) {
  ## check parameters
  checkData(data)
  ## end check

  newState <- stateToInteger(data$state)

  out <- statetable.msm(newState$state, data$id)

  # If there is at least 1 absorbing state, the matrix is not a square matrix
  out <- completeStatetable(out)

  colnames(out) <- newState$label$label[match(colnames(out), newState$label$code)]
  rownames(out) <- newState$label$label[match(rownames(out), newState$label$code)]

  if (removeDiagonal) {
    diag(out) <- 0
  }

  return(out)
}
