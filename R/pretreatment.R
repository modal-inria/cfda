#' Cut data to a maximal given time
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time},
#' time at which a change occurs and \code{state}, associated state.
#' @param Tmax max time considered
#' @param prolongLastState list of states to prolong (can be "all"). In the case where the last state of a trajectory is
#' lesser than \code{Tmax}, we can assume that this trajectory will be in the same state at time Tmax only if it is an
#' absorbing state. Otherwise it will add \code{NAstate} and throw a warning.
#' Set `prolongLastState = c()` to indicate there is no absorbing state.
#' @param NAstate state value used when the last state is not prolonged.
#' @param warning if TRUE, the function raises warnings when it has prolonged a trajectory with NAstate
#'
#' @return a data.frame with the same format as \code{data} where each individual has \code{Tmax} as last time entry.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' set.seed(42)
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#' tail(d_JK)
#'
#' # cut at Tmax = 8
#' d_JK2 <- cut_data(d_JK, Tmax = 8)
#' tail(d_JK2)
#'
#' try(d_JK2 <- cut_data(d_JK, Tmax = 12, prolongLastState = c()))
#'
#' @author Cristian Preda
#'
#' @export
cut_data <- function(data, Tmax, prolongLastState = "all", NAstate = "Not observable", warning = FALSE) {
  ## check parameters
  checkData(data)
  checkLogical(warning, "warning")
  if (length(NAstate) > 1) {
    stop("NAstate must have a length of 1")
  }
  if (any(is.na(Tmax)) || !is.numeric(Tmax) || (length(Tmax) != 1)) {
    stop("Tmax must be a real.")
  }
  ## end check

  d <- do.call(rbind, by(data, data$id, function(x) {
    cut_cfd(x, Tmax, prolongLastState, NAstate, warning)
  }))
  rownames(d) <- NULL

  return(d)
}

# @author Cristian Preda
cut_cfd <- function(data, Tmax, prolongLastState = "all", NAstate = NA, warning = FALSE) {
  l <- nrow(data)
  currTmax <- max(data$time)

  if (Tmax > currTmax) {
    if (((length(prolongLastState) > 0) && all(prolongLastState == "all")) || (data$state[l] %in% prolongLastState)) {
      return(rbind(data, data.frame(id = data$id[1], state = data$state[l], time = Tmax)))
    } else {
      if (warning) {
        warning(paste0("id ", data$id[1], " does not end with an absorbing state. Cannot impute the state until time ",
                       Tmax, ". Please, add more records or change the Tmax value."))
      }
      d <- data
      d$state[l] <- NAstate
      return(rbind(d, data.frame(id = data$id[1], state = NAstate, time = Tmax)))
    }
  } else {
    if (currTmax == Tmax) {
      return(data)
    } else {
      if (Tmax %in% data$time) {
        k <- which(data$time == Tmax)
        return(data[1:k, ])
      } else {
        k <- max(which(data$time <= Tmax))
        return(rbind(data[1:k, ], data.frame(state = data$state[k], time = Tmax, id = data$id[1])))
      }
    }
  }
}

# change the labels into integer
#
# @param state vector with labels
# @return a list with state containing the new formatted state and label,
# a data.frame containing the labels and corresponding integers
#
# @author Quentin Grimonprez
stateToInteger <- function(state) {
  lab <- data.frame(label = sort(unique(state)), code = seq_along(unique(state)))

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
refactorCategorical <- function(data, oldCateg = unique(data), newCateg = seq_along(oldCateg)) {
  ind <- match(data, oldCateg)

  if (any(is.na(ind[!is.na(data)]))) {
    warning("NA produced.")
  }

  return(newCateg[ind])
}


#' Remove duplicated states
#'
#' Remove duplicated consecutive states from data.
#' If for an individual there is two or more consecutive states that are identical, only the first is kept.
#' Only time when the state changes are kept.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#' @param keep.last if TRUE, keep the last state for every individual even if it is a duplicated state.
#'
#' @return \code{data} without duplicated consecutive states
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:3, c(10, 3, 8)), time = c(1:10, 1:3, 1:8),
#'   state = c(rep(1:5, each = 2), 1:3, rep(1:3, c(1, 6, 1)))
#' )
#' out <- remove_duplicated_states(data)
#' @author Quentin Grimonprez
#'
#' @export
remove_duplicated_states <- function(data, keep.last = TRUE) {
  out <- do.call(rbind, by(data, data$id, remove_duplicated_states.intern, keep.last))
  row.names(out) <- NULL

  out
}

# @author Quentin Grimonprez
remove_duplicated_states.intern <- function(data, keep.last = TRUE) {
  data <- data[order(data$time), ]

  outRle <- rle(as.character(data$state))
  indToKeep <- 1 + c(0, cumsum(outRle$lengths[-length(outRle$lengths)]))

  if (keep.last && indToKeep[length(indToKeep)] != nrow(data)) {
    indToKeep <- c(indToKeep, nrow(data))
  }

  return(data[indToKeep, ])
}


#' Convert a matrix to a cfda data.frame
#'
#' @param X matrix containing the states
#' @param times time values. If \code{NULL}, it uses a sequence of integers starting with 1
#' @param labels id labels. If \code{NULL}, it uses the matrix colnames
#' @param byrow if \code{FALSE}, one column = one trajectory
#'
#' @return a data.frame in the cfda format
#'
#' @examples
#' x <- matrix(c("a", "b", "c", "c",
#'               "c", "a", "a", "a",
#'               "b", "c", "a", "b"), ncol = 4, byrow = TRUE,
#'               dimnames = list(NULL, paste0("ind", 1:4)))
#' matrixToCfd(x)
#' @export
matrixToCfd <- function(X, times = NULL, labels = NULL, byrow = FALSE) {
  checkLogical(byrow, "byrow")
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or a data.frame")
  }
  nTimes <- ifelse(byrow, ncol(X), nrow(X))
  nInd <- ifelse(byrow, nrow(X), ncol(X))

  if (is.null(times)) {
    times <- seq_len(nTimes)
  } else {
    if (!is.numeric(times) || !((is.vector(times) && (length(times) == nTimes)) ||
                                (is.matrix(times) && (length(times) == nTimes)) ||
                                (is.matrix(times) && (length(times) == nTimes * nInd)))) {
      stop(paste0("times must be a numeric vector of length ", nTimes, " or a matrix of length ", nTimes, "x", nInd))
    }
  }

  if (byrow) {
    X <- t(X)
    if (is.vector(times)) {
      times <- matrix(times, ncol = 1)
    } else {
      times <- t(times)
    }
  }

  if (is.vector(times)) {
    times <- matrix(times, ncol = 1)
  }

  timesPerInd <- ncol(times) > 1

  if (is.null(labels)) {
    if (!is.null(colnames(X))) {
      labels <- colnames(X)
    } else {
      labels <- seq_len(ncol(X))
    }
  } else if (!is.vector(labels) || (length(labels) != nInd)) {
      stop(paste0("labels must be a vector of length ", nInd))
  }

  outData <- data.frame(id = c(), time = c(), state = c())

  for (ind in seq_len(nInd)) {
    indT <- ifelse(timesPerInd, ind, 1)
    outData <- rbind(outData, data.frame(id = labels[ind], time = times[1, indT], state = X[1, ind]))
    if (nTimes > 2) {
      # do not copy 2 consecutive time values with the same state
      for (time in 2:(nTimes - 1)) {
        if (X[time, ind] != X[time - 1, ind]) {
          outData <- rbind(outData, data.frame(id = labels[ind], time = times[time, indT], state = X[time, ind]))
        }
      }
    }
    outData <- rbind(outData, data.frame(id = labels[ind], time = times[nTimes, indT], state = X[nTimes, ind]))
  }

  rownames(outData) = NULL
  return(outData)
}

# convert a fd object to a categorical functional data frame (see convertToCfd)
fdToCfd <- function(fd, breaks, labels = NULL, include.lowest = FALSE, right = TRUE, times = NULL, idLabels = NULL, nx = 200) {
  if (!inherits(fd, "fd")) {
    stop("fd is not a fd object")
  }
  if (is.null(times)) {
    times <- seq(fd$basis$rangeval[1], fd$basis$rangeval[2], length = nx)
  }
  if (is.null(idLabels)) {
    idLabels <- fd$fdnames$reps
  }
  X <- eval.fd(times, fd)
  return(quantiMatrixToCfd(X, breaks, labels = labels, include.lowest = include.lowest, right = right,
                          idLabels = idLabels, times = times, byrow = FALSE))
}

# convert a qualitative matrix to a categorical functional data frame (see convertToCfd)
quantiMatrixToCfd <- function(X, breaks, labels = NULL, include.lowest = FALSE, right = TRUE,
                             times = NULL, idLabels = NULL, byrow = FALSE) {
  X <- matrix(cut(X, breaks = breaks, labels = labels, right = right, include.lowest = include.lowest),
              nrow = nrow(X), dimnames = dimnames(X))
  if (any(is.na(X))) {
    stop("The conversion has generated NA. Please, correct your breaks.")
  }
  return(matrixToCfd(X, times, idLabels, byrow))
}

#' Convert data to categorical functional data
#'
#' @param x matrix or fd object
#' @param breaks either a numeric vector of two or more unique cut points or a single number (greater than or equal to 2)
#' giving the number of intervals into which x is to be cut.
#' @param labels labels for the levels of the resulting category. By default, labels are constructed using "(a,b]"
#' interval notation. If labels = FALSE, simple integer codes are returned instead of a factor.
#' @param include.lowest logical, indicating if an ‘x[i]’ equal to the lowest (or highest, for right = FALSE) ‘breaks’ value should be included.
#' @param right logical, indicating if the intervals should be closed on the right (and open on the left) or vice versa.
#' @param times vector containing values at which \code{fd} is to be evaluated
#' @param idLabels vector containing id labels. If NULL it use the names found in the matrix or fd object
#' @param nx Only if \code{x} is a fd object. Number of points to evaluate \code{fd}
#' @param byrow Only if \code{x} is a matrix. If \code{FALSE}, one column = one trajectory
#'
#' @return a data.frame in the cfda format
#'
#' @examples
#' # fd object
#' data("CanadianWeather")
#' temp <- CanadianWeather$dailyAv[,, "Temperature.C"]
#' basis <- create.bspline.basis(c(1, 365), nbasis = 8, norder = 4)
#' fd <- smooth.basis(1:365, temp, basis)$fd
#'
#' # "Very Cold" = [-50:-10), "Cold" = [-10:0), ...
#' out <- convertToCfd(fd, breaks = c(-50, -10, 0, 10, 20, 50),
#'                     labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"),
#'                     times = 1:365)
#'
#' # matrix
#' out2 <- convertToCfd(temp, breaks = c(-50, -10, 0, 10, 20, 50),
#'                      labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"),
#'                      times = 1:365, byrow = FALSE)
#'
#' @export
convertToCfd <- function(x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE, times = NULL,
                         idLabels = NULL, nx = 200, byrow = FALSE) {
  if (inherits(x, "fd")) {
    return(fdToCfd(x, breaks, labels, include.lowest, right, times, idLabels = NULL, nx))
  } else if (is.matrix(x) || is.data.frame(x)) {
    return(quantiMatrixToCfd(x, breaks, labels, include.lowest, right, times, idLabels = NULL, byrow))
  }
}
