#' Cut data to a maximal given time
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
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
refactorCategorical <- function(data, oldCateg = unique(data), newCateg = 1:length(oldCateg)) {
  ind <- match(data, oldCateg)

  if (any(is.na(ind[!is.na(data)]))) {
    warning("NA produced.")
  }

  return(newCateg[ind])
}


#' Remove duplicated states
#'
#' Remove duplicated consecutive states from data. If for an individual there is two or more consecutive states that are identical,
#' only the first is kept. Only time when the state changes are kept.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
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
