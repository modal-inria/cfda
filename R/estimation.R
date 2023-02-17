
#' Estimate transition matrix and spent time
#'
#' Calculates crude initial values for transition intensities by assuming that the data represent the exact transition times of the Markov process.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state.
#'
#' @return list of two elements: \code{Q}, the estimated transition matrix, and \code{lambda}, the estimated time spent in each state
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 100, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # estimation
#' mark <- estimate_Markov(d_JK)
#' mark$P
#' mark$lambda
#' @seealso \link{plot.Markov}
#'
#' @author Cristian Preda
#'
#' @export
estimate_Markov <- function(data) {
  ## check parameters
  checkData(data)
  ## end check

  # il faut supprimer les sauts dans les mêmes états
  data <- remove_duplicated_states(data, keep.last = TRUE)

  P_est <- prop.table(statetable(data, removeDiagonal = TRUE), margin = 1)

  # estimation of the time spent in each state
  T_est <- estimateT(data)
  lambda_est <- rep(NaN, ncol(P_est))
  names(lambda_est) <- colnames(P_est)
  lambda_est[names(T_est)] <- 1 / T_est

  out <- list(P = P_est, lambda = lambda_est)
  class(out) <- "Markov"

  return(out)
}


# estimate the mean time spent in each state
#
# @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
# \code{state}, associated state.
#
estimateT <- function(data) {
  tSpentInState <- by(
    data, data$id,
    function(x) {
      t <- diff(x$time)
      s <- x$state[seq_along(t)]
      data.frame(duration = t, state = s)
    }
  )
  tSpentInState <- do.call(rbind, tSpentInState)

  t <- tapply(tSpentInState$duration, tSpentInState$state, mean)

  t
}



# identifier les états d'où on part jamais dans [0,T] (absorbants ou pas)
# add row (containing 0s) for the missing state
#
# @param aux output of statetable.msm
# @author Cristian Preda
completeStatetable <- function(aux) {

  aux1 <- matrix(0, nrow = length(colnames(aux)), ncol = length(colnames(aux)),
                 dimnames = list(colnames(aux), colnames(aux)))

  aux1[match(rownames(aux), colnames(aux)), ] <- aux

  return(aux1)
}


#' Plot the transition graph
#'
#' Plot the transition graph between the different states. A node corresponds to a state with the mean
#' time spent in this state. Each arrow represents the probability of transition between states.
#'
#' @param x output of \code{\link{estimate_Markov}} function
#' @param ... parameters of \code{plotmat} function from \code{diagram} package (see details).
#'
#' @return No return value, called for side effects
#'
#' @details
#' Some useful extra parameters:
#' \itemize{
#'   \item \code{main} main title.
#'   \item \code{dtext} controls the position of arrow text relative to arrowhead (default = 0.3).
#'   \item \code{relsize}	scaling factor for size of the graph (default = 1).
#'   \item \code{box.size} size of label box, one value or a vector with dimension = number of rows of \code{x$P}.
#'   \item \code{box.cex}	relative size of text in boxes, one value or a vector with dimension=number of rows of \code{x$P}.
#'   \item \code{arr.pos} relative position of arrowhead on arrow segment/curve.
#' }
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 100, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # estimation
#' mark <- estimate_Markov(d_JK)
#'
#' # transition graph
#' plot(mark)
#' @author Cristian Preda
#'
#' @export
plot.Markov <- function(x, ...) {
  extraParam <- list(...)
  defaultParam <- list(
    A = t(round(x$P, 2)), main = "The transition graph", box.prop = 0.3,
    box.col = "yellow", arr.length = 0.2, shadow.size = 0,
    name = paste0(
      colnames(x$P), rep(" (", ncol(x$P)), round(1 / x$lambda, 2),
      rep(")", ncol(x$P))
    )
  )

  param <- c(extraParam, defaultParam[which(!(names(defaultParam) %in% names(extraParam)))])

  do.call(plotmat, param)
}
