
#' Plot categorical functional data
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
#' @param group vector, of the same length as the number individuals of \code{data}, containing group index.
#' Groups are displayed on separate plots.
#' If \code{group = NA}, the corresponding individuals in \code{data} is ignored.
#' @param col a vector containing color for each state (can be named)
#' @param addId If TRUE, add id labels
#' @param addBorder If TRUE, add black border to each individual
#' @param sort If TRUE, id are sorted according to the duration in their first state
#' @param nCol number of columns when \code{group} is given
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#' On the plot, each row represents an individual over [0:Tmax].
#' The color at a given time gives the state of the individual.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # add a line with time Tmax at the end of each individual
#' d_JKT <- cut_data(d_JK, Tmax = 10)
#'
#' plotData(d_JKT)
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' plotData(d_JKT, col = c("red", "blue", "green", "brown")) +
#'   labs(title = "Trajectories of a Markov process")
#'
#'
#' # use the group variable: create a group with the 3 first variables and one with the others
#' group <- rep(1:2, c(3, 7))
#' plotData(d_JKT, group = group)
#'
#'
#' # use the group variable: remove the id number 5 and 6
#' group[c(5, 6)] <- NA
#' plotData(d_JKT, group = group)
#' @author Cristian Preda, Quentin Grimonprez
#' @family Descriptive statistics
#' @export
plotData <- function(data, group = NULL, col = NULL, addId = TRUE, addBorder = TRUE, sort = FALSE, nCol = NULL) {
  ## check parameters
  checkData(data)
  checkLogical(addId, "addId")
  checkLogical(addBorder, "addBorder")
  checkLogical(sort, "sort")
  if (!is.null(group) && ((!is.vector(group) && !is.factor(group)) || (length(group) != length(unique(data$id))))) {
    stop("group must be a vector with the same length than the number of ids of data.")
  }
  if (!is.null(nCol) && (!is.numeric(nCol) || (length(nCol) != 1) || !is.whole.number(nCol) || (nCol < 1))) {
    stop("nCol must be an integer > 0.")
  }
  ## end check

  if (!is.null(group)) {
    data$group <- rep(NA, nrow(data))
    idNames <- unique(data$id)
    for (i in seq_along(idNames)) {
      data$group[data$id == idNames[i]] <- group[i]
    }
    data <- data[!is.na(data$group), ]
  }

  d_graph <- rep_large_ind(data)

  if (!is.null(group)) {
    d_graph <- d_graph[order(d_graph$group), ]
  }

  # to be sure that state are considered a qualitative
  # if already a factor, we do not execute this in order to not drop unused levels
  if (!is.factor(d_graph$state)) {
    d_graph$state <- factor(d_graph$state)
  }

  nInd <- length(unique(d_graph$id))

  if (is.null(group)) {
    d_graph$position <- computePosition(data, d_graph$id, sort)
  } else {
    d_graph$position <- computePositionPerGroup(data, d_graph$id, d_graph$group, sort)
  }

  p <- ggplot() +
    geom_rect(data = d_graph,
              mapping = aes_string(xmin = "t_start", xmax = "t_end", ymin = "position - 0.5",
                                   ymax = "position + 0.5", fill = "state"),
              color = ifelse(addBorder, "black", NA)) +
    scale_x_continuous(name = "Time") +
    labs(fill = "State")

  if (!is.null(group)) {
    p <- p + facet_wrap("group", scales = "free_y", labeller = labeller(.default = createLabeller(group)), ncol = nCol)
  }


  if (addId) {
    p <- p + scale_y_continuous(name = "Id", breaks = seq_len(nInd), labels = unique(d_graph$id)[order(unique(d_graph$position))])
  } else {
    p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }

  if (!is.null(col)) {
    p <- p + scale_fill_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # do not remove unused labels

  return(p)
}


# transform the data format to a new format with 4 columns: id, t_stat, t_end, state.
# usefull for ggplot
# @author Cristian Preda
rep_large_ind <- function(data) {
  out <- by(data, data$id, function(x) {
    d <- data.frame(
      id = x$id[seq_len(nrow(x) - 1)],
      t_start = x$time[seq_len(nrow(x) - 1)],
      t_end = x$time[2:nrow(x)],
      state = x$state[seq_len(nrow(x) - 1)], stringsAsFactors = FALSE
    )

    if ("group" %in% names(data)) {
      d$group <- x$group[seq_len(nrow(x) - 1)]
    }

    return(d)
  })

  return(do.call(rbind, out[unique(data$id)]))
}


computePosition <- function(data, id, sort = FALSE) {
  if (sort) {
    # order according to first state duration
    b <- orderFirstState(data)

    ord <- match(id, b$id)
    position <- ord # position of id on the y-axis
  } else {
    position <- unclass(factor(id, levels = unique(id))) # return integers associated with the different ids (labels from a factor)
  }

  return(position)
}


computePositionPerGroup <- function(data, id, group, sort = FALSE) {
  pos <- list()
  groupName <- unique(group)
  for (i in seq_along(groupName)) {
    pos[[i]] <- computePosition(data[data$group == groupName[i], ], id[group == groupName[i]], sort = sort)
  }

  if (length(pos) > 1) {
    for (i in 2:length(pos)) {
      pos[[i]] <- pos[[i]] + max(pos[[i - 1]])
    }
  }

  unlist(pos)
}


# compute the duration of the first state per individual
# and return an ordered data.frame per first state and duration
# @author Quentin Grimonprez
orderFirstState <- function(data) {
  firstState <- do.call(rbind, by(data, data$id, function(x) {
    data.frame(id = x$id[1], time = ifelse(length(x$time) < 2, Inf, x$time[2] - x$time[1]),
               state = x$state[1], stringsAsFactors = FALSE)
  }))
  firstStateOrdered <- do.call(rbind, by(firstState, firstState$state, function(x) {
    x[order(x$time), ]
  }))
}


# function to label group in facet_wrap
# group name: n=effectif
createLabeller <- function(group) {
  part <- table(group)
  labelGroup <- as.list(paste0(names(part), ": n=", part))
  names(labelGroup) <- names(part)

  group_labeller <- function(variable, value) {
    return(labelGroup[value])
  }

  return(group_labeller)
}


#' @title Summary
#'
#' @description Get a summary of the data.frame containing categorical functional data
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
#' @param max.print maximal number of states to display
#'
#' @return a list containing:
#' \itemize{
#'   \item \code{nRow} number of rows
#'   \item \code{nInd} number of individuals
#'   \item \code{timeRange} minimal and maximal time value
#'   \item \code{uniqueStart} TRUE, if all individuals have the same time start value
#'   \item \code{uniqueEnd} TRUE, if all individuals have the same time start value
#'   \item \code{states} vector containing the different states
#'   \item \code{visit} number of individuals visiting each state
#' }
#'
#' @examples
#' data(biofam2)
#' summary_cfd(biofam2)
#' @author Quentin Grimonprez
#' @family Descriptive statistics
#' @export
summary_cfd <- function(data, max.print = 10) {
  checkData(data)

  nIndiv <- length(unique(data$id))
  if (is.factor(data$state)) {
    states <- levels(data$state)
  } else {
    states <- as.character(sort(unique(data$state)))
  }
  nState <- length(states)

  timeRange <- range(data$time)

  timeRangeInd <- do.call(rbind, tapply(data$time, data$id, range))
  sameStart <- (length(unique(timeRangeInd[, 1])) == 1)
  sameEnd <- (length(unique(timeRangeInd[, 2])) == 1)

  nRow <- nrow(data)

  nIndVisitingState <- tapply(data$id, data$state, function(x) length(unique(x)))

  cat("Number of rows:", nRow, "\n")
  cat("Number of individuals:", nIndiv, "\n")
  cat("Time Range:", timeRange[1], "-", timeRange[2], "\n")
  cat("Same time start value for all ids:", sameStart, "\n")
  cat("Same time end value for all ids:", sameEnd, "\n")
  cat("Number of states:", nState, "\n")
  cat("States:\n  ")
  cat(paste0(paste(head(states, n = 10), collapse = ", "), ifelse(nState > 10, ", ...", "")))
  cat("\n")
  cat("Number of individuals visiting each state:\n")
  print(head(nIndVisitingState, n = 10))

  invisible(list(
    nRow = nRow, nInd = nIndiv, timeRange = timeRange, uniqueStart = sameStart, uniqueEnd = sameEnd,
    states = states, visit = nIndVisitingState
  ))
}
