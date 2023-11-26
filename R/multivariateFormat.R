#' Convert a list of cfd into a multivariate cfd
#'
#' @param x list of cfd
#' @param stateColumns column names for multivariate states. By default, "state1", "state2", ...
#'
#' @return a multivariate cfd
#'
#' @examples
#' set.seed(42)
#' x1 <- generate_Markov(n = 10, K = 2)
#' x1 <- cut_data(x1, Tmax = 1)
#' x2 <- generate_Markov(n = 10, K = 2)
#' x2 <- cut_data(x2, Tmax = 1)
#'
#' x <- list(x1, x2)
#'
#' mvcfd <- convert2mvcfd(x)
#'
#' @export
convert2mvcfd <- function(x, stateColumns = NULL) {
  ### check data
  if (!inherits(x, "list")) {
    stop("data must be a list of data.frames")
  }

  for (elem in x) {
    checkData(elem)
    checkDataNoDuplicatedTimes(elem)
  }
  ### end check data

  nDim <- length(x)
  if (is.null(stateColumns)) {
    stateColumns <- paste0("state", seq_len(nDim))
  }

  # rename state columns to each data frame
  for (i in seq_len(nDim)) {
    x[[i]] <- x[[i]] %>% rename(!!stateColumns[i] := "state")
  }

  # merge data frames
  x <- Reduce(function(x1, x2) merge(x1, x2, by = c("id", "time"), all = TRUE), x)

  # order by id and time
  x <- arrange(x, .data[["id"]], .data[["time"]])

  # fill missing values with the state before
  x <- as.data.frame(x %>% group_by(.data[["id"]]) %>% fill(all_of(stateColumns), .direction = "downup"))

  return(distinct(x))
}

#' Convert a multivariate cfd to an univariate one
#'
#' @param x multivariate cfd: data.frame with id, time columns and several stateColumns (output of \code{\link{convert2mvcfd}})
#' @param sep separator used to concatenate states
#' @param stateColumns names of the state columns. If NULL, it looks for columns with state in their names
#'
#' @return univariate cfd to use in \code{\link{compute_optimal_encoding}}. A data.frame with is, time and state columns
#'
#' @examples
#' set.seed(42)
#' x1 <- generate_Markov(n = 10, K = 2)
#' x1 <- cut_data(x1, Tmax = 1)
#' x2 <- generate_Markov(n = 10, K = 2)
#' x2 <- cut_data(x2, Tmax = 1)
#'
#' x <- list(x1, x2)
#'
#' mvcfd <- convert2mvcfd(x)
#' cfd <- convertMvcfd2cfd(mvcfd)
#'
#' @export
convertMvcfd2cfd <- function(x, sep = "_", stateColumns = NULL) {
  checkData(x, requiredColNames = c("id", "time"))

  if (is.null(stateColumns)) {
    stateColumns <- colnames(x)[grep("state", colnames(x))]
  }

  xCfd <- x[c("id", "time")]
  xCfd$state <- as.factor(do.call(paste, c(x[stateColumns], sep = sep)))

  return(xCfd)
}

#' Convert a list of univariate cfds to an unique univariate cfd
#'
#' @param x list of cfd
#' @param sep separator used to concatenate states
#'
#' @return univariate cfd to use in \code{\link{compute_optimal_encoding}}. A data.frame with is, time and state columns
#'
#' @examples
#' set.seed(42)
#' x1 <- generate_Markov(n = 10, K = 2)
#' x1 <- cut_data(x1, Tmax = 1)
#' x2 <- generate_Markov(n = 10, K = 2)
#' x2 <- cut_data(x2, Tmax = 1)
#'
#' x <- list(x1, x2)
#'
#' cfd <- convertListCfd2Cfd(x)
#'
#' @export
convertListCfd2Cfd <- function(x, sep = "_") {
  xMcfd <- convert2mvcfd(x)
  return(convertMvcfd2cfd(xMcfd, sep = sep))
}
