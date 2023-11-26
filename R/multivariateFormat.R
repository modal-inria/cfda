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

