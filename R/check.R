# Check if the data.frame has the required format
# @author Quentin Grimonprez
checkData <- function(data, checkNrows = TRUE) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  requiredColNames <- c("id", "time", "state")
  missColNames <- !(requiredColNames %in% colnames(data))
  if (any(missColNames)) {
    stop(paste0("Missing columns in data: ", paste(requiredColNames[missColNames], collapse = ", "), "."))
  }

  if (checkNrows && (nrow(data) <= 1)) {
    stop("There is only one row or less.")
  }

  if (any(is.na(data))) {
    stop("There is some missing values.")
  }

  invisible(return(NULL))
}


# Check if all individuals end with the same time value
# @author Quentin Grimonprez
checkDataEndTmax <- function(data) {
  lastTime <- tapply(data$time, data$id, function(x) x[length(x)])

  nLastTime <- length(unique(lastTime))

  if (nLastTime != 1) {
    stop("Each individual must end with the same time value.")
  }

  invisible(NULL)
}


# Check if all individuals start with the same time value
# @author Quentin Grimonprez
checkDataBeginTime <- function(data) {
  firstTime <- tapply(data$time, data$id, function(x) x[1])

  nFirstTime <- length(unique(firstTime))

  if (nFirstTime != 1) {
    stop("Each individual must begin with the same time value.")
  }

  invisible(NULL)
}

# Check if all individual has different time values
# @author Quentin Grimonprez
checkDataNoDuplicatedTimes <- function(data) {
  duplicatedTimes <- any(tapply(data$time, data$id, function(x) any(duplicated(x))))

  if (duplicatedTimes) {
    warning("Some ids contain duplicated time values.")
  }

  invisible(NULL)
}

# Check if the given parameter is a single boolean
# @author Quentin Grimonprez
checkLogical <- function(x, paramName) {
  if (length(x) != 1) {
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  }

  if (is.na(x)) {
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  }

  if (!is.logical(x)) {
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  }

  invisible(NULL)
}

# Check if it is an integer (or vector of integer)
# @author Quentin Grimonprez
is.whole.number <- function(x) {
  x == as.integer(x)
}
