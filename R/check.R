# Check if the data.frame has the required format
# @author Quentin Grimonprez
checkData <- function(data, minSize = 1) {
  if (!is.data.frame(data) && !is_tibble(data)) {
    stop("data must be a data.frame.")
  }

  requiredColNames <- c("id", "time", "state")
  missColNames <- !(requiredColNames %in% colnames(data))
  if (any(missColNames)) {
    stop(paste0("Missing columns in data: ", paste(requiredColNames[missColNames], collapse = ", "), "."))
  }

  if (nrow(data) <= minSize) {
    stop(paste0("There is ", minSize, " row or less."))
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
  errorMessage <- paste0(paramName, " must be either TRUE or FALSE.")
  if (length(x) != 1) {
    stop(errorMessage)
  }

  if (is.na(x)) {
    stop(errorMessage)
  }

  if (!is.logical(x)) {
    stop(errorMessage)
  }

  invisible(NULL)
}

# Check if it is an integer (or vector of integer)
# @author Quentin Grimonprez
is.whole.number <- function(x) {
  is.numeric(x) & (x == as.integer(x))
}

checkInteger <- function(
    x, minValue = -Inf, maxValue = +Inf, acceptNULL = FALSE, minEqual = FALSE, maxEqual = FALSE,
    paramName = "x", customMessage = NULL) {
  errorMessage <- paste(paramName, "must be an integer")
  addAnd <- FALSE
  if (minValue > -Inf) {
    errorMessage <- paste(errorMessage, ifelse(minEqual, ">=", ">"), minValue)
    addAnd <- TRUE
  }
  if (maxValue < +Inf) {
    if (addAnd) {
      errorMessage <- paste(errorMessage, "and")
    }
    errorMessage <- paste(errorMessage, ifelse(minEqual, "<=", "<"), maxValue)
  }
  errorMessage <- paste0(errorMessage, ".")

  is_not_integer <- any(is.na(x)) || (length(x) != 1) || !is.numeric(x) || !is.whole.number(x)
  if (minEqual) {
    is_not_integer <- is_not_integer || (x < minValue)
  } else {
    is_not_integer <- is_not_integer || (x <= minValue)
  }
  if (maxEqual) {
    is_not_integer <- is_not_integer || (x > maxValue)
  } else {
    is_not_integer <- is_not_integer || (x >= maxValue)
  }

  if (acceptNULL) {
    is_not_integer <- !is.null(x) && is_not_integer
  } else {
    is_not_integer <- is.null(x) || is_not_integer
  }

  if (is_not_integer) {
    stop(ifelse(is.null(customMessage), errorMessage, customMessage))
  }
}


checkFmca <- function(x, paramName = "x") {
  if (!inherits(x, "fmca")) {
    stop(paste0(paramName, " must be a fmca object."))
  }
}
