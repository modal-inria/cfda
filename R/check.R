# Check if the data.frame has the required format
# @author Quentin Grimonprez
checkDataMsm <- function(data_msm)
{
  if(!is.data.frame(data_msm))
    stop("data_msm must be a data.frame.")
    
  requiredColNames <- c("id", "time", "state")
  missColNames <-  !(requiredColNames %in% colnames(data_msm))
  if(any(missColNames))
    stop(paste0("Missing columns in data_msm: ", paste(requiredColNames[missColNames], collapse = ", "), "."))
  
  if(nrow(data_msm) <= 1)
    stop("There is only one row or less.")
  
  if(any(is.na(data_msm)))
    stop("There is some missing values.")
  
  if(any(!is.whole.number(data_msm$state)) || any(data_msm$state <= 0))
    stop("state must be strictly positive integer.")
  
  invisible(return(NULL))
}


# Check if each individual ends by the same time value
# @author Quentin Grimonprez
checkDataEndTmax <- function(data_msm)
{
  lastTime <- tapply(data_msm$time, data_msm$id, function(x) x[length(x)])
  
  nLastTime <- length(unique(lastTime))
  
  if(nLastTime != 1)
    stop("Each individual must finished by the same time value.")
  
  invisible(return(NULL))
}

# Check if the given parameter is a single boolean
# @author Quentin Grimonprez
checkLogical <- function(x, paramName)
{
  if(length(x) != 1)
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  
  if(!is.logical(x))
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  
  invisible(return(NULL))
}

# Check if it is an integer (or vector of integer)
# @author Quentin Grimonprez
is.whole.number <- function(x)
{
  x == as.integer(x)
}
