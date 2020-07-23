

bootstrapEncoding <- function(data, basisobj, nBootstrap = 50, propBootstrap = 0.5, nCores = max(1, ceiling(detectCores()/2)), verbose = TRUE, ...)
{
  t1 <- proc.time()
  ## check parameters
  checkData(data)
  checkDataBeginTime(data)
  checkDataEndTmax(data)
  if(!is.basis(basisobj))
    stop("basisobj is not a basis object.")
  if(any(is.na(nCores)) || (length(nCores) != 1) || !is.whole.number(nCores) || (nCores < 1))
    stop("nCores must be an integer > 0.")
  if(any(is.na(nBootstrap)) || (length(nBootstrap) != 1) || !is.whole.number(nBootstrap) || (nBootstrap < 1))
    stop("nBootstrap must be an integer > 0.")
  if(any(is.na(propBootstrap)) || !is.numeric(propBootstrap) || (length(propBootstrap) != 1) || (propBootstrap > 1) || (propBootstrap < 0))
    stop("nCores must be an integer > 0.")
  ## end check
  
  if(verbose)
    cat("######### Compute encoding #########\n")
  
  # change state as integer
  out <- stateToInteger(data$state)
  data$state = out$state
  label <- out$label
  rm(out)
  
  # refactor labels as 1:nbId
  uniqueId <- unique(data$id)
  nId <- length(uniqueId)
  id2 <- refactorCategorical(data$id, uniqueId, seq_along(uniqueId)) 
  uniqueId2 <- unique(id2)
  
  sizeBootstrap <- floor(propBootstrap * nId)
  
  nCores <- min(max(1, nCores), detectCores()-1)
  
  Tmax <- max(data$time)
  K <- length(label$label)
  nBasis <- basisobj$nbasis
  
  if(verbose)
  {
    cat(paste0("Number of individuals: ", nId, "\n"))
    cat(paste0("Number of states: ", K, "\n"))
    cat(paste0("Basis type: ", basisobj$type, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
    cat(paste0("Number of cores: ", nCores, "\n"))
  }
  
  
  V <- computeVmatrix(data, uniqueId2, id2, basisobj, K, nCores, verbose, ...)
  
  Uval <- computeUmatrix(data, uniqueId2, id2, basisobj, K, nCores, verbose, ...)
  
  outEnc <- list()
  t3 <- proc.time()
  if(verbose)
    cat("---- Compute encoding:\n")
  for(i in 1:nBootstrap)
  {
    cat("*")
    idToKeep <- sample(nId, sizeBootstrap, replace = TRUE)
    
    outEnc[[i]] = computeEncoding(Uval[idToKeep, ], V[idToKeep, ], K, nBasis, idToKeep, label, verbose = FALSE)
    
    outEnc[[i]] = c(outEnc[[i]] , list(basisobj = basisobj))
    class(outEnc[[i]]) = "fmca"
  }
  
  # reorder alpha, pc such that representation have the same sign for each bootstrap sample
  outEnc = unifySign(outEnc)

  
  alpha <- lapply(seq_along(outEnc[[1]]$alpha), function(harm) Reduce("+", lapply(outEnc, function(x) x$alpha[[harm]]))/nBootstrap)
  invF05vec <- computeMeanInvF05vec(outEnc)
  
  t4 <- proc.time()
  if(verbose)
    cat(paste0("\nDONE in ", round((t4-t3)[3], 2), "s\n"))
  
  out <- c(list(out = outEnc, alpha = alpha, pc = V%*%invF05vec, V = V, invF05vec = invF05vec, basisobj = basisobj, nBootstrap = nBootstrap))
  class(out) = c("fmca", "fmcaBootstrap")
  
  t2 <- proc.time()
  
  if(verbose)
    cat(paste0("Run Time: ", round((t2-t1)[3], 2), "s\n"))
  
  return(out)
}

# eigenvectors are equivalent up to the sign
# we try to have the same sign for each bootstrap sample
unifySign <- function(out)
{
  signNegMax <- matrix(nrow = length(out[[1]]$alpha), ncol = length(out))
  for(i in seq_along(out))
  {
    for(j in seq_along(out[[i]]$alpha))
    {
      signNegMax[j, i] <- out[[i]]$alpha[[j]][which.max(abs(out[[i]]$alpha[[j]]))] < 0
    }
  }
  
  for(i in seq_along(out))
  {
    for(j in seq_along(out[[i]]$alpha))
    {
      if(signNegMax[j, i])
      {
        out[[i]]$alpha[[j]] = out[[i]]$alpha[[j]] * -1
        out[[i]]$pc[, j] = out[[i]]$pc[, j] * -1
        out[[i]]$invF05vec[, j] = out[[i]]$invF05vec[, j] * -1      
      }  
    }
  }

  return(out)
}


# compute the mean pc scores
#
# compute pc scores for each bootstrap sample then compute the mean
computeMeanInvF05vec <- function(encodingList)
{
  vec <- list()
  for(i in seq_along(encodingList))
  {
    vec[[i]] = encodingList[[i]]$invF05vec 
  }
  
  return(Reduce("+", vec)/length(encodingList))
}

