
# compute bootstrap estimates of encoding
#
# @param Uval output of computeUmatrix function
# @param V output of computeVmatrix function
# @param K number of states
# @param nBasis number of basis
# @param label label element from stateToInteger function
# @param nId number of individuals
# @param propBootstrap proportion of individuals used for bootstrap sample
# @param nBootstrap number of bootstrap samples
# @param output of getSignReference functionon the ful encoding
# @param verbose if TRUE print some information
#
# @return a list of computeEncoding function output
#
# @author Quentin Grimonprez
computeBootStrapEncoding <- function(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, signReference, verbose)
{
  outEnc <- list()
  t3 <- proc.time()
  if(verbose)
    cat("---- Compute Bootstrap Encoding:\n")
  
  for(i in 1:nBootstrap)
  {
    if(verbose)
      cat("*")
    idToKeep <- sample(nId, floor(propBootstrap * nId), replace = TRUE)
    
    try({outEnc[[i]] = computeEncoding(Uval[idToKeep, ], V[idToKeep, ], K, nBasis, idToKeep, label, verbose = FALSE)})
    
    
    # outEnc[[i]] = c(outEnc[[i]] , list(basisobj = basisobj))
    # class(outEnc[[i]]) = "fmca"
  }
  
  # reorder alpha, pc such that representation have the same sign for each bootstrap sample
  outEnc = unifySign(outEnc, signReference)
  
  # bootstrap estimate of alpha and invF05
  alpha <- lapply(1:(K * nBasis), function(harm) Reduce("+", lapply(outEnc, function(x) x$alpha[[harm]]))/nBootstrap)
  invF05vec <- computeMeanInvF05vec(outEnc)
  
  t4 <- proc.time()
  if(verbose)
    cat(paste0("\nDONE in ", round((t4-t3)[3], 2), "s\n"))
  
  return(outEnc)
}

# get the sign of the alpha of the full encoding to later ensure that bootstrap samples have the same sign
# @author Quentin Grimonprez
getSignReference <- function(alpha)
{
  pos  <- rep(0, length(alpha)) 
  isNeg <- rep(FALSE, length(alpha))
  for(i in seq_along(alpha))
  {
    pos[i] = which.max(abs(alpha[[i]]))
    isNeg[i] = alpha[[i]][pos[i]] < 0
  }
  
  return(list(position = pos, isNegative = isNeg))
}


# eigenvectors are equivalent up to the sign
# we try to have the same sign for each bootstrap sample
#
# @param out a list of computeEncoding function output
#
# @author Quentin Grimonprez
unifySign <- function(out, signReference)
{
  signNeg <- matrix(nrow = length(signReference$pos), ncol = length(out))
  for(i in seq_along(out))
  {
    if(!is.null(out[[i]])) # an output can be NULL due to inversion problem 
    {
      # we look if the element at the given position is negative
      for(j in seq_along(out[[i]]$alpha))
      {
        signNeg[j, i] <- out[[i]]$alpha[[j]][signReference$position[j]] < 0
        
        if(signNeg[j, i] != signReference$isNegative[j])
        {
          out[[i]]$alpha[[j]] = out[[i]]$alpha[[j]] * -1
          out[[i]]$pc[, j] = out[[i]]$pc[, j] * -1
          out[[i]]$invF05vec[, j] = out[[i]]$invF05vec[, j] * -1      
        }  
      }
    }
  }
  
  
  return(out)
}


# compute the mean pc scores
#
# compute pc scores for each bootstrap sample then compute the mean
#
# @param bootEncoding output of computeBootStrapEncoding function
#
# @author Quentin Grimonprez
computeMeanInvF05vec <- function(bootEncoding)
{
  vec <- list()
  for(i in seq_along(bootEncoding))
  {
    vec[[i]] = bootEncoding[[i]]$invF05vec 
  }
  
  return(Reduce("+", vec)/length(bootEncoding))
}


# Compute the variance of alpha
#
# @param bootEncoding output of computeBootStrapEncoding function
#
# @return a list (length number of harmonics) of list (length number of states) of variance matrix
#
# @author Quentin Grimonprez
computeVarianceAlpha <- function(bootEncoding)
{
  nHarm <- min(sapply(bootEncoding, function(x) length(x$alpha)))
  
  varAlpha <- list()
  for(harm in 1:nHarm)
  {
    nState <- ncol(bootEncoding[[1]]$alpha[[harm]])
    varAlpha[[harm]] <- list()
    for(iState in 1:nState)
    {
      varAlpha[[harm]][[iState]] <- cov(t(sapply(bootEncoding, function(x) x$alpha[[harm]][, iState])))
    }
  }
  
  return(varAlpha)
}


# compute the variance of encoding
#
# a_x(t) = sum_i alpha_ix * phi_i(t) 
# Var(a_x(t)) = sum_i var(alpha_ix) * phi_i^2(t) + sum_{i<j}  2 * phi_i(t) * phi_j(t) * cov(alpha_ix, alpha_jx)
#
# @author Quentin Grimonprez
computeVarianceEncoding <- function(varAlpha, basisobj, harm = 1, nx = 200)
{
  nBasis <- basisobj$nbasis
  phi <- fd(diag(nBasis), basisobj)
  nState <- length(varAlpha[[harm]])
  
  timeVal <- seq(basisobj$rangeval[1], basisobj$rangeval[2], length = nx)
  
  Phi <- matrix(nrow = length(timeVal), ncol = nBasis)
  for(i in 1:nBasis)
    Phi[, i] = eval.fd(timeVal, phi[i])
  
  funcVar <- list()
  for(iState in 1:nState)
  {
    funcVar[[iState]] = rep(NA, nx)
    for(j in seq_along(timeVal))
    {
      funcVar[[iState]][j] <- Phi[j, , drop = FALSE] %*% varAlpha[[harm]][[iState]] %*% t(Phi[j, , drop = FALSE])
    }
  }
  
  return(funcVar)
}

