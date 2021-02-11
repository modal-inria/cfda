#' Compute the optimal encoding for each state
#'
#' Compute the optimal encoding for categorical functional data using an extension of the multiple correspondence analysis to a stochastic process.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state. All individuals must begin at the same time T0 and end at the same time Tmax (use \code{\link{cut_data}}).
#' @param basisobj basis created using the \code{fda} package (cf. \code{\link{create.basis}}).
#' @param computeCI if TRUE, perform a bootstrap to estimate the variance of encoding's coefficients
#' @param nBootstrap number of bootstrap samples
#' @param propBootstrap size of bootstrap samples relative to the number of individuals: propBootstrap * number of individuals 
#' @param nCores number of cores used for parallelization. Default is the half of cores.
#' @param verbose if TRUE print some information
#' @param ... parameters for \code{\link{integrate}} function (see details).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{eigenvalues} eigenvalues
#'   \item \code{alpha} optimal encoding coefficients associated with each eigenvectors
#'   \item \code{pc} principal components
#'   \item \code{F} matrix containing the \eqn{F_{(x,i)(y,j)}}
#'   \item \code{V} matrix containing the \eqn{V_{(x,i)}}
#'   \item \code{G} covariance matrix of \code{V}
#'   \item \code{basisobj} \code{basisobj} input parameter
#'   \item \code{pt} output of \link{estimate_pt} function
#'   \item \code{bootstrap} Only if \code{computeCI = TRUE}. Output of every bootstrap run 
#'   \item \code{varAlpha} Only if \code{computeCI = TRUE}. Variance of alpha parameters
#'   \item \code{runTime} Total elapsed time
#' }
#'
#' @details 
#' See the vignette for the mathematical background: \code{RShowDoc("cfda", package = "cfda")}
#'
#' Extra parameters (\emph{...}) for the \code{\link{integrate}} function can be:
#' \itemize{
#'  \item \emph{subdivisions} the maximum number of subintervals.
#'  \item \emph{rel.tol} relative accuracy requested.
#'  \item \emph{abs.tol} absolute accuracy requested.
#' }	
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement 
#' K <- 4
#' Tmax <- 5
#' PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax,
#'                             labels = c("A", "C", "G", "T"))
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 5
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)
#' summary(encoding)
#' 
#' # plot the optimal encoding
#' plot(encoding)
#' 
#' # plot the two first components
#' plotComponent(encoding, comp = c(1, 2))
#' 
#' # extract the optimal encoding
#' get_encoding(encoding, harm = 1)
#'
#' @seealso \link{plot.fmca} \link{print.fmca} \link{summary.fmca} \link{plotComponent} \link{get_encoding}
#'
#' @author Cristian Preda, Quentin Grimonprez
#' 
#' @references 
#' \itemize{
#'   \item Deville J.C. (1982) Analyse de données chronologiques qualitatives : comment analyser des calendriers ?, Annales de l'INSEE, No 45, p. 45-104.
#'   \item Deville J.C. et  Saporta G. (1980) Analyse harmonique qualitative, DIDAY et al. (editors), Data Analysis and Informatics, North Holland, p. 375-389.
#'   \item Saporta G. (1981) Méthodes exploratoires d'analyse de données temporelles, Cahiers du B.U.R.O, Université Pierre et Marie Curie, 37-38, Paris.
#' }
#' 
#' 
#' @export
compute_optimal_encoding <- function(data, basisobj, computeCI = TRUE, nBootstrap = 50, propBootstrap = 1, nCores = max(1, ceiling(detectCores()/2)), verbose = TRUE, ...)
{
  t1 <- proc.time()
  ## check parameters
  checkData(data)
  checkDataBeginTime(data)
  checkDataEndTmax(data)
  if(!is.basis(basisobj))
    stop("basisobj is not a basis object.")
  if(any(is.na(nCores)) || !is.whole.number(nCores) || (nCores < 1))
    stop("nCores must be an integer > 0.")
  checkLogical(verbose, "verbose")
  checkLogical(computeCI, "computeCI")
  if(computeCI)
  {
    if(any(is.na(nBootstrap)) || (length(nBootstrap) != 1) || !is.whole.number(nBootstrap) || (nBootstrap < 1))
      stop("nBootstrap must be an integer > 0.")
    if(any(is.na(propBootstrap)) || !is.numeric(propBootstrap) || (length(propBootstrap) != 1) || (propBootstrap > 1) || (propBootstrap <= 0))
      stop("propBootstrap must be a real between 0 and 1.")
  }
  ## end check
  
  # used to determine the moments where the probability is 0 in plot.fmca
  pt <- estimate_pt(data)
  
  if(verbose)
    cat("######### Compute encoding #########\n")
  
  
  # change state as integer
  out <- stateToInteger(data$state)
  data$state = out$state
  label <- out$label
  rm(out)
  
  uniqueId <- unique(data$id)
  nId <- length(uniqueId)

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
  

  V <- computeVmatrix(data, basisobj, K, nCores, verbose, ...)
  
  Uval <- computeUmatrix(data, basisobj, K, nCores, verbose, ...)
  
  fullEncoding <- computeEncoding(Uval, V, K, nBasis, uniqueId, label, verbose, manage0 = TRUE)
    
  
  if(computeCI)
  {
    signRef <- getSignReference(fullEncoding$alpha)
      
    bootEncoding <- computeBootStrapEncoding(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, signRef, verbose)
    if(length(bootEncoding) == 0)
    {
      warning("All bootstrap samples return an error. Try to change the basis.")
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt))
    }
    else
    {
      varAlpha <- computeVarianceAlpha(bootEncoding, K, nBasis)
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt, bootstrap = bootEncoding, varAlpha = varAlpha))
    }
    
  }else{
    out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt))
  }
  
  class(out) = "fmca"
  t2 <- proc.time()

  out$runTime = as.numeric((t2-t1)[3])
  if(verbose)
    cat(paste0("Run Time: ", round(out$runTime, 2), "s\n"))
    
  return(out)
}



# return a matrix with nId rows and nBasis * nState columns
computeVmatrix <- function(data, basisobj, K, nCores, verbose, ...)
{
  nBasis <- basisobj$nbasis
  
  phi <- fd(diag(nBasis), basisobj) # les fonctions de base comme données fonctionnelles
  
  # declare parallelization
  if(nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }
  

  if(verbose)
  {
    cat("---- Compute V matrix:\n")
    pbo <- pboptions(char = "=")
  } else {
    pboptions(type = "none")
  }
  t2 <- proc.time()
  
  
  # on construit les variables V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  V <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Vxi, phi = phi, K = K))
  rownames(V) = NULL
  
  t3 <- proc.time()
  
  if(verbose)
  {
    cat(paste0("\nDONE in ", round((t3-t2)[3], 2), "s\n"))
  }
  
  # stop parallelization
  if (nCores > 1)
    stopCluster(cl)
  
  return(V)
}


# compute_Vxi  (plutot Vxj = \int_0^Tphi_j(t)X_t=x dt
#
# @param x one individual (id, time, state) 
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
# @param ... parameters for integrate function
#
# @return vector of size K*nBasis: V[(x=1,i=1)], 
# V[(x=1,i=2)],..., V[(x=1,i=nBasis)], V[(x=2,i=1)], V[(x=2,i=2)], ...
# 
# @examples
# K <- 4
# Tmax <- 6
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
# d_JK2 <- cut_data(d_JK, Tmax)
#
# m <- 6
# b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# I <- diag(rep(1, m))
# phi <- fd(I, b)
# compute_Vxi(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda
compute_Vxi <- function(x, phi, K, ...) 
{
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis)  #V11, V12,...V1m, V21, V22, ..., V2m, ... etc VK1... VKm
  
  for(u in seq_len(nrow(x)-1))
  {
    state = x$state[u]
    
    for(j in seq_len(nBasis))  # j = la base
    {
      aux[(state-1)*nBasis + j] = aux[(state-1)*nBasis + j] +
        integrate(function(t){
          eval.fd(t, phi[j])
        }, lower = x$time[u], upper = x$time[u+1], 
        stop.on.error = FALSE, ...)$value
    }
    
  }

  return(aux)
}

# return a matrix with nId rows and nBasis * nState columns
computeUmatrix <- function(data, basisobj, K, nCores, verbose, ...)
{
  nBasis <- basisobj$nbasis
  
  phi <- fd(diag(nBasis), basisobj) # les fonctions de base comme données fonctionnelles
  
  # declare parallelization
  if(nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }

  t3 <- proc.time()
  if(verbose)
  {
    cat(paste0("---- Compute U matrix:\n"))
    pbo <- pboptions(char = "=")
  } else {
    pbo <- pboptions(type = "none")
  }
  
  Uval <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Uxij, phi = phi, K = K))

  # stop parallelization
  if(nCores > 1)
    stopCluster(cl)
  
  
  t4 <- proc.time()
  if(verbose)
    cat(paste0("\nDONE in ", round((t4-t3)[3], 2), "s\n"))
  
  
  return(Uval)
}



# compute Uxij
#
# compute int_0^T phi_i(t) phi_j(t) 1_{X_{t}=x}  
#
# @param x one individual (id, time, state) 
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
# @param ... parameters for integrate function
#
# @return vector of size K*nBasis*nBasis: U[(x=1,i=1),(x=1,j=1)], 
# U[(x=1,i=1), (x=1,j=2)],..., U[(x=1,i=1), (x=1,j=nBasis)], U[(x=2,i=1), (x=2,j=1)], U[(x=2,i=1), (x=2,j=2)], ...
# 
# @examples
# K <- 4
# Tmax <- 6
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
# d_JK2 <- cut_data(d_JK, Tmax)
#
# m <- 6
# b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# I <- diag(rep(1, m))
# phi <- fd(I, b)
# compute_Uxij(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda, Quentin Grimonprez
compute_Uxij <- function(x, phi, K, ...)
{
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)
  
  for(u in seq_len(nrow(x)-1))
  {
    state <- x$state[u]
    for(i in 1:nBasis) 
    {
      for(j in i:nBasis) # symmetry between i and j  
      {
        
        integral <- integrate(function(t) { 
          eval.fd(t, phi[i]) * eval.fd(t, phi[j])
        }, lower = x$time[u], upper = x$time[u+1], 
        stop.on.error = FALSE, ...)$value
        
        aux[(state-1)*nBasis*nBasis + (i-1)*nBasis + j] = aux[(state-1)*nBasis*nBasis + (i-1)*nBasis + j] + integral
        
        
        # when i == j, we are on the diagonal of the matrix, no symmetry to apply
        if(i != j)
        {
          aux[(state-1)*nBasis*nBasis + (j-1)*nBasis + i] = aux[(state-1)*nBasis*nBasis + (j-1)*nBasis + i] + integral
        }
        
      }
      
    }
    
  }
  
  return(aux)     
}



# @author Cristian Preda, Quentin Grimonprez
computeEncoding <- function(Uval, V, K, nBasis, uniqueId, label, verbose, manage0 = FALSE)
{
  t4 <- proc.time()
  if(verbose)
    cat(paste0("---- Compute encoding: "))
  
  G <- cov(V)
  
  # create F matrix
  Fval = colMeans(Uval)
  Fmat <- matrix(0, ncol = K*nBasis, nrow = K*nBasis) # diagonal-block matrix with K blocks of size nBasis*nBasis
  for(i in seq_len(K))
  {
    Fmat[((i-1)*nBasis+1):(i*nBasis), ((i-1)*nBasis+1):(i*nBasis)] =
      matrix(Fval[((i-1)*nBasis*nBasis+1):(i*nBasis*nBasis)], ncol = nBasis, byrow = TRUE)
  }
  
  # save matrices before modifying them
  outMat <- list(F = Fmat, G = G)
  
  
  # manage the case where there is a column full of 0 (non invertible matrix)
  # if TRUE, we remove the 0-columns (and rows) and throw a warning
  # otherwise we throw an error and the process stops
  if(manage0)
  {
    # column full of 0
    ind0 <- (colSums(Fmat == 0) == nrow(Fmat))
    if(sum(ind0) > 0)
      warning("The F matrix contains at least one column of 0s. At least one state is not present in the support of one basis function. Corresponding coefficients in the alpha output will have a 0 value.")
    
    F05 <- t(mroot(Fmat[!ind0, !ind0]))#F  = t(F05)%*%F05
    G = G[!ind0, !ind0]
    V = V[, !ind0]
    
  }else{
    ind0 <- (colSums(Fmat == 0) == nrow(Fmat))
    
    # res = eigen(solve(F)%*%G)
    F05 <- t(mroot(Fmat)) #F  = t(F05)%*%F05
    
    if(any(dim(F05) != rep(K*nBasis, 2)))
    {
      cat("\n")
      if(any(colSums(Fmat) == 0))
        stop("F matrix is not invertible. In the support of each basis function, each state must be present at least once (p(x_t) != 0 for t in the support). You can try to change the basis.")
      
      stop("F matrix is not invertible. You can try to change the basis.")
    }
  }

  
  invF05 <- solve(F05)
  # res = eigen(F05%*%solve(F)%*%G%*%solve(F05))
  res <- eigen(t(invF05) %*% G %*% invF05)
  
  # les vecteurs propres (qui donnent les coeffs des m=nBasis codages, pour chaque val propre)
  # je les mets sous la forme d'une liste de matrices de taille m x K. La premiere matrice
  # correspond au premier vecteur propre. ce vecteur (1ere colonne dans res$vectors) contient
  # les coefs du codage pour l'état 1 sur les premières m (=nBasis) positions, ensuite pour l'état 2 sur
  # m positions et enfin pour le k-eme état. Je mets cette première colonne sous forme de
  # matrice et les coefs sont sous forme de colonnes de taille m
  
  # met la matrice de vecteurs propres comme une liste
  
  # aux1 = split(res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))
  invF05vec <- invF05 %*% res$vectors
  aux1 <- split(invF05vec, rep(seq_len(ncol(res$vectors)), each = nrow(res$vectors)))
  
  # on construit les matrices m x K pour chaque valeur propre : 1ere colonne les coefs pour etat 1,
  # 2eme col les coefs pour état 2, etc
  
  alpha <- lapply(aux1, function(w) {
    wb = rep(NA, nBasis * K)
    wb[!ind0] = w
    
    return(matrix(wb, ncol = K, dimnames = list(NULL, label$label)))
  })
  
  pc <- V %*% invF05vec
  rownames(pc) = uniqueId
  
  t5 <- proc.time()
  
  if(verbose)
    cat(paste0("\nDONE in ", round((t5-t4)[3], 2), "s\n"))

  return(list(eigenvalues = res$values, alpha = alpha, pc = pc, F = outMat$F, G = outMat$G))
}

