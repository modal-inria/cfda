#' Compute the optimal encoding for each state
#'
#' Compute the optimal encoding for categorical functional data using an extension of the multiple correspondence analysis to a stochastic process.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state. All individuals must begin at the same time T0 and end at the same time Tmax (use \code{\link{cut_data}}).
#' @param basisobj basis created using the \code{fda} package (cf. \code{\link{create.basis}}).
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
#' Tmax <- 6
#' PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax,
#'                             labels = c("A", "C", "G", "T"))
#' d_JK2 <- cut_data(d_JK, Tmax)
#'
#' # create basis object
#' m <- 6
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
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
#' @seealso \link{plot.fmca} \link{plotComponent} \link{get_encoding}
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
compute_optimal_encoding <- function(data, basisobj, nCores = max(1, ceiling(detectCores()/2)), verbose = TRUE, ...)
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
  
  nCores <- min(max(1, nCores), detectCores()-1)
  
  Tmax <- max(data$time)
  K <- length(label$label)
  nBasis <- basisobj$nbasis  # nombre de fonctions de base
  
  if(verbose)
  {
    cat(paste0("Number of individuals: ", nId, "\n"))
    cat(paste0("Number of states: ", K, "\n"))
    cat(paste0("Basis type: ", basisobj$type, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
    cat(paste0("Number of cores: ", nCores, "\n"))
  }
  

  I <- diag(rep(1, nBasis))
  phi <- fd(I, basisobj) #les fonctions de base comme données fonctionnelles
  
  
  
  V <- computeVmatrix(data, uniqueId2, id2, basisobj, K, nCores, verbose, ...)
  
  G <- cov(V)
  
  
  
  
  # declare parallelization
  if(nCores > 1)
  {
    cl <- makeCluster(nCores)
    registerDoSNOW(cl)
  }else{
    registerDoSEQ()
  }

  t3 <- proc.time()
  if(verbose)
  {
    cat(paste0("---- Compute F matrix:\n"))
    pb <- txtProgressBar(0, nId, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- list()
  }
  
  Fval <- foreach(i = uniqueId2, .combine = rbind, .options.snow = opts)%dopar%{
    res <- compute_Fxij(data[id2 == i, ], phi, K, ...)
    if((nCores == 1) && verbose)
      setTxtProgressBar(pb, i)
    return(res)
  } 
    
  if(verbose)
    close(pb)
  
  # stop parallelization
  if(nCores > 1)
    stopCluster(cl)
  
  # create F matrix
  Fval = colMeans(Fval)
  Fmat <- matrix(0, ncol = K*nBasis, nrow = K*nBasis) #matrice avec K blocs de taille nBasis*nBasis sur la diagonale
  for(i in 1:K)
  {
    Fmat[((i-1)*nBasis+1):(i*nBasis), ((i-1)*nBasis+1):(i*nBasis)] =
      matrix(Fval[((i-1)*nBasis*nBasis+1):(i*nBasis*nBasis)], ncol = nBasis, byrow = TRUE)
  }
  
  t4 <- proc.time()
  if(verbose)
    cat(paste0("\nDONE in ", round((t4-t3)[3], 2), "s\n---- Compute encoding: "))

  #res = eigen(solve(F)%*%G)
  F05 <- t(mroot(Fmat)) #F  = t(F05)%*%F05

  if(any(dim(F05) != rep(K*nBasis, 2)))
  {
    cat("\n")
    if(any(colSums(Fmat) == 0))
      stop("F matrix is not invertible. In the support of each basis function, each state must be present at least once (p(x_t) != 0 for t in the support). You can try to change the basis.")
    
    stop("F matrix is not invertible. You can try to change the basis.")
  }
  
  invF05 <- solve(F05)
  #res = eigen(F05%*%solve(F)%*%G%*%solve(F05))
  res <- eigen(t(invF05) %*% G %*% invF05)

  # les vecteurs propres (qui donnent les coeffs des m=nBasis codages, pour chaque val propre)
  # je les mets sous la forme d'une liste de matrices de taille m x K. La premiere matrice
  # correspond au premier vecteur propre. ce vecteur (1ere colonne dans res$vectors) contient
  # les coefs du codage pour l'état 1 sur les premières m (=nBasis) positions, ensuite pour l'état 2 sur
  # m positions et enfin pour le k-eme état. Je mets cette première colonne sous forme de
  # matrice et les coefs sont sous forme de colonnes de taille m

  # met la matrice de vecteurs propres comme une liste

  # aux1 = split(res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))
  aux1 = split(invF05 %*% res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))

  # on construit les matrices m x K pour chaque valeur propre : 1ere colonne les coefs pour etat 1,
  # 2eme col les coefs pour état 2, etc

  aux2 <- lapply(aux1, function(w){return(matrix(w, ncol = K, dimnames = list(NULL, label$label)))})

  pc <- V %*% (invF05 %*% res$vectors)
  rownames(pc) = uniqueId
  
  if(verbose)
    cat("DONE\n")
  
  out <- list(eigenvalues = res$values, alpha = aux2, pc = pc, F = Fmat, G = G, V = V, basisobj = basisobj)
  class(out) = "fmca"
  t6 <- proc.time()

  if(verbose)
    cat(paste0("Run Time: ", round((t6-t1)[3], 2), "s\n"))
    
  return(out)
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
# compute_Fxij(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda, Quentin Grimonprez
compute_Fxij <- function(x, phi, K, ...)
{
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)
  
  for(u in 1:(nrow(x)-1))
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



computeVmatrix <- function(data, uniqueId, id, basisobj, K, nCores, verbose, ...)
{
  nId <- length(uniqueId)
  nBasis <- basisobj$nbasis
  
  I <- diag(rep(1, nBasis))
  phi <- fd(I, basisobj) #les fonctions de base comme données fonctionnelles
  
  # declare parallelization
  if(nCores > 1)
  {
    cl <- makeCluster(nCores)
    registerDoSNOW(cl)
  }else{
    registerDoSEQ()
  }
  
  if(verbose)
  {
    cat("---- Compute V matrix:\n")
    pb <- txtProgressBar(0, nId, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- list()
  }
  t2 <- proc.time()
  
  
  # on construit les variables V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  V <- foreach(i = uniqueId, .combine = rbind, .options.snow = opts)%dopar%{
    res <- compute_Vxi(data[id == i, ], phi, K, ...)
    if((nCores == 1) && verbose)
      setTxtProgressBar(pb, i)
    return(res)
  }
  rownames(V) = NULL
  t3 <- proc.time()
  
  if(verbose)
  {
    close(pb)
    cat(paste0("\nDONE in ", round((t3-t2)[3], 2), "s\n"))
  }
  
  # stop parallelization
  if(nCores > 1)
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
  
  for(u in 1:(nrow(x)-1))
  {
    state = x$state[u]
    
    for(j in 1:nBasis)  # j = la base
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


