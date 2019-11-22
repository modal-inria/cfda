#' Compute the optimal encodings for each state
#'
#' Compute the optimal encodings for categorical functional data using an extension of the multiple correspondence analysis to a stochastic process.
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state. All individual must end at the same time Tmax (use \code{\link{msm2msmTmax}}).
#' @param basisobj basis created using the \code{fda} package.
#' @param nCores number of cores used for parallelization. Default is the half of cores.
#' @param verbose if TRUE print some information
#' @param ... parameters for \code{\link{integrate}} function.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{eigenvalues} eigenvalues
#'   \item \code{alpha} optimal encoding coefficients
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
#' @examples
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 10
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, 10)
#'
#' # create basis object
#' m <- 10
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
#' getEncoding(encoding)
#'
#' @seealso \link{plot.fmca} \link{plotComponent} \link{getEncoding}
#'
#' @author Cristian Preda, Quentin Grimonprez
#' 
#' @references 
#' Deville J.C. (1982) Analyse de données chronologiques qualitatives : comment analyser des calendriers ?, Annales de l'INSEE, No 45, p. 45-104.
#' 
#' @export
compute_optimal_encoding <- function(data_msm, basisobj, nCores = max(1, ceiling(detectCores()/2)), verbose = TRUE,  ...)
{
  t1 <- proc.time()
  ## check parameters
  checkDataMsm(data_msm)
  checkDataEndTmax(data_msm)
  if(!is.basis(basisobj))
    stop("basisobj is not a basis object.")
  if(!is.whole.number(nCores) || (nCores < 1))
    stop("nCores must be an integer > 0.")
  ## end check
  
  if(verbose)
    cat("######### Compute encodings #########\n")
  
  
  # change state as integer
  out <- stateToInteger(data_msm$state)
  data_msm$state = out$state
  label <- out$label
  rm(out)
  
  nCores <- min(max(1, nCores), detectCores()-1)
  
  Tmax <- max(data_msm$time)
  K <- length(label$label)
  nBasis <- basisobj$nbasis  # nombre de fonctions de base
  
  if(verbose)
  {
    cat(paste0("Number of individuals: ", length(unique(data_msm$id)), "\n"))
    cat(paste0("Number of states: ", K, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
    cat(paste0("Number of cores: ", nCores, "\n"))
  }
  

  I <- diag(rep(1, nBasis))
  phi <- fd(I, basisobj) #les fonctions de base comme données fonctionnelles
  
  # declare parallelization
  if(nCores > 1)
  {
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
  }else{
    registerDoSEQ()
  }


  if(verbose)
    cat("---- Compute V matrix:\n")
  t2 <- proc.time()
  
  # on construit les variables V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  V <- foreach(i = unique(data_msm$id), .combine = rbind)%dopar%{
    return(compute_Vxi(data_msm[data_msm$id == i, ], phi, K, ...))
  }
  rownames(V) = NULL
  G = cov(V)
  t3 <- proc.time()
  
  if(verbose)
    cat(paste0("DONE in ", round((t3-t2)[3], 2), "s\n---- Compute F matrix:\n"))
  
  Fval <- foreach(i = unique(data_msm$id), .combine = rbind)%dopar%{
    return(compute_Uxij(data_msm[data_msm$id == i, ], phi, K, ...))
  } 
    
  # stop parallelization
  if(nCores > 1)
  {
    stopCluster(cl)
  }
  
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
    cat(paste0("DONE in ", round((t4-t3)[3], 2), "s\n---- Compute encoding: "))
  
  #res = eigen(solve(F)%*%G)
  F05 <- t(mroot(Fmat)) #F  = t(F05)%*%F05
  invF05 <- solve(F05)
  #res = eigen(F05%*%solve(F)%*%G%*%solve(F05))
  res <- eigen(t(invF05) %*% G %*% invF05)

  # les vecteurs propres (qui donneent les coeffs des m=nBasis codages, pour chaque val propre)
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
# Tmax <- 10
# QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_QJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
# d_JK2 <- msm2msmTmax(d_JK, 10)
#
# m <- 10
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
# Tmax <- 10
# QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_QJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
# d_JK2 <- msm2msmTmax(d_JK, 10)
#
# m <- 10
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


#' Plot the optimal encoding
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param ... not used
#'
#' @author Quentin Grimonprez
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 10
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, 10)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # plot the optimal encoding
#' plot(encoding)
#' 
#' @seealso \link{plotComponent}
#' 
#' @export
plot.fmca <- function(x, ...)
{
  fdmat <- getEncoding(x, fdObject = FALSE)
  df <- data.frame(x = rep(fdmat$x, ncol(fdmat$y)), y = as.vector(fdmat$y), State = factor(rep(colnames(fdmat$y), each = nrow(fdmat$y)), levels = colnames(fdmat$y)))
  
  ggplot(df, aes_string(x = "x", y = "y", group = "State", colour = "State")) +
    geom_line() +
    labs(x = "Time", y = expression(paste("a"["x"], "(t)")), title = "First eigen optimal encoding")
}


#' Extract the computed encoding
#'
#' Extract the optimal encoding as an \code{fd} object or as a matrix
#' 
#' @param x Output of \code{\link{compute_optimal_encoding}}
#' @param fdObject If TRUE returns a \code{fd} object else a matrix
#' @param nx (Only if \code{fdObject = TRUE}) Number of points to evaluate the encoding
#'
#' @return a \code{fd} object or a list of two elements \code{y}, a matrix with \code{nx} rows containing the encoding of the state and \code{x}, the vector with time values.
#' 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 10
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, 10)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # extract the encoding
#' encodFd <- getEncoding(encoding, fdObject = TRUE)
#' encodMat <- getEncoding(encoding, nx = 200)
#' 
#' @author Cristian Preda
#'
#' @export
getEncoding <- function(x, fdObject = FALSE, nx = NULL)
{
  fdObj <- fd(x$alpha[[1]], x$basisobj)
  
  if(fdObject)
  {
    return(fdObj)
  }else{
    rangex <- fdObj$basis$rangeval
    nBasis <- fdObj$basis$nbasis
    
    if(is.null(nx))
      nx = max(501, 10 * nBasis + 1)
    
    timeVal <- seq(rangex[1], rangex[2], length = nx)
    
    fdmat <- eval.fd(timeVal, fdObj)
    
    return(list(x = timeVal, y = fdmat))
  }

}



#' Plot Components
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param comp a vector of two elements indicating the components to plot
#' @param addNames if TRUE, add the id labels on the plot
#'
#' @author Quentin Grimonprez
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 10
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, 10)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' plotComponent(encoding, comp = c(1, 2))
#' 
#' @seealso \link{plot.fmca}
#' 
#' @export
plotComponent <- function(x, comp = c(1, 2), addNames = TRUE)
{
  ## check parameters
  if(class(x) != "fmca")
    stop("x must be a fmca object.")
  checkLogical(addNames, "addNames")
  if(length(comp) != 2)
    stop("comp must be a vector of positive integers of length 2.")
  if(any(!is.whole.number(comp)) || any(comp < 0))
    stop("comp must be a vector of positive integers of length 2.")
  ##
  
  df <- as.data.frame(x$pc)
  df$name = 1:nrow(x$pc)
  
  p <- ggplot(df, aes_string(x = paste0("V", comp[1]), y = paste0("V", comp[2]))) +
    geom_point() + 
    labs(x = paste0("Comp ", comp[1]), y = paste0("Comp ", comp[2]))
  
  if(addNames)
    p = p + geom_text(aes_string(label = "name"), hjust = -0.15, vjust = -0.15)

  
  p
}
