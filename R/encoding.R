#' Compute the optimal encodings for each state
#'
#'
#' @param data_msm data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and \code{state}, associated state (integer starting at 1). All individual must end at the same time Tmax (use \code{\link{msm2msmTmax}}).
#' @param basisobj is a basis create using \code{fda} package.
#' @param nCores number of Cores used for parallelization. Default is the half of cores.
#' @param ... parameters for \code{\link{integrate}} function.
#'
#' @return A list containing:
#' \itemize{
#'   \item eigenvalues eigenvalues
#'   \item alpha optimal encoding coefficients
#'   \item pc principal components
#'   \item F matrix containing the \eqn{F_{(x,i)(y,j)}}
#'   \item V matrix containing the \eqn{V_{(x,i)}}
#'   \item G covariance matrix of \code{V}
#'   \item basisobj \code{basisobj} parameter
#' }
#'
#' @details 
#' see the vignette for the mathematical background: \code{vignette("cfda", package = "cfda")}
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
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 2)
#' 
#' # plot the optimal encoding
#' plot(encoding)
#' 
#' # plot the two first components
#' plotComponent(encoding, comp = c(1, 2))
#'
#' @seealso \link{plot.fmca} \link{plotComponent}
#'
#' @author Cristian Preda
#' 
#' @references 
#' Deville J.C. (1982) Analyse de donn\'ees chronologiques qualitatives: comment analyser des calendriers ?, Annales de l'INSEE, No 45, p. 45-104.
#' 
#' @export
compute_optimal_encoding <- function(data_msm, basisobj, nCores = max(1, ceiling(detectCores()/2)), ...)
{

  ## check parameters
  checkDataMsm(data_msm)
  checkDataEndTmax(data_msm)
  if(!is.basis(basisobj))
    stop("basisobj is not a basis object.")
  if(!is.whole.number(nCores) || (nCores < 1))
    stop("nCores must be an integer > 0.")
  ## end check
  
  nCore <- min(max(1, nCores), detectCores())
  
  Tmax <- max(data_msm$time)
  
  K <- length(unique(data_msm$state))
  nBasis <- basisobj$nbasis  # nombre de fonctions de base
  I <- diag(rep(1, nBasis))
  phi <- fd(I, basisobj) #les fonctions de base comme données fonctionnelles
  
  # declare parallelization
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  
  # on construit les variables V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  res <- foreach(i = unique(data_msm$id))%dopar%{
    return(compute_Vxi(data_msm[data_msm$id == i, ], phi, K, ...))
  }
  V <- do.call(rbind, res)
  G = cov(V)
  #print("matrix G : computed !")
  res <- foreach(i = unique(data_msm$id))%dopar%{
    return(compute_Uxij(data_msm[data_msm$id == i, ], phi, K, ...))
  } 
    
  # stop parallelization
  stopCluster(cl)
  
  Fval <- do.call(rbind, res)
  Fval = colMeans(Fval)
  Fmat <- matrix(0, ncol = K*nBasis, nrow = K*nBasis) #matrice avec K bloks de taille nBasis*nBasis sur la diagonale
  for(i in 1:K)
  {
    Fmat[((i-1)*nBasis+1):(i*nBasis), ((i-1)*nBasis+1):(i*nBasis)] =
      matrix(Fval[((i-1)*nBasis*nBasis+1):(i*nBasis*nBasis)], ncol = nBasis, byrow = TRUE)
  }
  #print("matrix F : computed !")

  #res = eigen(solve(F)%*%G)
  F05 <- t(mroot(Fmat)) #F  = t(F05)%*%F05
  #res = eigen(F05%*%solve(F)%*%G%*%solve(F05))
  res <- eigen(solve(t(F05)) %*% G %*% solve(F05))

  # les vecteurs propres (qui donneent les coeffs des m=nBasis codages, pour chaque val propre)
  # je les mets sous la forme d'une liste de matrices de taille m x K. La premiere matrice
  # correspond au premier vecteur propre. ce vecteur (1ere colonne dans res$vectors) contient
  # les coefs du codage pour l'état 1 sur les premières m positions, ensuite pour l'état 2 sur
  # m positions et enfin pour le k-eme état. Je mets cette première colonne sous forme de
  # matrice et les coefs sont sous forme de colonnes de taille m

  # met la matrice de vecteurs propres comme une liste

  # aux1 = split(res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))
  aux1 = split(solve(F05) %*% res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))

  # on construit les matrices m x K pour chaque valeur propre : 1ere colonne les coefs pour etat 1,
  # 2eme col les coefs pour état 2, etc

  aux2 <- lapply(aux1, function(w){return(matrix(w, ncol = K))})

  pc <- V %*% (solve(F05) %*% res$vectors)

  out <- list(eigenvalues = res$values, alpha = aux2, pc = pc, F = Fmat, G = G, V = V, basisobj = basisobj)
  class(out) = "fmca"
  
  
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
    state <- x[u, "state"]
    for(i in 1:nBasis) 
    {
      for(j in i:nBasis) # symmetry between i and j  
      {
        
        integral <- integrate(function(t) { 
          eval.fd(t, phi[i]) * eval.fd(t, phi[j])
        }, lower = x[u, "time"], upper = x[u+1, "time"], 
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
    state = x[u, "state"]
    
    for(j in 1:nBasis)  # j = la base
    {
      aux[(state-1)*nBasis + j] = aux[(state-1)*nBasis + j] +
        integrate(function(t){
          eval.fd(t, phi[j])
        }, lower = x[u, "time"], upper = x[u+1, "time"], 
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
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 2)
#' 
#' plot(encoding)
#' 
#' @seealso \link{plotComponent}
#' 
#' @export
plot.fmca <- function(x, ...)
{
  fdObj <- fd(x$alpha[[1]], x$basisobj)
  rangex <- fdObj$basis$rangeval
  nBasis <- fdObj$basis$nbasis
  nx <- max(c(501, 10 * nBasis + 1))
  y <- seq(rangex[1], rangex[2], length = nx)
  
  fdmat <- eval.fd(y, fdObj)
  df <- data.frame(x = rep(y, ncol(fdmat)), y = as.vector(fdmat), State = factor(rep(1:ncol(fdmat), each = nx), levels = 1:ncol(fdmat)))
  
  ggplot(df, aes_string(x = "x", y = "y", group = "State", colour = "State")) +
    geom_line() +
    labs(x = "Time", y = "a_x(t)", title = "First eigen optimal encoding")
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
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 2)
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
