#' Compute the optimal encodings f each state
#'
#'
#' @param data_msm data.frame containing \code{id}, \code{time} and \code{state} (see \code{\link{generate_Markov_cfd}}). All individual must end at the same time Tmax (use \code{\link{msm2msmTmax}}).
#' @param basisobj is a basis create using \code{fda} package
#'
#' @return a list containing:
#' \itemize{
#'   \item vp eigenvalues
#'   \item alpha optimal encoding coefficients
#'   \item pc principal components
#'   \item F 
#'   \item G  
#'   \item V 
#' }
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
#' encoding <- compute_optimal_encoding(d_JK2, b)
#'
#' @author Cristian Preda
#' 
#' @export
compute_optimal_encoding <- function(data_msm, basisobj)
{
  Tmax <- max(data_msm$time)
  
  K <- length(unique(data_msm$state))
  nBasis <- basisobj$nbasis  # nombre de fonctions de base
  I <- diag(rep(1, nBasis))
  phi <- fd(I, basisobj) #les fonctions de base comme données fonctionnelles
  
  # on construit les variables V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  res <- by(data_msm, data_msm$id, function(x){compute_Vxi(x, phi, K)})
  V <- do.call(rbind, res)
  G = cov(V)
  #print("matrix G : computed !")
  res <- by(data_msm, data_msm$id, function(x){compute_Uxij(x, phi, K)})
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

  return(list(vp = res$values, alpha = aux2, pc = pc, F = Fmat, G = G, V = V))

}



# compute Uxij
#
# compute int_0^T phi_i(t) phi_j(t) 1_{X_{t}=x}  
#
# @param x one individual (id, time, state) 
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
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
# @author Cristian Preda
compute_Uxij <- function(x, phi, K)
{
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)
  
  for(state in 1:K) 
    for(i in 1:nBasis) 
      for(j in 1:nBasis)  
      {
        idx <- which(x$state == state)
        for(u in idx)
        {
          if(u < nrow(x))
          {
            aux[(state-1)*nBasis*nBasis + (i-1)*nBasis + j] = aux[(state-1)*nBasis*nBasis + (i-1)*nBasis + j] + 
              integrate(function(t) {
                eval.fd(t, phi[i]) * eval.fd(t, phi[j])
              }, lower = x[u, "time"], upper = x[u+1, "time"],
              stop.on.error = FALSE)$value
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
compute_Vxi <- function(x, phi, K) 
{
  nBasis <- phi$basis$nbasis
  aux = rep(0, K * nBasis)  #V11, V12,...V1m, V21, V22, ..., V2m, ... etc VK1... VKm
  for(state in 1:K) 
  {
    idx <- which(x$state == state)
    
    for(j in 1:nBasis)  # j = la base
    {
      for(u in idx)
      {
        if(u < nrow(x))
        {
          aux[(state-1)*nBasis + j] = aux[(state-1)*nBasis + j] +
            integrate(function(t){
              eval.fd(t, phi[j])
            }, lower = x[u, "time"], upper = x[u+1, "time"], 
            stop.on.error = FALSE)$value
        }
      }
    }
  }

  return(aux)
}
