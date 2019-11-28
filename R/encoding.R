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
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, Tmax)
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
    cat("######### Compute encoding #########\n")
  
  
  # change state as integer
  out <- stateToInteger(data_msm$state)
  data_msm$state = out$state
  label <- out$label
  rm(out)
  
  # refactor labels as 1:nbId
  uniqueId <- unique(data_msm$id)
  nId <- length(uniqueId)
  id2 <- refactorCategorical(data_msm$id, uniqueId, seq_along(uniqueId)) 
  uniqueId2 <- unique(id2)
  
  nCores <- min(max(1, nCores), detectCores()-1)
  
  Tmax <- max(data_msm$time)
  K <- length(label$label)
  nBasis <- basisobj$nbasis  # nombre de fonctions de base
  
  if(verbose)
  {
    cat(paste0("Number of individuals: ", nId, "\n"))
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
  V <- foreach(i = uniqueId2, .combine = rbind, .options.snow = opts)%dopar%{
    res <- compute_Vxi(data_msm[id2 == i, ], phi, K, ...)
    if((nCores == 1) && verbose)
      setTxtProgressBar(pb, i)
    return(res)
  }
  rownames(V) = NULL
  G = cov(V)
  t3 <- proc.time()
  
  if(verbose)
  {
    close(pb)
    cat(paste0("\nDONE in ", round((t3-t2)[3], 2), "s\n---- Compute F matrix:\n"))
    pb <- txtProgressBar(0, nId, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{
    opts <- list()
  }

  
  Fval <- foreach(i = uniqueId2, .combine = rbind, .options.snow = opts)%dopar%{
    res <- compute_Uxij(data_msm[id2 == i, ], phi, K, ...)
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
    stop("In the support of each basis function, each state must be present at least once (p(x_t) != 0 for t in the support).")
  
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
# QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_QJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
# d_JK2 <- msm2msmTmax(d_JK, Tmax)
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
# Tmax <- 6
# QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_QJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
# d_JK2 <- msm2msmTmax(d_JK, Tmax)
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
#' @param harm harmonic to use for the encoding
#' @param ... not used
#'
#'
#' @details 
#' The encoding is \eqn{a_{x} \approx \sum_{i=1}^m \alpha_{x,i}\phi_i}.
#'
#' @author Quentin Grimonprez
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, Tmax)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding 
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # plot the encoding using 1 harmonic
#' plot(encoding)
#' 
#' @seealso \link{plotComponent} \link{plotEigenvalues}
#' 
#' @export
plot.fmca <- function(x, harm = 1, ...)
{
  fdmat <- getEncoding(x, harm = harm, fdObject = FALSE)
  df <- data.frame(x = rep(fdmat$x, ncol(fdmat$y)), y = as.vector(fdmat$y), State = factor(rep(colnames(fdmat$y), each = nrow(fdmat$y)), levels = colnames(fdmat$y)))
  
  ggplot(df, aes_string(x = "x", y = "y", group = "State", colour = "State")) +
    geom_line() +
    labs(x = "Time", y = expression(paste("a"["x"], "(t)")), title = paste0("Encoding with harmonic number ", harm))
}


#' Extract the computed encoding
#'
#' Extract the encoding as an \code{fd} object or as a matrix
#' 
#' @param x Output of \code{\link{compute_optimal_encoding}}
#' @param harm harmonic to use for the encoding
#' @param fdObject If TRUE returns a \code{fd} object else a matrix
#' @param nx (Only if \code{fdObject = TRUE}) Number of points to evaluate the encoding
#'
#' @return a \code{fd} object or a list of two elements \code{y}, a matrix with \code{nx} rows containing the encoding of the state and \code{x}, the vector with time values.
#' 
#' @details 
#' The encoding is \eqn{a_{x} \approx \sum_{i=1}^m \alpha_{x,i}\phi_i}.
#' 
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, Tmax)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # extract the encoding using 1 harmonic
#' encodFd <- getEncoding(encoding, fdObject = TRUE)
#' encodMat <- getEncoding(encoding, nx = 200)
#' 
#' @author Cristian Preda
#'
#' @export
getEncoding <- function(x, harm = 1, fdObject = FALSE, nx = NULL)
{
  ## check parameters
  if(class(x) != "fmca")
    stop("x must be a fmca object.")
  checkLogical(fdObject, "fdObject")
  if(!is.null(nx))
  {
    if((length(nx) > 1) || !is.whole.number(nx) || (nx < 0))
      stop("nx must be a positive integer.")
  }
  if((length(harm) > 1) || !is.whole.number(harm) || (harm < 1) || (harm > length(x$alpha)))
    stop("harm must be an integer between 1 and the number of components.")
  ##
  
  fdObj <- fd(x$alpha[[harm]], x$basisobj)

  
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
#' @param nudge_x,nudge_y horizontal and vertical adjustment to nudge labels by
#' @param size size of labels
#'
#' @author Quentin Grimonprez
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, Tmax)
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
#' @seealso \link{plot.fmca} \link{plotEigenvalues}
#' 
#' @export
plotComponent <- function(x, comp = c(1, 2), addNames = TRUE, nudge_x = 0.1, nudge_y = 0.1, size = 4)
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
  df$name = rownames(x$pc)
  
  p <- ggplot(df, aes_string(x = paste0("V", comp[1]), y = paste0("V", comp[2]))) +
    geom_point() + 
    labs(x = paste0("Comp ", comp[1]), y = paste0("Comp ", comp[2]))
  
  if(addNames)
    p = p + geom_text(aes_string(label = "name"), nudge_x = nudge_x, nudge_y = nudge_y, size = size)

  p
}


#' Plot Eigenvalues
#'
#' @param x output of \code{\link{compute_optimal_encoding}} function
#' @param cumulative if TRUE, plot the cumualtive eigenvalues
#' @param normalize if TRUE eigenvalues are normalized for summing to 1
#'
#' @author Quentin Grimonprez
#' 
#' @examples 
#' # simulate the Jukes Cantor models of nucleotides replacement. 
#' K <- 4
#' Tmax <- 6
#' QJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
#' lambda_QJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov_cfd(n = 10, K = K, Q = QJK, lambda = lambda_QJK, Tmax = Tmax)
#' d_JK2 <- msm2msmTmax(d_JK, Tmax)
#'
#' # create basis object
#' m <- 10
#' b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
#' 
#' # compute encoding
#' encoding <- compute_optimal_encoding(d_JK2, b, nCores = 1)
#' 
#' # plot eigenvalues
#' plotEigenvalues(encoding, cumulative = TRUE, normalize = TRUE)
#' 
#' @seealso \link{plot.fmca} \link{plotComponent}
#' 
#' @export
plotEigenvalues <- function(x, cumulative = FALSE, normalize = FALSE)
{
  ## check parameters
  if(class(x) != "fmca")
    stop("x must be a fmca object.")
  checkLogical(cumulative, "cumulative")
  checkLogical(normalize, "normalize")
  ##
  
  if(normalize)
    eigenv <- x$eigenvalues/sum(x$eigenvalues)
  else
    eigenv <- x$eigenvalues
  
  if(cumulative)
    eigenv = cumsum(eigenv)

  
  df <- data.frame(eigenvalues = eigenv, component = seq_along(eigenv))

  p <- ggplot(df, aes_string(x = "component", y = "eigenvalues")) +
    geom_point() + geom_step() +
    labs(title = ifelse(cumulative, "Cumulative eigenvalues", "Eigenvalues"))
  
  p
}
