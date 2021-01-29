# Author: Quentin Grimonprez

context("Manage column fills with 0 in compute_optimal_encoding")


oldcompute_optimal_encoding <- function(data, basisobj, nCores = max(1, ceiling(detectCores()/2)), verbose = TRUE, ...)
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
  
  # declare parallelization
  cl <- makeCluster(nCores)
  parallel::clusterExport(cl, "eval.fd")
  if(verbose)
  {
    cat("---- Compute V matrix:\n")
  }
  t2 <- proc.time()
  
  
  # on construit les variables V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  V <- do.call(rbind, pblapply(cl=cl, split(data, data$id), cfda:::compute_Vxi, phi = phi, K = K))
  rownames(V) = NULL
  G = cov(V)
  t3 <- proc.time()
  
  if(verbose)
  {
    cat(paste0("\nDONE in ", round((t3-t2)[3], 2), "s\n---- Compute F matrix:\n"))
  }
  
  Fval <- do.call(rbind, pblapply(cl = cl, split(data, data$id), cfda:::compute_Uxij, phi = phi, K = K))
  
  # stop parallelization
  stopCluster(cl)
  
  # create F matrix
  Fval = colMeans(Fval)
  Fmat <- matrix(0, ncol = K*nBasis, nrow = K*nBasis) # matrice avec K blocs de taille nBasis*nBasis sur la diagonale
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


test_that("compute_optimal_encodings has the same results as before when there is no 0-column", {
  skip_on_cran()
  
  set.seed(42)
  n <- 25
  Tmax <- 1
  K <- 2
  m <- 10
  d <- generate_2State(n)
  dT <- cut_data(d, Tmax)
  row.names(dT) = NULL
  
  b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
  expect_silent(fmcaNew <- compute_optimal_encoding(dT, b, computeCI = FALSE, nCores = 1, verbose = FALSE))
  
  expect_silent(fmcaOld <- oldcompute_optimal_encoding(dT, b, nCores = 1, verbose = FALSE))
  
  expect_equal(fmcaOld$eigenvalues, fmcaNew$eigenvalues)
  expect_equal(fmcaOld$alpha, fmcaNew$alpha)
  expect_equal(fmcaOld$pc, fmcaNew$pc)
  expect_equal(fmcaOld$F, fmcaNew$F)
  expect_equal(fmcaOld$G, fmcaNew$G)
  expect_equal(fmcaOld$V, fmcaNew$V)
  expect_false(any(is.na(fmcaNew$alpha[[1]])))
})


test_that("compute_optimal_encodings works when there is some 0-column", {
  skip_on_cran()
  
  data(biofam2)
  d <- biofam2[biofam2$id <= 100, ]
  nState <- length(unique(d$state))
  nbasis <- 8 
  b <- create.bspline.basis(c(15, 30), nbasis = nbasis, norder = 4)
  
  expect_warning({fmcaNew <- compute_optimal_encoding(d, b, computeCI = FALSE, nCores = 1, verbose = FALSE)}, 
                 regexp = "The F matrix contains at least one column of 0s. At least one state is not present in the support of one basis function. Corresponding coefficients in the alpha output will have a 0 value.",
                 fixed = TRUE)
  expect_error({fmcaOld <- oldcompute_optimal_encoding(d, b, nCores = 1, verbose = FALSE)}, 
               regexp = "F matrix is not invertible. In the support of each basis function, each state must be present at least once (p(x_t) != 0 for t in the support). You can try to change the basis.",
               fixed = TRUE)
  
  expect_equal(dim(fmcaNew$F), c(nbasis * nState, nbasis * nState))
  expect_equal(dim(fmcaNew$G), c(nbasis * nState, nbasis * nState))
  expect_equal(dim(fmcaNew$V), c(length(unique(d$id)), nbasis * nState))
  expect_true(any(is.na(fmcaNew$alpha[[1]])))
  
  out <- get_encoding(fmcaNew, fdObject = FALSE)
  expect_true(any(is.na(out$y)))
  
  expect_warning({fmcaNewBootstrap <- compute_optimal_encoding(d, b, computeCI = TRUE, nCores = 1, verbose = FALSE)}, 
                 regexp = "The F matrix contains at least one column of 0s. At least one state is not present in the support of one basis function. Corresponding coefficients in the alpha output will have a 0 value.",
                 fixed = TRUE)
  
  expect_true("bootstrap" %in% names(fmcaNewBootstrap))
  expect_true("varAlpha" %in% names(fmcaNewBootstrap))
})


# test_that("oldcompute_optimal_encoding throws an error when the basis is not well suited", {
# 
#   data_msm <- data.frame(id = rep(1:2, each = 3), time = c(0, 3, 5, 0, 4, 5), state = c(1, 2, 2, 1, 2, 2))
#   b <- create.bspline.basis(c(0, 5), nbasis = 3, norder = 2)
# 
#   expect_error({fmca <- oldcompute_optimal_encoding(data_msm, b, nCores = 1)},
#                regexp = "F matrix is not invertible. In the support of each basis function, each state must be present at least once (p(x_t) != 0 for t in the support). You can try to change the basis.",
#                fixed = TRUE)
# })

test_that("removeTimeAssociatedWithNACoeff works", {
  fdmat <- matrix(rnorm(80), nrow = 20, ncol = 4, dimnames = list(NULL, c("b", "d", "a", "e")))
  timeVal <- 1:20
  
  pt <- list(t = seq(1, 20, length = 10), pt = matrix(runif(50, 0, 1), nrow = 5, ncol = 10, dimnames = list(c(letters[1:5]), NULL)))
  class(pt) = "pt"
  pt$pt[1, 1:5] = 0
  pt$pt[4, 3:5] = 0
  pt$pt[4, 9:10] = 0
  pt$pt[5, 8:10] = 0
  
  out <- removeTimeAssociatedWithNACoeff(fdmat, timeVal, pt)
  expect_equal(dim(out), dim(fdmat))
  expect_equal(colnames(out), colnames(fdmat))
  expect_equal(out[!is.na(out)], fdmat[!is.na(out)])
  expect_true(all(is.na(out[1:11, 3])))
  expect_true(all(is.na(out[6:11, 2])))
  expect_true(all(is.na(out[18:20, 2])))
  expect_true(all(is.na(out[16:20, 4])))
})



test_that("predict works when there is some 0-column", {
  skip_on_cran()
  
  data(biofam2)
  d <- biofam2[biofam2$id <= 100, ]
  nState <- length(unique(d$state))
  nbasis <- 8 
  b <- create.bspline.basis(c(15, 30), nbasis = nbasis, norder = 4)

  fmcaNewBootstrap <- compute_optimal_encoding(d, b, computeCI = TRUE, nCores = 1, verbose = FALSE)
  
  out <- predict(fmcaNewBootstrap, d, nCores = 1, verbose = FALSE)

  expect_equivalent(out, fmcaNewBootstrap$pc)
})

