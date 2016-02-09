BCfit <-
function(y, X, covlist, R, z, mu, updateR, iters, thin = 1, burn = 0, priW = c(nrow(z) + 2 * ncol(z), 2 * ncol(z)), verbose = 0) {
  
  broken = which(diag(var(y)) == 0)
  
  if (any(diag(var(y)) == 0)) {
    stop("No variation in y for species", paste(broken, collapse = ", "), 
         "; correlation not defined")
  }
  
  stopifnot(all(is.finite(X)), all(is.finite(y)))
    
  spnames <- colnames(y)
  y <- t(y)         # NOTE: y is transposed relative to the input!
  nsp <- dim(y)[1]  # number of species
  n <- dim(y)[2]    # number of sites
  iR <- solve(R)    # Inverse correlation matrix
  e <- z - mu       # component of z that can't be explained by mu
  
  nsamp <- (iters - burn) %/% thin
  
  chain_info <- list(burn = burn, thin = thin)
  
  trace <- list(
    R = array(NA, dim = c(nsamp, ((nsp * nsp - nsp) / 2))),
    B = NULL, 
    z = array(NA, dim = c(nsamp, n, nsp))
  )
  
  for (i in 1:nsp) {
    temp <- matrix(NA, nsamp, length(covlist[[i]]))
    colnames(temp) <- colnames(X)[covlist[[i]]]
    trace$B[[spnames[i]]] <- temp
  }
  rm(temp)
  
  nam <- rep(NA, n * n)
  for (i in 1:nsp) {
    for (j in 1:nsp) {
      nam[(i - 1) * nsp + j] <- paste(spnames[i], "_", spnames[j], sep="")
    }
  }
  
  colnames(trace$R) <- nam[which(upper.tri(diag(nsp)))]
  dimnames(trace$z)[[3]] <- spnames
  
  # start sampler
  rec <- 0 # counter for record number after burn-in and thinning
  for (iter in 1:iters) {
    
    # get truncation values and sample z
    trunc <- find_trunc(mu, y)
    e <- sample_e(e, trunc, iR)
    z <- mu + e
    
    # sample mu and calculate e
    mulis<- sample_mu(z, X, covlist)
    mu <- mulis[[1]]
    e <- z - mu
    
    # sample R
    if (updateR) {
      R <- sample_R(z - mu, priW)
      iR <- chol2inv(chol(R))
    }
    
    if(verbose == 2){
      message(iter)
    }
    if(verbose > 0 & iter == burn){
      message("burn-in complete")
    }
    
    # record parameters
    if (iter %% thin == 0 & iter > burn) {
      if(verbose == 1){
        message(iter)
      }
      rec <- rec + 1
      trace$R[rec, ] <- R[upper.tri(R)]
      trace$z[rec, , ] <- z
      for (i in 1:nsp) {
        trace$B[[i]][rec, ] <- mulis[[2]][[i]] 
      }
    }
  }  # sampler
  
  out = list(
    trace = trace,
    call = list(model = NULL, Y = y, X = X, covlist = covlist, its = iters,
                     start = chain_info$burn + 1, thin = chain_info$thin)
  )
  
  class(out) <- "bayescomm"
  
  out
}
