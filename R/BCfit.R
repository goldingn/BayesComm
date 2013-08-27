BCfit <-
function(y, X, covlist, R, z, mu, updateR, iters, thin = 1, burn = 0, priW = c(nrow(z) + 2 * ncol(z), 2 * ncol(z))) {
  
  spnames <- colnames(y)
  y <- t(y)
  nsp <- dim(y)[1]
  n <- dim(y)[2]
  iR <- solve(R)
  e <- z - mu
  
  # set up trace arrays
  nsamp <- (iters - burn) %/% thin
  trace_R <- array(NA, dim = c(nsamp, ((nsp * nsp - nsp) / 2))) #w /out redundant columns
  trace_z <- array(NA, dim = c(nsamp, n, nsp))
  trace_mu <- trace_z
  trace_B <- NULL
  for (i in 1:nsp) {
    temp <- matrix(NA, nsamp, length(covlist[[i]]))
    colnames(temp) <- colnames(X)[covlist[[i]]]
    trace_B[[spnames[i]]] <- temp
  }
  rm(temp)
  
  nam <- rep(NA, n * n)
  for (i in 1:nsp) {
    for (j in 1:nsp) {
      nam[(i - 1) * nsp + j] <- paste(spnames[i], "_", spnames[j], sep="")
    }
  }
  
  colnames(trace_R) <- nam[which(upper.tri(diag(nsp)))]
  dimnames(trace_z)[[3]] <- spnames
  dimnames(trace_mu)[[3]] <- spnames
  
  # start sampler
  start <- Sys.time()
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
    
    # record parameters    
    if (iter %% thin == 0 & iter > burn) {
      rec <- (iter - burn) %/% thin
      trace_R[rec, ] <- R[upper.tri(R)]
      trace_z[rec, , ] <- z
      trace_mu[rec, , ] <- mu
      for (i in 1:nsp) {
        trace_B[[i]][rec, ] <- mulis[[2]][[i]] 
      }
    }
  }  # sampler
  
  list(R = trace_R, B = trace_B, z = trace_z, burn = burn, thin = thin)
}