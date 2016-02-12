simulate.bayescomm <- function(object, nsim = NULL, seed = NULL, 
                               replace = FALSE, ...){
  
  ## Begin code from `simulate.lm`:
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ## End code from `simulate.lm`

  if (is.null(nsim)) {
    nsim <- nrow(object$trace$R)
  }
   
  indices <- sample.int(nrow(object$trace$R), nsim, replace = replace)
   
  if (!is.null(object$other$mu)){
    stop("predictions are not supported for non-null mu")
  }
  
  B <- bindSpeciesCoefficients(object)
  X <- object$call$X
  
  val <- array(NA, dim = c(nsim, dim(object$call$Y)))
  
  for (i in 1:length(indices)) {
    index = indices[i]
    R <- upper2cor(object$trace$R[index, ])
    z <- X %*% B[index, , ] + mvtnorm::rmvnorm(nrow(X), sigma = R)
    val[i, , ] <- ifelse(z > 0, 1, 0)
  }

  
  ## Begin code from `simulate.lm`
  attr(val, "seed") <- RNGstate
  val
  ## End code from `simulate.lm`
}
