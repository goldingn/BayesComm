predict.bayescomm <- function(object, newdata, ...) {
  # Return a matrix of marginal probabilitiesat new sites.  
  
  # Probably only works with the "full" model type and no `covlist`, `mu`, or
  # `condition` specified.
  
  if(!is.null(object$other$mu)){
    stop("predictions are not supported for non-null mu")
  }
  
  X <- cbind(intercept = 1, newdata)
  B <- bindSpeciesCoefficients(object)
  
  predictions <- array(
    NA, 
    c(nrow(X), dim(B)[3], dim(B)[1]),
    dimnames = list(row.names(X), dimnames(B)[[3]], NULL)
  )
  
  
  for (i in 1:nrow(B)) {
    predictions[ , , i] <- pnorm(X %*% B[i, , ])
  }
  
  predictions
}
