# Helper function to bind the coefficient lists in B into an array.
# Currently assumes that (i.e. no `covlist` was specified)
bindSpeciesCoefficients <- function(object) {
  pre.binding <- lapply(object$trace$B, function(x) {
    dim(x) <- c(dim(x), 1)
    x
  })
  
  out <- do.call(abind::abind, pre.binding)
  colnames(out) <- colnames(object$trace$B[[1]])
  
  out
} 
