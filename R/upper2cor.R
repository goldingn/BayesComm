upper2cor <- function(x){
  # Turn a numeric vector with upper triangle values and
  # complete the correlation matrix with this information
  
  # Dimension of the correlation matrix
  # dim defined by solving length(x) == (dim * (dim - 1)) / 2
  dim = (sqrt(8 * length(x) + 1) + 1) / 2
  
  out = matrix(0, dim, dim)  # Initialize
  out[upper.tri(out)] = x    # Fill in upper triangle
  out = out + t(out)         # Fill in lower triangle
  diag(out) = 1              # Diagonal of correlation matrix is 1
  
  out
}  
