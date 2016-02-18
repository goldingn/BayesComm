logLik.bayescomm <- function(model, newdata, y, thin = 1){
  if (missing(newdata)) {
    # Already includes intercept column
    x <- model$call$X
  } else {
    # Need to add intercept column
    x <- cbind(1, newdata)
  }
  
  if (missing(y)) {
     y <- model$call$Y
  }
  
  B <- bindSpeciesCoefficients(model)
  
  lower = ifelse(y, 0, -Inf)
  upper = ifelse(y, Inf, 0)
  
  indices = which(1:nrow(B) %% thin == 0)
  
  stopifnot(length(indices) > 0)

  
  pb <- progress::progress_bar$new(
    format = ":current / :total  [:bar] :percent eta: :eta",
    total = nrow(newdata) * length(indices), clear = FALSE, width = 60)
  
  likelihood = sapply(
    indices,
    function(i){
      mean.matrix <- x %*% B[i, , ]
      R <- upper2cor(model$trace$R[i, ])
      
      sapply(
        1:nrow(mean.matrix),
        function(j){
          pb$tick()
          mvtnorm::pmvnorm(
            lower = lower[j, ], 
            upper = upper[j, ], 
            mean = mean.matrix[j, ], 
            corr = R
          )
        }
      )
    }
  )
  log(likelihood)
}

