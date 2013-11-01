predict.bayescomm <- function(model, newData, species = NULL,
                              niche = c('realised', 'fundamental'),
                              trace = FALSE) {
  
  if (is.null(species)) {
    # if they want all species, do them individually by recursion
    
    # number of species to do
    k <- length(model$trace$B)
    
    # loop through the species
    pred_list <- lapply(1:k,
                        function(sp, model, newData, niche, trace) {
                          predict.bayescomm(model, newData, sp, niche, trace)
                        }, model, newData, niche, trace)
    
    #     combine the columns
    #     pred <- do.call(cbind, pred_list)
    
    # add their names back in
    names(pred_list) <- colnames(model$call$Y)
    
    # and return this list (ignores rest of this function)
    return(pred_list)
  }
  
  # otherwise just do the one they asked for
  
  # match the niche argument
  niche <- match.arg(niche)
  
  # if the first column of newData isn't all ones
  if (!all(newData[, 1] == rep(1, nrow(newData)))) {

    # add a new column which is
    newData <- cbind(intercept = rep(1, nrow(newData)),
                     newData)
  }
  
  # now check it's the right number of covariates
  if (ncol(newData) != ncol(model$call$X)) {
    stop('number of columns in newData does not match number of columns used to fit model')
  }
  
  #  predict probability of presence of species i in the absence of competition
  # at environmental conditions given by newData
  mu <- model$trace$B[[species]] %*% t(newData[, model$call$covlist[[species]]])
  
  # if they wanted the realised niche (in presence of competition) add on
  # the correlated error terms
  if (niche == 'realised') {
    # if realised
    mu <- mu + model$trace$z[, , species]
  }
  
  # transpose and switch these back to probability scale 
  # if trace = TRUE rows are now records, columns are iterations
  pred <- t(mu)
  
  if (!trace) {
    # if they don't want the full trace posterior, just give them the MAP estimatez
    pred <- rowMeans(pred)
  }
  
  # return the prediction
  return(pnorm(pred)) 
}