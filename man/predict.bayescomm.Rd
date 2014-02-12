\name{predict.bayescomm}
\alias{predict.bayescomm}

\title{
  Predict from a fitted bayescomm model to a new dataset.
}
\description{
  \code{predict.bayescomm} predicts the posterior probability of presence of each species at new sites, given a matrix \code{newData} giving the values of covariates at those sites.
}
\usage{
  \method{predict}{bayescomm}(object, newdata = NULL, species = NULL,
                              type = c('link', 'response'),
                              niche = c('realised', 'fundamental'),
                              trace = FALSE)
}

\arguments{
  \item{object}{
    a \code{bayescomm} object
  }
  \item{newdata}{
    optionally, a matrix giving the values of covariates at sites at which to predict (rows are sites). If omitted the covariates used to fit \code{object} will be reused.
  }
  \item{species}{
    optionally, a vector giving indices for the species for which the user would like to predict. If \code{NULL} predictions will be made for all species used to fit \code{object}
  }
  \item{type}{
  the type of prediction required. The default of \code{'link'} is on the scale of the linear predictors, the alternative \code{'response'} is on the scale of the response variable - i.e. probability of presence.
  }
  \item{niche}{
  whether to make predictions in the presence (\code{'realised'}) or the absence (\code{'fundamental'}) of interactions with other species.
  }
  \item{trace}{
  whether to return all of the posterior samples for the predictions. If \code{FALSE} only the mean of these samples (the posterior mean prediction) is returned.  
  }
  
  
}
\value{
 If predictions for only one species are requested (i.e. \code{species} is of length one) the function returns either a vector (if \code{trace = FALSE}) or a matrix (if \code{trace = TRUE}; with rows giving records to predict to and columns giving samples from the posterior).
}

\seealso{
  \code{\link{BC}}, \code{\link{summary.mcmc}}
}
\examples{
  # TO DO
}
