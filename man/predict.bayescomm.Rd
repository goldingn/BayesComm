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
=======
\alias{predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make predictions at new locations
}
\description{
For each set of parameter values sampled by the model (including values of Z),
simulate the occurrence probabilities for each species at each new location.
}
\usage{
predict(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A bayescomm object
}
  \item{newdata}{
A data.frame with the same columns as X from the original BC model
}
  \item{...}{
further arguments passed to or from other methods.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An array of occurrence probabilities.
Rows index locations.  Columns index species.  Slices index MCMC samples.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David J. Harris (http://davharris.github.io)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
>>>>>>> 21d58ba0a786d6ef57b86427119745bbdf534795
