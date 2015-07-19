\name{predict.bayescomm}
\alias{predict.bayescomm}
\title{
Function to make predictions at new locations
}
\description{
For each set of parameter values sampled by the model (including values of Z),
simulate the occurrence probabilities for each species at each new location.
}
\usage{
\method{predict}{bayescomm}(object, newdata, ...)
}
\arguments{
  \item{object}{
A bayescomm object
}
  \item{newdata}{
A data.frame with the same columns as X from the original BC model
}
  \item{...}{
Further arguments passed to or from other methods.
}
}
\value{
An array of occurrence probabilities.
Rows index locations.  Columns index species.  Slices index MCMC samples.
}
\author{
David J. Harris (http://davharris.github.io)
}
\examples{
# load model from first example
m1 <- example(BC)[[1]]

# use the first five sites of the training data as newdata
newdata <- X[1:5, -1]

# get predictions
prob <- predict(m1, newdata)
}
