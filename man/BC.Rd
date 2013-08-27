\name{BC}
\alias{BC}

\title{
Run a BayesComm model
}
\description{
\code{BC} is the main function for running BayesComm models.
It is a wrapper function for \code{BCfit}; it checks inputs, sets up the model types and specifies a number of default \code{BCfit} settings.
}

\usage{
BC(Y, X = NULL, model = "null", covlist = NULL, condition = NULL, its = 100, ...)
}

\arguments{
  \item{Y}{
matrix of species presence/absence data
}
  \item{X}{
matrix of environmental covariates
}
  \item{model}{
type of model to run
}
  \item{covlist}{
optional list of which covariates to assign to each species
}
  \item{condition}{
matrix of conditioning variables
}
  \item{its}{
number of iterations for sampling phase
}
  \item{\dots}{
further arguments to pass to \code{BCfit}
}
}
\details{
\code{Y} must be a matrix with records as rows and species as columns and \code{X} a matrix with records as rows and covariates as columns.
\code{model} must be one of: \code{"null"} (intercept only), \code{"environment"} (intercept & covariates), \code{"community"} (intercept & community matrix) or \code{"full"} (intercept, covariates & community matrix).
\code{covlist} must have the same length as the number of species with, each element a vector of column indices for \code{X}. \code{covlist} defaults to \code{NULL}, which includes all covariates for all species.
For more details of arguments for model fitting see \code{\link{BCfit}}. \code{condition} is an optional matrix of conditioning variables.
These are fitted in the same way as \code{X} but are not removed in null and community models.
}
\value{
An object of class \code{bayescomm} containing the model call and parameter chains which can be viewed and manipulated using \code{window}, \code{plot}, \code{print} and \code{summary}.
}

\seealso{
\code{\link{BCfit}}
}
\examples{
# create fake data
n <- 100
nsp <- 4
k <- 3

X <- matrix(c(rep(1, n), rnorm(n * k)), n)  # covariate matrix
W <- matrix(rnorm(nsp * nsp), nsp)
W <- W \%*\% t(W) / 2  # true covariance matrix
B <- matrix(rnorm(nsp * (k + 1), 0, 3), nsp)  # true covariates
mu <- apply(B, 1, function(b, x) x \%*\% b, X)  # true mean
e <- matrix(rnorm(n * nsp), n) \%*\% chol(W)  # true e
z <- mu + e  # true z
Y <- ifelse(z > 0, 1, 0)  # true presence/absence

# run BC (after removing intercept column from design matrix)
m1 <- BC(Y, X[, -1], model = "full", its = 100)
}
