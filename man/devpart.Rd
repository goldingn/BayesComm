\name{devpart}
\alias{devpart}

\title{
Deviance partitioning
}
\description{
Runs a deviance partitioning procedure on a set of four \code{bayescomm} objects.
}
\usage{
devpart(null, environment, community, full)
}

\arguments{
  \item{null}{
a \code{bayescomm} object containing a 'null' model
}
  \item{environment}{
a \code{bayescomm} object containing an 'environment' model
}
  \item{community}{
a \code{bayescomm} object containing a 'community' model
}
  \item{full}{
a \code{bayescomm} object containing a 'full' model
}
}
\details{
The deviance partitioning procedure determines the proportion of the null deviance explained by each of the other three model types.
The four model types are those created by \code{BC}.
}
\value{
A list containing elements
	\item{devpart}{matrix containing the proportion of the null deviance explained by each model for each species}
	\item{null}{a matrix containing the mean and 95\% credible intervals for the deviance for each species in the null model}
	\item{environment}{a matrix containing the mean and 95\% credible intervals for the deviance for each species in the evironment model}
	\item{community}{a matrix containing the mean and 95\% credible intervals for the deviance for each species in the community model}
	\item{full}{a matrix containing the mean and 95\% credible intervals for the deviance for each species in the full model}
	}
\seealso{
\code{\link{BC}}
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
null <- BC(Y, X[, -1], model = "null", its = 100)
comm <- BC(Y, X[, -1], model = "community",its = 100)
envi <- BC(Y, X[, -1], model = "environment", its = 100)
full <- BC(Y, X[, -1], model = "full", its = 100)

devpart(null, envi, comm, full)
}


