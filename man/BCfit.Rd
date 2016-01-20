\name{BCfit}
\alias{BCfit}

\title{
Fit a BayesComm model
}
\description{
\code{BCfit} is the workhorse function for the BayesComm model.
It is highly recommended to use the wrapper function \code{\link{BC}} which checks inputs and sets up different model types and initial values.
\code{BCfit} arguments can be accessed through \code{BC} using the \code{...} argument.
}
\usage{
BCfit(y, X, covlist, R, z, mu, updateR, iters, thin = 1, burn = 0,
	priW = c(nrow(z) + 2 * ncol(z), 2 * ncol(z)), updateMu = TRUE, 
	verbose = 0)
}

\arguments{
  \item{y}{
matrix of species presence/absence data
}
  \item{X}{
matrix of environmental covariates
}
  \item{covlist}{
optional list of which covariates to assign to each species
}
  \item{R}{
initial values for correlation matrix
}
  \item{z}{
initial values for z
}
  \item{mu}{
initial values for mu
}
  \item{updateR}{
logical; if true the correlation matrix is updated, if false it is fixed at \code{R}
}
  \item{iters}{
total number of iterations
}
  \item{thin}{
amount to thin the posterior chains. Defaults to 1 (no thinning)
}
  \item{burn}{
number of iterations to discard at the beginning of the chain
}
  \item{priW}{
prior specification for correlation matrix W
}
  \item{updateMu}{
logical; if true, the mean matrix is updated; if false, it is fixed
at \code{mu} (and \code{B} is not estimated). Defaults to \code{TRUE}. Note that 
this argument cannot be passed from \link{\code{BC}}; it can only used in 
\code{BCfit} directly
}
  \item{verbose}{
how often to print updates to the console.  \code{0} for no updates. \code{1} for
updates every \code{thin} iterations after burnin.  \code{2} for updates every iteration.
}
}
\details{
\code{priW} specifies the inverse Wishart prior on the unknown and unidentifiable covariance matrix W from which the correlation matrix R is derived.
\code{priW} is a vector of length two, the first element specifies the degrees of freedom, the second element is multiplied by an identity matrix to form the scale matrix.
The default for \code{priW} is c(n + 2p, 2p), where n is the number of records and p is the number of species in the community; this therefore forms the prior: iW(n + 2p, 2pI).
This prior was determined to exert minimal influence on the posterior of R whilst limiting dependence of R on the unidentifiable variance parameters of W.\cr
\cr
For further details on how to specify \code{Y}, \code{X} and \code{covlist} see \code{\link{BC}}.
}
\value{
A list containing elements:
	\item{R}{samples from posteriors of the correlation matrix}
	\item{B}{samples from posteriors of regression coefficients (a list of matrices)}
	\item{z}{samples from posteriors of latent variables z}
}
\seealso{
\code{\link{BC}}
}

