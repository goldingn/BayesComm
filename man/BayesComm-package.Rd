\name{BayesComm-package}
\alias{BayesComm-package}
\alias{BayesComm}
\docType{package}
\title{
Bayesian community ecology analysis
}
\description{
BayesComm fits Bayesian multivariate binary (probit) regression models for analysis of ecological communities.
These models can be used to make inference about underlying inter-species interactions in communities and to separate the effects of environmental covariates and inter-species interactions on community assembly.
This package accompanies the paper (in preparation) by Golding et al. (2013) and is based on a model described by Edwards and Allenby (2003).
}
\details{
\tabular{ll}{
Package: \tab BayesComm\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2014-02-12\cr
License: \tab GPL (>=2)\cr
}

BayesComm models take as input a matrix of species presence/absence records and optionally a matrix of environmental covariates.
\code{\link{BC}} is the main function for setting up models.
It is a wrapper function to \code{\link{BCfit}} and returns a \code{bayescomm} object.
\code{bayescomm} objects have associated \code{summary}, \code{plot}, \code{print}, \code{window} and \code{residuals} functions.
Functions are also provided to calculate Deviance Information Criteria (\code{\link{DIC}}) and run a deviance partitioning procedure on model outputs (\code{\link{devpart}}).
\cr
\cr
Full details of formulation of the model are given in Golding et al. (2013).
}
\author{
Nick Golding <nick.golding@zoo.ox.ac.uk>
\&
Dave Harris
}
\references{
Golding (2013) Mapping and understanding the distributions of potential vector mosquitoes in the UK: New methods and applications. (Chapter 3) \url{http://dx.doi.org/10.6084/m9.figshare.767289}\cr
Edwards, Y., Allenby, G. (2003) Multivariate analysis of multiple response data. \emph{Journal of Marketing Research}, 40 (3) 321-34.
}

\keyword{ package }
\seealso{
\code{\link{BC}}, \code{\link{BCfit}}, 
\code{\link{window.bayescomm}}, \code{\link{plot.bayescomm}}, \code{\link{print.bayescomm}}, \code{\link{summary.bayescomm}}, \code{\link{residuals.bayescomm}}, 
\code{\link{DIC}}, \code{\link{devpart}},
}