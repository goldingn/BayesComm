\name{summary.bayescomm}
\alias{summary.bayescomm}

\title{
Summarise bayescomm parameter chains
}
\description{
\code{summary.bayescomm} creates summaries of a subset of the parameter chains in a \code{bayescomm} object using the \code{coda} package.
}
\usage{
\method{summary}{bayescomm}(object, chain, ...)
}

\arguments{
  \item{object}{
a \code{bayescomm} object
}
  \item{chain}{
a character string of the parameter chains to plot
}
  \item{\dots}{
further arguments to pass to \code{summary.mcmc}
}
}
\details{
\code{chain} should be one of \code{'R'} (for correlation coefficients) or \code{'B$sp'} where \code{sp} is the species of interest (for regression coefficients).
}

\seealso{
\code{\link{BC}}, \code{\link{summary.mcmc}}
}
\examples{
m1 <- example(BC)[[1]]
summary(m1, 'R')
summary(m1, 'B$sp1')
}

