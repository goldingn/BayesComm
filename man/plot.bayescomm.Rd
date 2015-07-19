\name{plot.bayescomm}
\alias{plot.bayescomm}

\title{
Plot bayescomm parameter chains
}
\description{
\code{plot.bayescomm} creates summary plots of a subset of the parameter chains in a \code{bayescomm} object using the \code{coda} package.
}
\usage{
\method{plot}{bayescomm}(x, chain, ...)
}

\arguments{
  \item{x}{
a \code{bayescomm} object
}
  \item{chain}{
a character string of the parameter chains to plot
}
  \item{\dots}{
further arguments to pass to \code{plot.mcmc}
}
}
\details{
\code{chain} should be one of \code{'R'} (for correlation coefficients) or \code{'B$sp'} where \code{sp} is the species of interest (for regression coefficients).
}
\seealso{
\code{\link{BC}}, \code{\link{plot.mcmc}}
}
\examples{
m1 <- example(BC)[[1]]
plot(m1, 'R')
plot(m1, 'B$sp1')
}