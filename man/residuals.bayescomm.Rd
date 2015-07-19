\name{residuals.bayescomm}
\alias{residuals.bayescomm}

\title{
Extract bayescomm model residuals
}
\description{
\code{residuals.bayescomm} extracts model residuals from a \code{bayescomm} object.
Residuals are calculated based on the mean of the posterior probability of presence.
}
\usage{
\method{residuals}{bayescomm}(object, ...)
}

\arguments{
  \item{object}{
a \code{bayescomm} object
}
  \item{\dots}{
other arguments
}
}

\seealso{
\code{\link{BC}}
}
\examples{
m1 <- example(BC)[[1]]
m1.res <- residuals(m1)
}

