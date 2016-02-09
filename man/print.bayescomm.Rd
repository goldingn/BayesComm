\name{print.bayescomm}
\alias{print.bayescomm}

\title{
Print a bayescomm object
}
\description{
\code{print.bayescomm} prints a brief summary of a \code{bayescomm} object.
}
\usage{
\method{print}{bayescomm}(x, ...)
}

\arguments{
  \item{x}{
a \code{bayescomm} object
}
  \item{\dots}{
further arguments to pass to \code{print}
}
}

\seealso{
\code{\link{BC}}
}
\examples{
m1 <- example(BC)[[1]]
print(m1)
m1
}
