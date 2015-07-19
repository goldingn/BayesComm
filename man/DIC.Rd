\name{DIC}
\alias{DIC}

\title{
Deviance Information Criterion
}
\description{
Calculates Deviance Information Criteria for \code{bayescomm} objects.
}
\usage{
DIC(BC)
}

\arguments{
  \item{BC}{
a \code{bayescomm} object
}
}

\references{
Spiegelhalter, D.J.. Best, N.G., Carlin, B.P., van der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society, Series B}, 64 (4): 583-639.
}

\seealso{
\code{\link{BC}}
}

\examples{
m1 <- example(BC)[[1]]
DIC(m1)
}
