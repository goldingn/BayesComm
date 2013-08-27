\name{window.bayescomm}
\alias{window.bayescomm}

\title{
Window bayescomm objects
}
\description{
\code{window.bayescomm} is window function for \code{bayescomm} objects, it calls \code{window.mcmc} from the \code{coda} package.
Parameter chains are subsetted by \code{start} and \code{end} and thinned by \code{thin}.
}
\usage{
\method{window}{bayescomm}(x, start = NULL, end = NULL, thin = 1, ...)
}

\arguments{
  \item{x}{
a \code{bayescomm} object
}
  \item{start}{
start iteration
}
  \item{end}{
end iteration
}
  \item{thin}{
thinning interval
}
  \item{\dots}{
further arguments to pass to \code{window.mcmc}
}
}

\details{
If \code{start = NULL} (default) the start is taken as the first iteration.
If \code{end = NULL} (default) the end is taken as the  final iteration.
If \code{thin = 1} (default) all iterations within the window are retained.
}

\value{
A \code{bayescomm} object with windowed parameter chains.
}
\seealso{
\code{\link{BC}}, \code{\link{window.mcmc}}
}
\examples{
m1 <- example(BC)[[1]]
m2 <- window(m1, 51, 150, 10) 
}

