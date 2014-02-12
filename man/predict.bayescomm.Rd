\name{predict.bayescomm}
\alias{predict.bayescomm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make predictions at new locations
}
\description{
For each set of parameter values sampled by the model (including values of Z),
simulate the occurrence probabilities for each species at each new location.
}
\usage{
predict.bayescomm(object, newdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A bayescomm object
}
  \item{newdata}{
A data.frame with the same columns as X from the original BC model
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An array of occurrence probabilities.
Rows index locations.  Columns index species.  Slices index MCMC samples.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David J. Harris (http://davharris.github.io)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
