\name{mixFit}
\alias{mixFit}
\alias{mixPar}
\title{
  Fit univariate 2-components mixture model
}
\description{
  Fit univariate 2-components mixture model
}
\usage{
mixFit(x)
mixPar(m)
}
\arguments{
  \item{x}{vector of values on which to fit the mixture model}
  \item{m}{already fitted mixture model, as returned by \code{mixFit}}
}
\details{
  Fit univariate 2-components mixture model.
  Uses the \code{Mclust} function in the \code{mclust} package.
}
\value{
  \code{mixFit} returns the fitted model object.
  \code{mixPar} extract the mixture parameters from the mixture model object.
}
\author{
  Antonio, Fabio Di Narzo.
}
\examples{
set.seed(1234)
x <- rnorm(100, sample(c(-2, 2), size=100, replace=TRUE))
mixPar(mixFit(x))$mean
#        1         2 
#-2.018164  2.039135 
}
\keyword{ cluster }
\keyword{ univar }
\keyword{ models }
