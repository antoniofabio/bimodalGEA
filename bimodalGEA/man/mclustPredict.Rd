\name{mclustPredict}
\alias{mclustPredict}
\alias{mclustLik}
\title{
  Predict clusters from univariate mixture fit
}
\description{
  Predict clusters from univariate mixture fit
}
\usage{
mclustPredict(mc, y)
mclustLik(mc, y)
}
\arguments{
  \item{mc}{\code{\link{Mclust}} mixture fit}
  \item{y}{vector of observations for which to predict the cluster}
}
\details{
  Predict clusters from univariate Gaussian mixture fit, as computed by
  the \code{\link{Mclust}} function.
}
\value{
  \code{mclustPredict} returns a vector of group indexes.

  \code{mclustLik} returns a matrix with one row for each element of
  \code{y}, and one column for each fitted mixture component, containing
  the likelihood of each mixture component for each observation.
}
\author{
  Antonio, Fabio Di Narzo.
}
\examples{
if(require(Mclust)) {
  set.seed(1234)
  y <- rnorm(20)
  mc <- mclust(y, G = 2)
  mclustPredict(mc, y)
}
}
\seealso{ \code{\link{Mclust}} }
\keyword{ cluster }
\keyword{ univar }
\keyword{ models }
