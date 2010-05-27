\name{fisherCorrelation}
\alias{fisherCorrelation}
\alias{fisherCorrelations}
\title{
Fisher test between two categorical variables
}
\description{
  Fisher test between two categorical variables
}
\usage{
fisherCorrelation(xa, xb)
fisherCorrelations(X.discrete)
}
\arguments{
  \item{xa, xb}{two variables vectors}
  \item{X.discrete}{vertical matrix of discrete variables}
}
\details{
  \code{fisherCorrelations} returns the square matrix of the pairwise
  tests p-values.
}
\value{
  \code{fiherCorrelation} returns the Fisher test p-value.
  \code{fiherCorrelations} returns a square matrix containing all the
  pairwise tests p-values (\code{NA} on the main diagonal).
}
\author{
  Antonio, Fabio Di Narzo
}

\examples{
set.seed(1234)
a <- sample(letters[1:2], size=100, replace=TRUE)
b <- sample(letters[3:4], size=100, replace=TRUE)
fisherCorrelation(a, b)
# [1] 0.6879004

fisherCorrelations(cbind(a, b))
#          a         b
#a        NA 0.6879004
#b 0.6879004        NA
}
\keyword{ htest }
