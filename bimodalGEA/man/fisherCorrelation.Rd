\name{fisherCorrelation}
\alias{fisherCorrelation}
\alias{fisherCorrelations}
\title{
Fisher test between two categorical variables
}
\description{
  Fisher test between two categorical variables (or groups of them)
}
\usage{
fisherCorrelation(x, y)
fisherCorrelations(X, Y)
}
\arguments{
  \item{x, y}{two discrete variables vectors}
  \item{X}{vertical matrix of discrete variables}
  \item{Y}{vertical matrix of discrete variables}
}
\value{
  \code{fisherCorrelation} returns the Fisher test p-value between two
  categorical variables.

  \code{fisherCorrelations} returns a matrix containing all the
  pairwise tests p-values. If \code{Y} is missing, returns the symmetric
  square matrix of all the pairwise tests p-values (\code{NA} on the
  main diagonal).
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

fisherCorrelations(a, cbind(a, b))
#           a            b 
#1.627380e-29 6.879004e-01 
}
\keyword{ htest }
