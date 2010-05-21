\name{getPowerScore}
\alias{getPowerScore}
\title{
  proportion-penalized bi-modality score
}
\description{
  proportion-penalized bi-modality score
}
\usage{
getPowerScore(pow)
}
\arguments{
  \item{pow}{ penalization exponent }
}
\details{
  Returns a new power score function, obtained by penalizing the
  \code{\link{Rsquare}} function for estimated proportions which are far
  from the perfect equilibrium \code{c(0.5,0.5)}.
  
}
\value{
  A new bimodality score function, with the desired power penalty.
}
\author{
Antonio, Fabio Di Narzo.
}

\seealso{
  \code{\link{geneScores}}, \code{\link{scoreFuns}}
}
\examples{
set.seed(1234)
x1 <- rnorm(100, sample(c(-1, 1), size=100, replace=TRUE, prob=c(0.8,
  0.2)))
x2 <- rnorm(100, sample(c(-1.05, 1.05), size=100, replace=TRUE))
X <- cbind(x1, x2)
geneScores(X)
#       x1        x2 
#0.3902909 0.3451000 
geneScores(X, getPowerScore(2))
#       x1        x2 
#0.2108463 0.2310068 
geneScores(X, getPowerScore(4))
#       x1        x2 
#0.1139052 0.1546338 
}
\keyword{ math }
\keyword{ univar }
\keyword{ cluster }
