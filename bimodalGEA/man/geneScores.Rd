\name{geneScores}
\alias{geneScores}
\alias{subsetScores}
\alias{strataScores}
\title{
  Compute gene-specific bimodality scores
}
\description{
    Compute gene-specific bimodality scores
}
\usage{
geneScores(X, F=Rsquare)
subsetScores(X, F, subset, strata)
strataScores(X, F, strata)
}
\arguments{
  \item{X}{ gene expression matrix (genes are in columns) }
  \item{F}{ score function (see the \code{\link{scoreFuns}} man page) }
  \item{strata}{stratification variable}
  \item{subset}{subsetting variable}
}
\details{
  Compute gene-specific bimodality scores.

  \code{strataScores} does it separately for each sample strata, and
  then averages the results.

  \code{subsetScores} does it separately for each sample strata and subsets, and
  then averages the results within each subset.
}
\value{
  The vector of computed bimodality scores.
}
\author{
  Antonio, Fabio Di Narzo
}
\seealso{
  \code{\link{scoreFuns}}, for a list of the available score functions.
}
\examples{
set.seed(1234)
x1 <- rnorm(100, sample(c(-2, 2), size=100, replace=TRUE))
x2 <- rnorm(100)
X <- cbind(x1, x2)
geneScores(X)
#        x1         x2 
#0.81370903 0.04887823 
}
\keyword{ cluster }
\keyword{ univar }
\keyword{ manip }
