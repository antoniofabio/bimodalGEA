\name{geneScores}
\alias{bimodalityScores}
\alias{bimodalityFilter}
\alias{geneScores}
\alias{subsetScores}
\alias{strataScores}
\title{
  Compute gene-specific bimodality scores
}
\description{
  Compute gene-specific bimodality scores, and filter genes basing on
  custom thresholds.
}
\usage{
bimodalityScores(X)
bimodalityFilter(X, meansDistance=2, proportion=0.01, Rsquare=0)

geneScores(X, ...)
subsetScores(X, F, subset, strata)
strataScores(X, F, strata)
}
\arguments{
  \item{X}{ gene expression matrix (genes are in columns) }
  \item{meansDistance}{ min allowed distance between mixture means }
  \item{proportion}{ min allowed mixture proportion }
  \item{Rsquare}{ min allowed Rsquare value }
  \item{...}{ comma separated list of score functions (see the \code{\link{scoreFuns}} man page) }
  \item{F}{ score function (see the \code{\link{scoreFuns}} man page) }
  \item{strata}{stratification variable}
  \item{subset}{subsetting variable}
}
\details{
  Compute gene-specific bimodality scores, and filter genes basing on
  custom thresholds over distance between means and over minimum
  estimated proportion for a mixture component.

  \code{bimodalityScores} computes \code{Rsquare}, distance between
  mixture means and minimum proportion of the mixture components, for
  each gene.

  \code{bimodalityFilter} calls \code{bimodalityScores} and then selects
  the gene names which met the specified requirements.

  \code{strataScores} does it separately for each sample strata, and
  then averages the results.

  \code{subsetScores} does it separately for each sample strata and subsets, and
  then averages the results within each subset.
}
\value{
  \code{bimodalityScores} returns a \code{data.frame} containing gene
  names, \code{Rsquare} index, distance between means and minimum
  mixture proportion as columns.

  \code{bimodalityFilter} filters the genes in \code{X} removing those
  with too small distance beween mixture means, too small proportion
  of one of the mixture components or too small explained variance
  proportion (Rsquare).
  
  \code{geneScores}, \code{subsetScores} and \code{strataScores} return
  the vector of computed bimodality scores.
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
#            F
#x1 0.81370903
#x2 0.04887823

## Compute multiple indexes at once:
geneScores(X,
  Rsq=Rsquare,
  mD=meansDistanceScore,
  mP=minProportion)
##           Rsq         mD        mP
## x1 0.81370903 0.02097167 0.4677786
## x2 0.04887823 0.23773018 0.4630958

## bimodality-based filtering of genes:
bimodalityScores(X)
bimodalityFilter(X, 0.1, 0.1)
}
\keyword{ cluster }
\keyword{ univar }
\keyword{ manip }
