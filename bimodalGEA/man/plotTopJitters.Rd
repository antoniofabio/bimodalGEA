\name{plotTopJitters}
\title{ Plot top-scores distributions }
\alias{plotTopJitters}
\alias{plotMixJitter}
\description{ Plot top-scores distributions }
\usage{
plotMixJitter(x, score, col="black", method="jitter",
              ylim=range(x), vertical=TRUE, ...)

plotTopJitters(X, score, NR=4, NC=7,
              PAGG=ceiling(length(score)/(NR*NC)),
              ylim=range(X))
}
\arguments{
  \item{x}{gene expression vector}
  \item{score}{pre-computed score (for \code{plotMixDensity}) or a named
    vector of score values (for \code{plotTopDensities}). Vector names must match
    \code{X} column names}
  \item{col, method, ylim, vertical, ...}{further arguments to \code{\link{stripchart}}}
  \item{X}{gene expression matrix (genes in columns)}
  \item{NR}{number of figure rows}
  \item{NC}{number of figure columns}
  \item{PAGG}{number of pages to be produced}
}
\details{
  dotchart plots of the top-scoring genes, with symbol basing on the
  mixture classification.
}
\value{ Nothing interesting. }
\author{
  Antonio, Fabio Di Narzo
}
\seealso{
  \code{\link{geneScores}}, \code{\link{plotTopDensities}}
}
\examples{
set.seed(1234)
X <- matrix(rnorm(400), 100)
colnames(X) <- letters[seq_len(NCOL(X))]
scores <- geneScores(X)
plotTopJitters(X, scores, NR=1, NC=2)
}
\keyword{ hplot }
