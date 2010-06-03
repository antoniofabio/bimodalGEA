\name{plotTopDensities}
\alias{plotTopDensities}
\alias{plotMixDensity}
\title{plotTopDensities}
\description{ Plot top-scores distribution densities }
\usage{
plotMixDensity(x, score, lwd=2, col="darkgray", ...)

plotTopDensities(X, score,
  NR = 4, NC = 7,
  PAGG = ceiling(length(score)/(NR * NC)),
  xlim=range(X), addRug=FALSE)
}
\arguments{
  \item{x}{gene expression vector}
  \item{score}{pre-computed score (for \code{plotMixDensity}) or a named
    vector of score values (for \code{plotTopDensities}). Vector names must match
    \code{X} column names}
  \item{lwd, col, xlim, ...}{further arguments to \code{plot}}
  \item{X}{gene expression matrix (genes in columns)}
  \item{NR}{number of figure rows}
  \item{NC}{number of figure columns}
  \item{PAGG}{number of pages to be produced}
  \item{addRug}{calls \code{\link{histSpike}} on every plot}
}
\details{
  Plot the kernel-estimated distribution density of the top-scoring
  genes in \code{X}.
}
\value{ Nothing interesting. }
\author{
  Antonio, Fabio Di Narzo
}
\seealso{
  \code{\link{geneScores}}, \code{\link{plotTopJitters}},
  \code{\link{crossBarPlot}}
}
\examples{
set.seed(1234)
X <- matrix(rnorm(400), 100)
colnames(X) <- letters[seq_len(NCOL(X))]
scores <- sort(geneScores(X), decreasing=TRUE)
plotTopDensities(X, scores, NR=1, NC=2, addRug=TRUE)

## A more fine-grained control:
par(mar=c(2,2,4.2,1), mfrow=c(2,1))
plotMixDensity(X[,'a'], scores['a'], main="gene 'a'", xlab="", ylab="",
  xlim=range(X))
plotMixDensity(X[,'b'], scores['b'], main="gene 'b'", xlab="", ylab="",
  col="pink", xlim=range(X))
}
\keyword{ hplot }
