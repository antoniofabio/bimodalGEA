\name{plotTopDensities}
\alias{plotTopDensities}
\title{plotTopDensities}
\description{ Plot top-scores distribution densities }
\usage{
plotTopDensities(X, score,
  NR = 4, NC = 7,
  PAGG = ceiling(length(score)/(NR * NC)))
}
\arguments{
  \item{X}{gene expression matrix (genes in columns)}
  \item{score}{named vector of score values. Vector names must match
    \code{X} column names}
  \item{NR}{number of figure rows}
  \item{NC}{number of figure columns}
  \item{PAGG}{number of pages to be produced}
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
  \code{\link{geneScores}}
}
\examples{
set.seed(1234)
X <- matrix(rnorm(400), 100)
colnames(X) <- letters[seq_len(NCOL(X))]
scores <- geneScores(X)
plotTopDensities(X, scores, NR=1, NC=2, PAGG=1)
}
\keyword{ hplot }
