\name{plotTopMosaics}
\alias{plotTopMosaics}
\title{
  Plot most 'associated' mosaics plots from an expression matrix
}
\description{
    Plot most 'associated' mosaics plots from an expression matrix
}
\usage{
plotTopMosaics(X, NR = 4, NC = 7, PAGG = ceiling(NCOL(X)/(NR * NC)),
  color = TRUE, ...)
}
\arguments{
  \item{X}{ gene expression matrix (genes in columns) }
  \item{NR, NC, PAGG}{num. rows, columns and pages}
  \item{color, ...}{extra arguments to \code{\link{mosaicplot}}}
}
\value{
  None interesting.

  A plot is produced, with a number of panels.
  In each panel, there is the mosaic plot between the ith gene and the
  most associated (according to a Fisher test p-value) other gene.
  The panel title displays \code{-log10(p.value)} (Fisher test p-value).
  The expression profile of each gene is discretized using \code{\link{mixDiscretize}}.
}
\author{
  Antonio, Fabio Di Narzo
}
\examples{
set.seed(1234)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
plotTopMosaics(cbind(x1, x2, x3), 1, 3)
}
\keyword{ hplot }
