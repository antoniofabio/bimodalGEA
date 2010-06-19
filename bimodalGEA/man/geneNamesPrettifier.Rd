\name{geneNamesPrettifier}
\alias{geneNamesPrettifier}
\alias{disambiguate}
\title{fix genes naming schemes}
\description{
  fix genes naming schemes
}
\usage{
geneNamesPrettifier(geneSymbols, geneIds)
}
\arguments{
  \item{geneSymbols}{vector of human-friendly gene symbols}
  \item{geneIds}{vector of more unique and stable genes ids}
}
\value{
  the modified vector of gene names
}
\author{
  Antonio, Fabio Di Narzo
}
\examples{
geneNames <- c('a', 'b', 'a', 'c')
geneIds <- seq_along(geneNames)
geneNamesPrettifier(geneNames, geneIds)
## [1] "a (1)" "b"     "a (3)" "c"
}
\keyword{ util }
\keyword{ misc }
