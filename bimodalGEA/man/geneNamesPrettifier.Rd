\name{geneNamesPrettifier}
\alias{geneNamesPrettifier}
\alias{disambiguate}
\title{fix genes naming schemes}
\description{
  fix genes naming schemes
}
\usage{
geneNamesPrettifier(geneSymbols, geneIds)
disambiguate(x, disambiguator=seq_along, sep=" ")
}
\arguments{
  \item{geneSymbols}{vector of human-friendly gene symbols}
  \item{geneIds}{vector of more unique and stable genes ids}
  \item{x}{vector of names}
  \item{disambiguator}{func producing suffixes to resolve ambiguities}
  \item{sep}{passed to \code{\link{paste}}}
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
## With disambiguator, we give different suffixes to each replica:
disambiguate(geneNames)
## [1] "a 1" "b"   "a 2" "c"
disambiguate(geneNames, sep="/")
## [1] "a/1" "b"   "a/2" "c"
}
\keyword{ util }
\keyword{ misc }
