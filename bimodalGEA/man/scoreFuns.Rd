\name{Bi-modality score functions}
\alias{Rsquare}
\alias{BI.wang}
\alias{Mscore}
\alias{meansDistanceScore}
\alias{minProportion}
\alias{scoreFuns}
\title{
  Bi-modality score functions
}
\description{
  Bi-modality score functions
}
\usage{
Rsquare(pp)
BI.wang(pp)
Mscore(pp)
meansDistanceScore(pp)
scoreFuns
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pp}{fitted mixture parameters, as returned by \code{\link{mixFit}}}
}
\details{
    Bi-modality score functions. Look at the functions definitions for
    their respective formulas.

    \code{scoreFuns} is a list of the bimodality score functions
    currently available.
}
\value{
  The computed bi-modality score
}
\author{
Antonio, Fabio Di Narzo
}
\seealso{
  \code{\link{geneScores}}, \code{\link{getPowerScore}}
}
\examples{
## see 'geneScores' man page
}
\keyword{ math }
\keyword{ univar }
\keyword{ cluster }
