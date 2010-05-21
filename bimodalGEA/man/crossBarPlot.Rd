\name{crossBarPlot}
\alias{crossBarPlot}
\alias{crossBarPlotStrata}
\title{crossBarPlot}
\description{Barplots of discretized gene expression vs a factor variable}
\usage{
crossBarPlot(x, fac, xlab = "", legendTitle = NULL, legend.text = TRUE,
  args.legend = list(title = legendTitle), counts = crossTable(x, fac),
  beside = TRUE, ...)
crossBarPlotStrata(strata, x, fac, main, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{gene expression values}
  \item{fac}{factor}
  \item{xlab, legendTitle, legend.text, beside, args.legend, main, ...}{
    arguments passed to \code{\link{barplot}}
  }
  \item{counts}{
    crossed table frequencies
  }
  \item{strata}{stratification variable}
}
\value{ Nothing useful. }
\author{
  Antonio, Fabio Di Narzo
}

\examples{
## TODO
}
\keyword{ hplot }
