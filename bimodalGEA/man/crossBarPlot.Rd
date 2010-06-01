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
\seealso{
  \code{\link{plotTopMosaics}}, \code{\link{plotTopDensities}},
  \code{\link{plotTopJitters}}
}
\examples{
set.seed(1234)
f <- sample(c('m','f'), size=100, replace=TRUE)
x <- rnorm(100, sample(c(-2, 2), size=100, replace=TRUE))
crossBarPlot(x, f)

crossBarPlot(x, f, legendTitle="sex")
}
\keyword{ hplot }
