\name{Compose}
\alias{Compose}
\title{Compose two or more functions in a new function}
\description{Compose two or more functions in a new function}
\usage{
Compose(...)
}
\arguments{
  \item{\dots}{comma separated list of the functions to be composed}
}
\details{
  Compose a list of functions together into a new function, which is returned.
}
\value{
  If called with \code{Compose(f, g, h)}, returns the function which
  computes \code{f(g(h(x)))}.
}
\author{
Antonio, Fabio Di Narzo
}
\examples{
## define the 'sum of squares' function as the composition
## of the 'square' and 'sum' functions:
square <- function(x) x^2
sumOfSquares <- Compose(sum, square)
sumOfSquares(1:3)
## [1] 14

Compose("!", any, is.na)(c(1, NA, 2))
## [1] FALSE

\dontshow{
stopifnot(sumOfSquares(1:3) == 14)
}

}
\keyword{ programming }
\keyword{ utilities }
