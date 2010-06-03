\name{isUnique}
\alias{isUnique}
\title{Check if a vector is made of unique values only}
\description{Check if a vector is made of unique values only}
\usage{
isUnique(x)
}
\arguments{
  \item{x}{a vector}
}
\value{
  \code{TRUE} if the vector is made of unique values only and contains
  no missing values, \code{FALSE} otherwise.
}
\author{
Antonio, Fabio Di Narzo
}
\examples{
isUnique(1:3)
## [1] TRUE
isUnique(c(1,2,NA))
## [1] FALSE
isUnique(c(1,2,2))
## [1] FALSE

\dontshow{
stopifnot(isUnique(1:3))
stopifnot(!isUnique(c(1,2,NA)))
stopifnot(!isUnique(c(1,2,2)))
}

}
\keyword{ manip }
\keyword{ programming }
\keyword{ utilities }
