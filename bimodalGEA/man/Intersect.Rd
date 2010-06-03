\name{setOperations}
\alias{Intersect}
\alias{Union}
\title{Intersection and Union of list elements}
\description{Intersection and Union of list elements}
\usage{
Intersect(lst)
Union(lst)
}
\arguments{
  \item{lst}{a \code{list}}
}
\value{
  \code{Intersect} returns the vector of elements present in all the
  list slots.

  \code{Union} returns the vector of elements present in at least one the
  list slots.
}
\author{
Antonio, Fabio Di Narzo
}
\examples{
lst <- list(letters[1:5], letters[4:8], letters[4:9])
Intersect(lst)
## [1] "d" "e"
Union(lst)
## [1] "a" "b" "c" "d" "e" "f" "g" "h" "i"

\dontshow{
stopifnot(identical(Union(list(1:2, 3:4)), 1:4))
stopifnot(Intersect(list(1:3, 3:5)) == 3)
}

}
\keyword{ programming }
\keyword{ utilities }
