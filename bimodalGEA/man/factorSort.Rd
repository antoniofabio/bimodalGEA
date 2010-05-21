\name{factorSort}
\alias{factorSort}
\title{re-sort a factor levels by frequency}
\description{
  re-sort a factor levels by frequency
}
\usage{
factorSort(fac)
}
\arguments{
  \item{fac}{the factor variable to be recoded}
}
\value{
  the recoded factor variable
}
\author{
  Antonio, Fabio Di Narzo
}
\examples{
sex <- rep(c('m','f'), c(10, 5))
table(sex)
#sex
# f  m 
# 5 10
table(factorSort(sex))
#
# m  f 
#10  5
}
\keyword{ util }
\keyword{ misc }
