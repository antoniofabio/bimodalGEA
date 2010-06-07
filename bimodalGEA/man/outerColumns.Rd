\name{outerColumns}
\alias{outerColumns}
\title{outer product between two vertical matrices}
\description{column by column operations on two vertical matrices}
\usage{
outerColumns(A, B=A, FUN)
}
\arguments{
  \item{A, B}{two vertical matrices. They must have the same number of rows}
  \item{FUN}{a two-arguments function}
}
\value{
  Returns a matrix in which the element (i,j) equals \code{FUN(A[,i], B[,j])}.
}
\author{
Antonio, Fabio Di Narzo
}
\examples{
A <- matrix(1:2, 1, dimnames=list(NULL, letters[1:2]))
B <- matrix(1:3, 1, dimnames=list(NULL, letters[3:5]))
outerColumns(A, B, "*")
#   c d e
# a 1 2 3
# b 2 4 6

outerColumns(1, B, "*")
# c d e 
# 1 2 3

## arguments can be dataframes too:
outerColumns(as.data.frame(A), as.data.frame(B), "*")
#   c d e
# a 1 2 3
# b 2 4 6
}
\keyword{ programming }
\keyword{ utilities }
