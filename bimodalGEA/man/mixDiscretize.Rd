\name{mixDiscretize}
\alias{mixDiscretize}
\alias{classifyBy}
\alias{mixDiscretizeStrata}
\title{
Discretize gene expression levels
}
\description{
  Discretize gene expression levels into 'high'/'low'.
}
\usage{
mixDiscretize(X)
mixDiscretizeStrata(X, strata)
}
\arguments{
  \item{X}{ gene expression matrix (genes in columns) }
  \item{strata}{stratification variable}
}
\details{
  Discretize gene expression levels into 'high'/'low', basing on
  2-components mixture fit.

  \code{mixDiscretizeStrata} does that separately for each level of
  'strata', and finally merges all the values back together.
}
\value{
  A matrix of factor variables, one per gene, with levels 'high' and 'low'.
}
\author{
  Antonio, Fabio Di Narzo
}
\examples{
set.seed(1234)
X <- matrix(ifelse(sample(0:1, size=200, repl=TRUE),
rnorm(200, -2), rnorm(200, 2)), 100)
Y <- mixDiscretize(X)

tapply(X[,1], Y[,1], range)
% $high
% [1] -0.0162678  3.9130594

% $low
% [1] -4.8557587 -0.3160362

\dontshow{
mm <- tapply(X[,1], Y[,1], mean)
stopifnot(all(names(mm) \%in\% c("high", "low")))
stopifnot(mm['high'] >= mm['low'])
}

}
\keyword{ manip }
\keyword{ cluster }
