\name{bmix}
\alias{bmix}
\title{
  Bayesian estimation of R-squared bimodality index
}
\description{
  Bayesian estimation of R-squared bimodality index,
  with weakly informative priors on model parameters.
}
\usage{
bmix(y, n, n.burnin = 0, thin = 1)
}
\arguments{
  \item{y}{observed values}
  \item{n}{number of MCMC batches}
  \item{n.burnin}{number of MCMC iterations to be discarded}
  \item{thin}{thinning interval}
}
\details{
  Bayesian estimation of R-squared bimodality index,
  with weakly informative priors on model parameters.

  The parameters posterior distribution is estimated using
  an adaptive random walk MCMC algorithm.
}
\value{
  Trace of the simulated R-squared values, after dropping the burn-in iterations.
}
\author{
  Antonio, Fabio Di Narzo.
}
\examples{
set.seed(1234)
y <- rnorm(100, sample(c(-2, 2), size=100, replace=TRUE))
x <- bmix(y, n = 10000)
plot(ts(x))
abline(h = 0.8, lty = 2)
}
\keyword{ cluster }
\keyword{ univar }
\keyword{ models }
