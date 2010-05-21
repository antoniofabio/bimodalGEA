\name{lmClean}
\alias{lmClean}
\title{
  Remove uninteresting covariates effects from a gene expression matrix
}
\description{
  Remove uninteresting covariates effects from a gene expression matrix
  by linear model fitting.
}
\usage{
lmClean(X, formula, data, X.design = model.matrix(formula, data))
}
\arguments{
  \item{X}{
    gene expression matrix (genes in columns)
  }
  \item{formula}{ right hand side of the linear model formula to be used
  to remove the covariates effects }
\item{data}{
  data frame containing the variables used in the \code{formula}
  arguments, with the same number of rows as \code{X}.
}
  \item{X.design}{linear model design matrix}
}
\details{
  Remove uninteresting covariates effects from a gene expression matrix
  by linear model fitting.

  The linear model is fitted separately for each gene. The model
  residuals are extracted, and re-shifted to have the original data
  overall mean.
}
\value{
  The covariate-effects corrected gene expression matrix, in the
  original data scale.
}
\author{
  Antonio, Fabio Di Narzo.
}
\examples{
df <- data.frame(sex = factor(rep(c('m','f'), c(10,10))))
set.seed(1234)
X.m <- matrix(rnorm(20, -2), 10)
X.f <- matrix(rnorm(20, 2), 10)
X <- rbind(X.m, X.f)
Y <- lmClean(X, ~ sex, data=df)
range(X)
# [1] -4.345698  3.102298
range(Y)
# [1] -2.348092  2.091824
coef(lm(X[,1] ~ df$sex))
#(Intercept)     df$sexm 
#   1.612053   -3.995211 
coef(lm(Y[,1] ~ df$sex))
#  (Intercept)       df$sexm 
#-3.855521e-01  9.930137e-17 
}
\keyword{ manip }
\keyword{ univar }
