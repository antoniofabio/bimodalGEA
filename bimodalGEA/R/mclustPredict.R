mclustLik <- function(mc, y) {
  pp <- mc$parameters
  G <- pp$variance$G
  N <- length(y)
  yy <- rep(y, G)
  mm <- rep(pp$mean, each = N)
  ssd <- sqrt(pp$variance$sigmasq)
  if(length(ssd) == 1) {
    ssd <- rep(ssd, G)
  }
  ss <- rep(ssd, each = N)
  return(matrix(dnorm(y, mm, ssd), nrow = N))
}

mclustPredict <- function(mc, y) {
  z <- sweep(mclustLik(mc, y), 2, mc$parameters$pro, "*")
  return(apply(z, 1, which.max))
}
