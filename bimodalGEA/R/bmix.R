bmix <- function(y, n, n.burnin = 0, thin = 1) {
  y <- y - median(y)
  x0 <- with(mixPar(mixFit(y)), c(beta = pro[1],
                                  mL = mean[1],
                                  mH = mean[2],
                                  s2 = variance$sigmasq))
  reproject <- function(x, l, u) max(min(x, u), l)
  x0[1] <- reproject(x0[1], 0.1, 0.9)
  x0[2] <- reproject(x0[2], -10, 10)
  x0[3] <- reproject(x0[3], -10, 10)
  x0[4] <- reproject(x0[4], 0.001, 10000)
  S0 <- diag(c(0.1, 1.0, 1.0, 1.0))
  .C("bmix_init",
     as.double(y),
     as.integer(length(y)),
     as.double(x0),
     as.double(S0),
     as.integer(n.burnin),
     PACKAGE = "bimodalGEA")
  on.exit(.C("bmix_free", PACKAGE = "bimodalGEA"))
  out <- rep(-1.0, n + n.burnin)
  out <- .C("bmix_update",
            as.integer(n + n.burnin),
            as.integer(thin),
            out = as.double(out),
            PACKAGE = "bimodalGEA")$out
  return(if(n.burnin > 0) tail(out, -n.burnin) else out)
}
