plotTopMosaics <- function(X, NR=4, NC=7, PAGG=ceiling(NCOL(X)/(NR*NC)),
                           color=TRUE, ...) {
  X.discrete <- mixDiscretize(X)
  par(mfrow=c(NR,NC), mar=c(2,2,4.2,1))
  nms <- colnames(X)
  for(il in seq(1, min(length(nms), NR*NC*PAGG))) {
    l <- nms[il]
    corrs <- rep(NA, length(nms))
    names(corrs) <- nms
    for(m in setdiff(nms, l))
      corrs[m] <- fisherCorrelation(X.discrete[,l], X.discrete[,m])
    m <- nms[which.min(corrs)]
    tab <- table(X.discrete[,l], X.discrete[,m])
    tab <- tab[c('low', 'high'), c('high', 'low')]
    mosaicplot(tab,
               main=format(-log10(corrs[m]), digits=3),
               xlab=l, ylab=m, color=color)
  }
}
