##
## Bi-modality based gene expression analysis
## Auxiliary jitter plotting functions
##

plotMixJitter <- function(x, score, col="black", method="jitter",
                          ylim=range(x), vertical=TRUE, ...) {
  cl <- mixDiscretize(x)
  stripchart(x[cl == 'low'], col=col, method=method, pch='l', ylim=ylim,
             vertical=vertical, ...)
  stripchart(x[cl == 'high'], col=col, method=method, pch='h',
             ylim=ylim, vertical=vertical, add=TRUE, ...)
  .plotMixParameters(mixFit(x))
  if(!missing(score)) {
    score.txt <- paste("score=",format(score, digits=2))
    mtext(score.txt, line=-1, adj=0.99, cex=0.7)
  }
}

plotTopJitters <- function(X, score, NR=4, NC=7,
                           PAGG=ceiling(length(score)/(NR*NC)),
                           ylim=range(X)) {
  if(is.matrix(score)) {
    nms <- rownames(score)
    score <- as.vector(score)
    names(score) <- nms
  }
  X <- X[,names(score)]
  xi <- seq_along(score)
  par(mfrow=c(NR,NC), mar=c(2,2,4.2,1))
  for(i in seq_len(min(NCOL(X), PAGG*NR*NC))) {
    symb <- colnames(X)[xi[i]]
    plotMixJitter(X[,symb],
                  score=score[symb],
                  main=symb,
                  xlab="", ylab="", ylim=ylim)
  }
}
