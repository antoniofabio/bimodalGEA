##
## Bi-modality based gene expression analysis
## Auxiliary density plotting functions
##

.formatMixParameters <- function(m) {
  list(p.txt = paste(paste("p={",
         paste(format(m$parameters$pro, digits=2),collapse=";"),
         "}", sep=""), collapse="; "),
       mu.txt = paste(paste("m={",
         paste(format(m$parameters$mean, digits=2),collapse=";"),
         "}", sep=""), collapse="; "),
       sigma.txt = paste("s=",
         format(sqrt(m$parameters$variance$sigmasq[1]), digits=2),
         sep=""))
}

.plotMixParameters <- function(m, adj=0.01, cex=0.7, ...) {
  with(.formatMixParameters(m), {
    mtext(mu.txt, line=-1, adj=adj, cex=cex, ...)
   mtext(p.txt, line=-2, adj=adj, cex=cex, ...)
    mtext(sigma.txt, line=-3, adj=adj, cex=cex, ...)
  })
}

## Plot densities of gene expressions (columns of X)
## sorted by scores. Superimpose mixture estimates.
## Scores and columns of X are matched by 'names'.
plotTopDensities <- function(X, score, NR=4, NC=7,
                             PAGG=ceiling(length(score)/(NR*NC))) {
  X <- X[,names(score)]
  xi <- order(score, decreasing=TRUE)
  par(mfrow=c(NR,NC), mar=c(2,2,4.2,1))
  for(i in seq_len(PAGG*NR*NC)) {
    yi <- X[,xi[i]]
    symb <- colnames(X)[xi[i]]
    di <- density(yi)
    plot(di, main=symb, lwd=2, col="darkgray")
    .plotMixParameters(mixFit(yi))
    score.txt <- paste("score=",format(score[xi[i]], digits=2))
    mtext(score.txt, line=-1, adj=0.99, cex=0.7)
  }
}
