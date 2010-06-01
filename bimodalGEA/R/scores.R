##
## Bi-modality based gene expression analysis
## Antonio F. Di Narzo, 20 May 2010
##
mixFit <- function(x) Mclust(data=x, G=2, modelNames="E")
mixPar <- function(m) m$parameters
.mixDiscretize <- function(x) {
  m <- mixFit(x)
  ix <- sort(m$parameters$mean, index=TRUE, decreasing=TRUE)$ix
  return(factor(m$classification, levels=1:2, labels=c("high", "low")[ix]))
}
delta <- function(pp) abs(diff(pp$mean))/sqrt(pp$variance$sigmasq)

BI.wang <- function(pp) delta(pp) * sqrt(prod(pp$pro))

Mscore <- function(pp) sqrt(delta(pp)) * prod(pp$pro)

varBetween <- function(pp) prod(pp$pro) * (diff(pp$mean)^2)
varWithin <- function(pp) sum(pp$pro * pp$variance$sigmasq)
varTot <- function(pp) varBetween(pp) + varWithin(pp)
Rsquare <- function(pp) varBetween(pp) / varTot(pp)

getPowerScore <- function(pow) function(pp) Rsquare(pp) * (prod(pp$pro)/0.25)^pow

meansDistanceScore <- function(pp) diff(abs(pp$mean))

minProportion <- function(pp) min(pp$pro)

scoreFuns <- list(BI=BI.wang,
                  Mscore=Mscore,
                  R2=Rsquare,
                  mDS=meansDistanceScore)

## Compute bimodality scores for genes in the expression matrix 'X'.
## Each column is a gene.
## The function 'F' is a score formula, which can be one of (see scoreFuns):
## - BI.wang: 'Bimodality Index', by Wang et al.
## - Mscore: modified BI index, by Delorenzi
## - Rsquare: R square index
## - meansDistanceScore: absolute distance between mixture means
## Returns an score value for each gene
geneScores <- function(X, ...) {
  funs <- list(...)
  if(length(funs) == 0) {
    funs$F <- Rsquare
  }
  fits <- apply(X, 2, Compose(mixPar, mixFit))
  ans <- sapply(funs, function(F) sapply(fits, F))
  if(NCOL(ans)==1) {
    ans <- as.vector(ans)
    names(ans) <- colnames(X)
  } else {
    rownames(ans) <- colnames(X)
  }
  return(ans)
}

bimodalityScores <- function(X) {
  ans <- as.data.frame(geneScores(X, Rsquare=Rsquare,
                                  meansDistance=meansDistanceScore,
                                  proportion=minProportion))
  ans <- cbind(gene=rownames(ans), ans)
  ans <- transform(ans, gene=as.character(gene))
  return(ans)
}

bimodalityFilter <- function(X,
                             meansDistance=2,
                             proportion=0.01,
                             Rsquare=0) {
  mD <- meansDistance
  pr <- proportion
  rsq <- Rsquare
  with(bimodalityScores(X),
       gene[(meansDistance >= mD) & (proportion >= pr) & (Rsquare >= rsq)])
}

## remove uninteresting covariates effects from gene expression data
## X: gene expression matrix
## formula: rhs of a linear model formula
## data: dataframe containing the covariates specified in 'formula'
lmClean <- function(X, formula, data,
                    X.design=model.matrix(formula, data)) {
  means <- apply(X, 2, mean)
  E <- apply(X, 2, function(x) lm.fit(X.design, x)$residuals)
  return(sweep(E, 2, means, "+"))
}

## compute score 'F' for each stratum, then compute the average score
strataScores <- function(X, F, strata) {
  ans <- do.call(rbind, by(X, strata, function(Y) geneScores(Y, F),
                           simplify=FALSE))
  return(apply(ans, 2, mean))
}

## compute score 'F' for each data subset. Each time, stratify by 'strata'
subsetScores <- function(X, F, subset, strata) {
  ans <- list()
  for(k in unique(subset))
    ans[[k]] <- strataScores(X[subset==k,], F, strata[subset==k])
  return(do.call(cbind, ans))
}

mixDiscretize <- function(X) {
  if(is.matrix(X)) {
    return(apply(X, 2, .mixDiscretize))
  } else {
    return(.mixDiscretize(X))
  }
}
classifyBy <- mixDiscretize

## for each level of 'x', summarize df
mixTable <- function(x, df) {
  by(df, mixDiscretize(x), summary)
}

crossTable <- function(x, fac) {
  group <- mixDiscretize(x)
  return(table(fac, group))
}

## univariate mixture clustering for each column of X, separately for each stratification variable
mixDiscretizeStrata <- function(X, strata) {
  ans <- list()
  for(k in unique(strata))
    ans[[k]] <- mixDiscretize(X[strata==k,])
  ans <- do.call(rbind, ans)[rownames(X),]
  return(ans)
}
