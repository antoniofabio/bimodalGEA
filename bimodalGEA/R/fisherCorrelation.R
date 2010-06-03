## evaluating association between categorical variables

fisherCorrelation <- function(x, y) fisher.test(table(x, y))$p.value

fisherCorrelations <- function(X, Y) {
  if(missing(Y))
    return(.fisherCorrelationMatrix(X))
  if(!is.matrix(X))
    return(apply(Y, 2, fisherCorrelation, x=X))
  return(t(apply(X, 2, fisherCorrelations, Y=Y)))
}

.fisherCorrelationMatrix <- function(X) {
  stopifnot(is.matrix(X))
  stopifnot(!is.null(colnames(X)))
  nms <- colnames(X)
  ans <- matrix(NA, length(nms), length(nms))
  colnames(ans) <- nms
  rownames(ans) <- nms
  for(i1 in seq(1, length(nms)-1))
    for(i2 in seq(i1+1, length(nms)))
      ans[i1,i2] <- ans[i2,i1] <- fisherCorrelation(X[, nms[i1]], X[, nms[i2]])
  return(ans)
}
