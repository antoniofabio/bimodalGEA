## evaluating association between categorical variables

fisherCorrelation <- function(xa, xb) fisher.test(table(xa, xb))$p.value

fisherCorrelations <- function(X.discrete) {
  nms <- colnames(X.discrete)
  ans <- matrix(NA, length(nms), length(nms))
  colnames(ans) <- nms
  rownames(ans) <- nms
  for(i1 in seq(1, length(nms)-1))
    for(i2 in seq(i1+1, length(nms)))
      ans[i1,i2] <- ans[i2,i1] <- fisherCorrelation(X.discrete[, nms[i1]], X.discrete[, nms[i2]])
  return(ans)
}
