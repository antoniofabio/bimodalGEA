## evaluating association between categorical variables

fisherCorrelation <- function(x, y) fisher.test(table(x, y))$p.value

fisherCorrelations <- function(X, Y) {
  args <- list(A=X, FUN=fisherCorrelation)
  if(!missing(Y)) {
    args$B <- Y
  }
  return(do.call(outerColumns, args))
}
