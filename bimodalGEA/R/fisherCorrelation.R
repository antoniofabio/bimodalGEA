## evaluating association between categorical variables

fisherCorrelation <- function(x, y) {
  tryCatch({
    fisher.test(table(x, y))$p.value
  }, error=function(e) NA)
}

fisherCorrelations <- function(X, Y) {
  args <- list(A=X, FUN=fisherCorrelation)
  if(!missing(Y)) {
    args$B <- Y
  }
  return(do.call(outerColumns, args))
}
