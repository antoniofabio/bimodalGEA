##utility function: create a factor variable
## with levels sorted by frequency
factorSort <- function(fac) {
  cc <- sort(table(fac), decreasing=TRUE)
  return(factor(as.character(fac), levels=names(cc)))
}

.compose <- function(f,g) {
  force(f)
  force(g)
  function(x) f(g(x))
}
Compose <- function(...) Reduce(.compose, list(...))
