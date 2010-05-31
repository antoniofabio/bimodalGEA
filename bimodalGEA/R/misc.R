##utility function: create a factor variable
## with levels sorted by frequency
factorSort <- function(fac) {
  cc <- sort(table(fac), decreasing=TRUE)
  return(factor(as.character(fac), levels=names(cc)))
}

.compose <- function(f,g) {
  f <- match.fun(f)
  g <- match.fun(g)
  function(x) f(g(x))
}
Compose <- function(...) Reduce(.compose, list(...))
