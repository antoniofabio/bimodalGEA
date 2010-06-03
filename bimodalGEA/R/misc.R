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

isUnique <- function(x) {
  all(!is.na(x)) && (length(unique(x)) == length(x))
}

Intersect <- function(lst) Reduce(intersect, lst)
Union <- function(lst) Reduce(union, lst)

outerColumns <- function(A, B=A, FUN) {
  FUN <- match.fun(FUN)
  if(is.vector(A)) {
    if(NROW(A) == NROW(B)) {
      ans <- outerColumns(A=matrix(A), B=B, FUN=FUN)
      return(structure(as.vector(ans), names=colnames(ans)))
    } else {
      stop("arguments 'A' and 'B' should be matrices")
    }
  }
  stopifnot(is.matrix(A) && is.matrix(B))
  ans <- t(apply(A, 2, function(a) {
    apply(B, 2, function(b) {
      FUN(a, b)
    })
  }))
  rownames(ans) <- colnames(A)
  colnames(ans) <- colnames(B)
  return(ans)
}
