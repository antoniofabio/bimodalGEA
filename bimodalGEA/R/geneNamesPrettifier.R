## Antonio, Fabio Di Narzo, 2010
##
## Make duplicated (and empty) gene symbols unique
## by adding the geneId identifier when necessary
##
geneNamesPrettifier <- function(geneSymbols, geneIds) {
  ans <- geneSymbols
  isEmpty <- ans==""
  ans[isEmpty] <- geneIds[isEmpty]
  smb <- table(ans)
  duplId <- which(ans %in% names(smb)[smb>1])
  ans[duplId] <- sprintf("%s (%s)", ans[duplId],
                         as.character(geneIds[duplId]))
  return(ans)
}

disambiguate <- function(x, disambiguator=seq_along, sep=" ") {
  stopifnot(is.vector(x))
  stopifnot(all(!is.na(x)))
  disambiguator <- match.fun(disambiguator)
  x <- as.character(x)
  tbl <- table(x)
  duplicated <- names(tbl)[tbl>1]
  for(nm in duplicated) {
    inm <- x==nm
    x[inm] <- paste(x[inm], disambiguator(x[inm]), sep=sep)
  }
  return(x)
}
