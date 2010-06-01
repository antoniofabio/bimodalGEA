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
