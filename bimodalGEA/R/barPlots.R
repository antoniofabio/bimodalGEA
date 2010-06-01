crossBarPlot <- function(x, fac,
                         xlab="", legendTitle=NULL,
                         legend.text=TRUE,
                         args.legend=list(title=legendTitle),
                         counts=crossTable(x,fac),
                         beside=TRUE,
                         ...) {
  barplot(counts, beside=beside, xlab=xlab,
          legend.text=legend.text, args.legend=args.legend,
          ...)
}

##do a crossbarplot for each strata level
crossBarPlotStrata <- function(strata, x, fac, main, ...) {
  lvls <- levels(strata)
  for(lvl in lvls) {
    dev.new()
    xi <- x[strata == lvl]
    fi <- fac[strata == lvl]
    if(missing(main)) {
      crossBarPlot(xi, fi, main=lvl, ...)
    } else {
      crossBarPlot(xi, fi, ...)
    }
  }
}
