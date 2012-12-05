plotPvsRank <- function(pvals, s.legend=NULL, colos=NULL, ylim=NULL, log="y", plot.title=NULL, zero.sub=0.00001, legend.pos="topright", lty=1, lwd=2){
  stopifnot(class(pvals)=="list")
  if(!all(sapply(pvals, function(xx){
    dim(xx)[2]==2
  })))
     stop("pvals must contain data.frames or arrays with two columns")
  if(is.null(names(pvals)) && is.null(s.legend))
    stop("You have to provide legend names for unnamed lists")
  if(is.null(names(pvals)))
    if(length(pvals) != length(s.legend))
      stop("Legend must have as many entries as the length of pvals")
  s.legend <- unlist(ifelse(is.null(names(pvals)), list(s.legend), list(names(pvals))))

  if(!is.null(ylim) && (any(is.na(ylim)) || length(ylim)!=2 || min(ylim)<0 || max(ylim)>1))
    stop("bad ylim provided")
  if(is.null(ylim)){
    ylim <- c(min(sapply(pvals, function(x) min(x[,2]))), 1)
    ylim <- sort(ylim, decreasing=TRUE)
  }
  
  if(length(pvals)>11 & is.null(colos))
    stop("You have to provide the color list yourself if there are more than 11 entries in pvals")
  if(!is.null(colos) && length(pvals)!=length(colos))
    stop("The length of provided colors does not fit the number of entries")
  if(is.null(colos)){
    nn <- max(length(pvals),3)
    colos <- RColorBrewer::brewer.pal(nn,"BrBG")
  }

  if(length(lty)!=length(pvals) && length(lty)!=1){
    stop("lty must by either 1 or as long as pvals")
  }
  if(length(lty)==1)
    lty <- rep(lty, length(pvals))
  
  if(length(lwd)!=length(pvals) && length(lwd)!=1){
    stop("lwd must by either 1 or as long as pvals")
  }
  if(length(lwd)==1)
    lwd <- rep(lwd, length(pvals))

  
  ##replace zeros with zero.sub
  pvals <- lapply(1:length(pvals), function(i){
    xx <- pvals[[i]]
    xx[xx<zero.sub] <- zero.sub
    xx
  })
  
  ## ymin <- min(sapply(pvals, min, na.rm=na.rm), na.rm=na.rm)
  ## ymax <- max(sapply(pvals, max, na.rm=na.rm), na.rm=na.rm)
  xmin <- 0
  xmax <- max(sapply(pvals, nrow))


  plot(NA, log=log, xlim=c(1, xmax), ylim=ylim, main=plot.title, xlab="index", ylab="p", las=1)
  for(i in 1:length(pvals)){
    temp <- pvals[[i]][order(pvals[[i]][,1]) ,]
    lines(rank(temp[,1], ties.method="first"), temp[,2], type="l", lty=lty[i], lwd=lwd[i], col=colos[i])
  }
  grid()
  legend(legend.pos, legend=s.legend, bty="n", lwd=lwd, lty=lty, col=colos)
}
