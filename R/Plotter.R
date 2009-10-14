## Class for plotting ranklists and FDR images

## function for plotting ranklist heatmap
## rankList: data.frame containing ranks, rows features columns study, methods...
## title: image title
## top.col: vector of colors used for plotting column headers (default= NULL)
## side.col: vector indicating how to print row side colors (default= NULL)
## col: colors to be used for plotting ranks (default: rainbow(nrow(rankList))
plotHeatMap <- function(rankList, title="", top.col=NULL, side.col=NULL, col=rainbow(nrow(rankList))){
  cat("plotting heatmaps\n")
  if(is.null(top.col) && !is.null(side.col)){
  heatmap.2(as.matrix(rankList), distfun = function(m) dist(m, method="euclidian"), scale="none", col=col, RowSideColors=side.col, trace="none", main=title)
  return()
  }
  if(is.null(side.col) && !is.null(top.col)){
  heatmap.2(as.matrix(rankList), distfun = function(m) dist(m, method="euclidian"), scale="none", col=col, ColSideColors=top.col, trace="none", main=title)
  return()
  }
  if(is.null(side.col) && is.null(top.col)){
  heatmap.2(as.matrix(rankList), distfun = function(m) dist(m, method="euclidian"), scale="none", col=col, trace="none", main=title)
  return()
  }
  if(!is.null(side.col) && !is.null(top.col)){
  heatmap.2(as.matrix(rankList), distfun = function(m) dist(m, method="euclidian"), scale="none", col=col, ColSideColors=top.col, RowSideColors=side.col, trace="none", main=title)
  }	
}


##plot results FDRs in one image
## x: matrix, rownames beig features and colnames being studies/methods
## title: title for image
## types: 2 level factor indicating if a x comes from a normal 
##method or a ensemble approach, difference will be in color and line type
## allowed levels are "s", "m", and "e"
plotFDR <- function(x, title="", types=factor(rep("s", ncol(x)))){
  cat("plotting FDR...\n")
  allowed.levels <- c("s","m","e")
  if(!all(levels(types) %in% allowed.levels)){
    stop(paste("types have to be in", allowed.levels))
  }

  #prepare colors + line types
  sum.s <- sum(types == "s")
  sum.m <- sum(types == "m")
  sum.e <- sum(types == "e")
  
  if(sum.s < 3){
    cols.s <- rainbow(sum.s)
  } else {
    cols.s <- brewer.pal(sum.s, "Set1")
  }
  if(sum.m < 3){
    cols.m <- heat.colors(sum.m)
  } else {
    cols.m <- brewer.pal(sum.m, "Set2")
  }
  if(sum.e < 3){
    cols.e <- topo.colors(sum.e)
  } else {
    cols.e <- brewer.pal(sum.e, "Set3")
  }

  cols <- rep("", length(types))
  ltys <- rep("", length(types))

  cols[which(types== "s")] <- cols.s
  cols[which(types== "m")] <- cols.m
  cols[which(types== "e")] <- cols.e

  ltys[which(types== "s")] <- 1
  ltys[which(types== "m")] <- 2
  ltys[which(types== "e")] <- 3

  for(m in 1:ncol(x)){ 
    data <- x[rank(order(x[,m]), ties="random"),m]
    if(m==1){
      plot(rank(data, ties="first"), data, type="l", lty=as.numeric(ltys[m]), col=cols[m], ylim=c(1,0), ylab="Significance", xlab="Rank", main=title)
    } else {
      lines(rank(data, ties="first"), data, type="l", lty=as.numeric(ltys[m]), col=cols[m])
    }
  }
  if(is.null(colnames(x))){
    coln <- 1:ncol(x)
  } else {
    coln <- colnames(x)
  }
  legend("topright", as.character(coln), lty=as.numeric(ltys), col=cols)
}

## function to plot a volcano plot using information from Effect Size 
## based MA and significance based MA
## todo: provide points to be plotted in there as well
## x: matrix with two columns (first column containing ES 
## title: title for plot
## col: color for volcano plot
## pch: character type for volcano plot
## points: matrix with two columns: specifing if certain points are to be drawn separately
## points.pch: character type for plotting points
## points.col: color for plotting points
## and second significance information
plotMAVolcano <- function(x, title="", col="black", pch=".", points=NULL, points.col="red", points.pch=1){
  plot(x[,1], x[,2], ylim=c(max(x[,2]), min(x[,2])), col=col, pch=pch, xlab="Effect Size", ylab="Significance", main=title)
  if(!is.null(points)){
    points(points, pch=points.pch, col=points.col)
  }
}
