library(genefilter)

## Class to wrap parameters for Effect Size Meta Analysis

## function hiding callof ES method by Choi
## A: list o ExpressionSets
## cl: list classes to be tested for first class 
## should be in the first sample of the first study)
## useREM: boolean indictaing wether to use REM ES model
## nperm: number of permutations to calculate FDR (default is 5000)
doES <- function(A, cl, val= c(0,1), useREM=TRUE, nperm= 5000){
  if(!all(levels(factor(val)) == levels(factor(unlist(cl))))){
    stop("factor levels not the same as val")
  }
  if(length(levels(unlist(cl))) != 2){
    stop("ES wrapper works with two factors only")
  }
  if(length(A) == 0){
    stop("called ES with an empty list of ExpressionSets")
  }
  if(cl[[1]][1] != val[1]){
    stop(paste("first sample in list of ExpressionSets has to be of class", val[1]))
  }
  res <- list()
  res.temp <- zScoreFDR(A, cl, useREM=useREM, nperm=nperm)
  res$g2up <- data.frame(res.temp$neg)
  res$g2down <- data.frame(res.temp$pos)

  #add rank information
  res$g2up$rank <- rank(res$g2up$FDR, ties="random")
  res$g2down$rank <- rank(res$g2down$FDR, ties="random")

  #truncate FDR to be between [0,1]
  res$g2down$FDR[res$g2down$FDR > 1] <- 1
  res$g2up$FDR[res$g2up$FDR > 1] <- 1
  #report FDR as p-value for this method (will be resued in the FM ensemble)
  res$g2up$p.value <- res$g2up$FDR
  res$g2up$p.value <- res$g2up$FDR


  return(res)
}
