## Class for performing classical Rhodes MA method
## by deriving significance from p-values and permutating labels

##default substitute for pvalues of zero
default.zero.subst <- 1*10^-6

## perform method as proposed by Rhodes
## data: data frame or matrix containing p-values from various mehtods
## B: number of permutations
## zero.subst: value to replace zeros with
## cluster: snow cluster object to perform permutations with

doRhodesFDR <- function(data, B=10000, zero.subst= default.zero.subst, cluster=NULL){
  cat("starting method by Rhodes...\n")
  if(is.null(rownames(data))){
    stop("rownames must not be NULL")
  }
  if(is.null(cluster)){
    warning("you should use a cluster to get your results in bearable time!")
    #create a pseuo custe rwith one node only
    require(snow)
    cluster <- makeCluster(1, type = "MPI")
  }
   #substitute p-values of 0
  data[data == 0] <- zero.subst
  #calculate classical fisher sum S
  res <- data.frame(S=apply(data, 1, fishersum), row.names=rownames(data), p.value=rep(0, nrow(data)))
  #vector containing S p-values
  for(id in 1:nrow(data)){
    counter <- parSapply(cluster, 1:B, .doRandomS, data, res$S[id])   
    res$p.value[id] <- sum(counter)/B
    #output showing to be still alive
    if(id %% ceiling(nrow(data)/100) - 1 == 0){
	    cat("=")
	    flush.console()
    }
  }
  cat("\n")
  flush.console()
  #get corrected p-values
  res$q.value <- qvalue(res$p.value, lambda=0)$qvalues
  #rank result, solve ties by random assignment
  res$rank <- rank(res$q.value, ties="random")
  return(res)
}

## function for sampling a S-statistic
## and comparing it to the expected one
## b: permutation counter (is ignored)
## data: data containing studies p-values
## S: expected S
## Value: 1 if random S <= expected S, 0 else
.doRandomS <- function(b, data, S){
  rS <- sum(-2*log(apply(data, 2, sample, size=1)))
  if(rS >= S){
    cbind(1)
  } else {
    cbind(0)
  }
}
