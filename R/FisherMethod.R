## Class for performing classical Fisher method
## by deriving significance from chi square distribution

## perform classical fisher method
## rownames must not be NULL
## data: data.frame with features in rows and studies in columns
## zero.subst: since p-values of 0 cause problems they are substituted
## by a very small number (default is default.zero.subst)
fisherMethod <- function(data, zero.subst= default.zero.subst,cluster=NULL){
	cat("starting Fishers method...\n")
	if(is.null(rownames(data))){
		stop("rownames must not be NULL")
	}
	#substitute p-values of 0
	data[data == 0] <- zero.subst
	#calculate classical fisher sum S
	res <- data.frame(S = apply(data, 1, fishersum), row.names = rownames(data))
	#derive p-values from chi square
	res$p.value <- 1-pchisq(res$S, df=(ncol(data)*2))
	#get corrected p-values
	res$q.value <- qvalue(res$p.value, lambda=0)$qvalues
	#rank result, solve ties by random assignment
	res$rank <- rank(res$q.value, ties="random")
	return(res)
}

## function to calculate fisher sum
## P: vector of p-values
fishersum <- function(P){
	return(sum(-2*log(P)))
}