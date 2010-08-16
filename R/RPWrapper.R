## Class to wrap parameters for Rank Product Meta Analysis

## function hiding callof RankProduct method by Breitling
## A: list of ExpressionSets
## cl: list classes to be tested for first class 
## should be in the first sample of the first study)
## nperm: number of permutations to calculate FDR (default is 1000)
doRP <- function(A, cl, cl.val= c(0,1), nperm= 1000, gene.names=NULL, cluster=NULL){
	if(require(RankProd)){
		if(!all(levels(factor(cl.val)) == levels(factor(unlist(cl))))){
			stop("factor levels not the same as cl.val")
		}
		if(length(levels(unlist(cl))) != 2){
			stop("multiTTest works with two factors only")
		}
		if(length(A) == 0){
			stop("called ES with an empty list of ExpressionSets")
		}
		if(cl[[1]][1] != cl.val[1]){
			stop(paste("first sample in list of ExpressionSets has to be of class", cl.val[1]))
		}
		if(is.null(gene.names)){
			gene.names= featureNames(A[[1]])
		}
		
		#convert ExpressionSets to exprs matrix
		if (is.null(cluster)) {
			temp.rp <- lapply(A, exprs)
		} else{
			temp.rp <- parLapply(cluster,A, exprs)
		}
		#container for data
		rp.data <- c()
		#container for origin of data study 1 = 1, study 2= 2, etc.
		rp.orig <- c()
		for(i in 1:length(temp.rp)){
			rp.data <- cbind(rp.data, temp.rp[[i]])
			rp.orig <- c(rp.orig, rep(i, ncol(temp.rp[[i]])))
		}
		#convert cl to numeric
		if (is.null(cluster)) {
			cl.num <- lapply(cl, function(cc) {as.numeric(cc) - 1})
		} else{
			cl.num <- parLapply(cluster,cl, function(cc) {as.numeric(cc) - 1})
		}
		
		res.rp <- RPadvance(rp.data, unlist(cl.num), rp.orig, num.perm=nperm, gene.names=gene.names)
		
		#get rank information from method
		res.rp.ranked <- res.rp$RPrank
		colnames(res.rp.ranked)[grep("<", colnames(res.rp.ranked))] <- "g2up"
		colnames(res.rp.ranked)[grep(">", colnames(res.rp.ranked))] <- "g2down"
		#get unadjusted p-values
		res.rp.sig <- res.rp$pval
		colnames(res.rp.sig)[grep("<", colnames(res.rp.sig))] <- "g2up"
		colnames(res.rp.sig)[grep(">", colnames(res.rp.sig))] <- "g2down"
		#get FDR info
		res.rp.pfp <- res.rp$pval
		colnames(res.rp.pfp)[grep("<", colnames(res.rp.pfp))] <- "g2up"
		colnames(res.rp.pfp)[grep(">", colnames(res.rp.pfp))] <- "g2down"
		#crop fdr bigger than 1 to 1
		res.rp.pfp[res.rp.pfp > 1] <- 1
		res <- list()
		res$g2up <- data.frame(p.value= res.rp.sig$g2up, q.value=res.rp.pfp$g2up, rank=res.rp.ranked$g2up, row.names=gene.names)
		res$g2down <- data.frame(p.value= res.rp.sig$g2down, q.value=res.rp.pfp$g2down, rank=res.rp.ranked$g2down, row.names=gene.names)
		
		return(res)
	}
}
