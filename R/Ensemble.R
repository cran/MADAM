## Class provdiding infterface for Ensemble Approach
## combining classical Fisher method, effect size (Choi) and Rank Prod approach

## TODO: make single methods optional according to paramter settings
## TODO: include Rhodes method

##method starting the Ensemble Approach
## A: list of ExpressionSets, sharing the same identifiers
## and the same ordering of them
## cl: list of factors, for discriminating between the two groups
## val: vector of objects indicating group1 and group2 (default: 0, 1)
## feature.names: names of features common to all sets. 
## if NULL featureNames(A[[1]]) is used. 
## Order has to be the same acros the sets
## do.FM: do meta analysis based on significance (Fisher method)
## perm.ES: number of permutation for calculating FDR for ES Choi
## perm.RP: number of permutaionts for rank product approach
## perm.ENS.RP: number of permutations for Ensemble rank product approach
## perm.ENS.ES: number of oermutations for Ensemble rank sum approach
## do.ENS.FM: do Ensemble Fisher Method
## cluster: snow object (for parallel computing)
## useREM: boolean parameter indicating wether to use ES REM model (default is TRUE)
## results.dir: directory to which results are written as csv. 
## The script will create a subfolde with date and time, where everything will be stored
## write.all: if TRUE all results are written as csv a subfolder of results.dir
## plot.fdr: plot FDR figures for MA & Ensemble methods
## plot.ranks: plot heatmaps for single MA methods
##
## Value: list object with all the results
## TODO: write value details
doEnsemble <- function(A, cl, cl.val=c(0,1), feature.names=NULL, do.FM= TRUE, perm.ES=1000, perm.RP=1000, perm.ENS.RP= 1000, perm.ENS.RS= 1000, do.ENS.FM= TRUE, cluster= NULL, useREM= TRUE, results.dir=getwd(), write.all=FALSE, plot.fdr= FALSE, plot.ranks=FALSE){
	require(Biobase)
	#validy checks
	if(!all(levels(factor(cl.val)) == levels(factor(unlist(cl))))){
		stop("flist factor levels not the same as cl.val!")
	}
	if(length(levels(unlist(cl))) != 2){
		stop("Ensemble Approach works with two factors only!")
	}
	if(perm.ES < 1){
		stop("number of permutations has to be at least 1 for ES approach!")
	}
	if(perm.RP < 1){
		stop("number of permutations has to be at least 1 for RP approach!")
	}
	if(length(A) == 0){
		stop("called method with an empty list of ExpressionSets!")
	}
	if(cl[[1]][1] != cl.val[1]){
		stop(paste("first sample in list of ExpressionSets has to be of class", cl.val[1], "!"))
	}
	if(is.null(feature.names)){
		feature.names= featureNames(A[[1]])
	}
	if(write.all){
		results.dir <- paste(results.dir, format(Sys.time(), "%Y%m%d-%H%M%S"), sep="/")
		dir.create(results.dir)
	}
	#check for featurenames to be consistent across sets of A
	if (is.null(cluster)) {
		warning("cluster should not be NULL, computing time may be long")
		
		if (!all(sapply(A, function(a) {
							identical(featureNames(a), feature.names)
						}))) {
			stop("featureNames are not consistent accros sets!")
		}
		
	} else{
		
		if (!all(parSapply(cluster,A, function(a) {
							identical(featureNames(a), feature.names)
						}))) {
			stop("featureNames are not consistent accros sets!")
		}
		
	}
	
	cat("Starting Meta Analysis Ensemble Method...\n")
	
	#list for storing results
	res <- list()
	
	#perform multi t-test in both direction
	if(do.FM){
		cat("perform multi t-test in both direction...\n")
		res$ttest.g2up <- multiTtest(A, cl, cl.val=cl.val, alternative="less",cluster=cluster)
		res$ttest.g2down <- multiTtest(A, cl, cl.val=cl.val, alternative="greater",cluster=cluster)
	}
	
	#perform fisher method for both directions
	if(do.FM){
		cat("perform fisher method for both directions...\n")
		res$MAFM.g2up <- fisherMethod(res$ttest.g2up,cluster=cluster)
		res$MAFM.g2down <- fisherMethod(res$ttest.g2down,cluster=cluster)
	}
	
	#perform ES (Choi)
	if(perm.ES > 0){
		cat("performing effect size method (Choi)...\n")
		res.es <- doES(A, cl, cl.val, useREM=useREM, nperm=perm.ES,cluster=cluster)
		res$MAES.g2up <- res.es$g2up
		res$MAES.g2down <- res.es$g2down
		res.es <- NULL
	}
	
	#perform RP
	if(perm.RP > 0){
		cat("performing rank product method...\n")
		res.rp <- doRP(A, cl, cl.val, nperm=perm.RP, gene.names=feature.names,cluster=cluster)
		res$MARP.g2up <- res.rp$g2up
		res$MARP.g2down <- res.rp$g2down
		res.rp <- NULL
	}
	
	if(do.FM && perm.ES > 0 && perm.RP > 0){
		# make list of ranks for each group
		cat("creating MA rank lists...\n")
		res$ranks.g2up <- data.frame(MAFM=res$MAFM.g2up$rank[match(rownames(res$MAFM.g2up), feature.names)])
		res$ranks.g2up$MAES <- res$MAES.g2up$rank[match(rownames(res$MAES.g2up), feature.names)]
		res$ranks.g2up$MARP <- res$MARP.g2up$rank[match(rownames(res$MARP.g2up), feature.names)]
		rownames(res$ranks.g2up) <- feature.names
		res$ranks.g2down <- data.frame(MAFM=res$MAFM.g2down$rank[match(rownames(res$MAFM.g2down), feature.names)])
		res$ranks.g2down$MAES <- res$MAES.g2down$rank[match(rownames(res$MAES.g2down), feature.names)]
		res$ranks.g2down$MARP <- res$MARP.g2down$rank[match(rownames(res$MARP.g2down), feature.names)]
		rownames(res$ranks.g2down) <- feature.names
		
		#make list of raw p-values for each method
		cat("creating MA p-value lists...\n")
		res$pvalues.g2up <- data.frame(MAFM=res$MAFM.g2up$p.value[match(rownames(res$MAFM.g2up), feature.names)])
		res$pvalues.g2up$MAES <- res$MAES.g2up$p.value[match(rownames(res$MAES.g2up), feature.names)]
		res$pvalues.g2up$MARP <- res$MARP.g2up$p.value[match(rownames(res$MARP.g2up), feature.names)]
		rownames(res$pvalues.g2up) <- feature.names
		res$pvalues.g2down <- data.frame(MAFM=res$MAFM.g2down$p.value[match(rownames(res$MAFM.g2down), feature.names)])
		res$pvalues.g2down$MAES <- res$MAES.g2down$p.value[match(rownames(res$MAES.g2down), feature.names)]
		res$pvalues.g2down$MARP <- res$MARP.g2down$p.value[match(rownames(res$MARP.g2down), feature.names)]
		rownames(res$pvalues.g2down) <- feature.names
		
		#calculate correlations between the methods based on ranks 
		cat("calculating correlations between MA methods based on ranks...\n")
		res$corr.m.g2up <- cor(res$ranks.g2up)
		res$corr.m.g2down <- cor(res$ranks.g2down)
		
		cat("starting Ensemble Approaches...\n")
		
		if(perm.ENS.RP > 0){
			#RP as ensemble
			res$ENSRP.g2up <- calculateRankProduct(res$ranks.g2up, B=perm.ENS.RP, cluster=cluster)
			res$ENSRP.g2down <- calculateRankProduct(res$ranks.g2down, B=perm.ENS.RP, cluster=cluster)
		}
		
		if(perm.ENS.RS > 0){
			#RS as ensemble
			res$ENSRS.g2up <- calculateRankSum(res$ranks.g2up, B=perm.ENS.RS, cluster=cluster)
			res$ENSRS.g2down <- calculateRankSum(res$ranks.g2down, B=perm.ENS.RS, cluster=cluster)
		}
		
		if(do.ENS.FM){
			#fisher as ensemble
			res$ENSFM.g2up <- fisherMethod(res$pvalues.g2up)
			res$ENSFM.g2down <- fisherMethod(res$pvalues.g2down)
		}
		
		#put all ensemble q-values together
		cat("collecting ensemble q-values...\n")
		res$ENSQ.g2up <- data.frame(ENSRP=res$ENSRP.g2up$q.value, ENSRS=res$ENSRS.g2up$q.value, ENSFM=res$ENSFM.g2up$q.value, row.names=feature.names)
		res$ENSQ.g2down <- data.frame(ENSRP=res$ENSRP.g2down$q.value, ENSRS=res$ENSRS.g2down$q.value, ENSFM=res$ENSFM.g2down$q.value, row.names=feature.names)
		
		#put all single MA Q-values together
		cat("collecting MA q-values...\n")
		res$MAQ.g2up <- data.frame(MAFM=res$MAFM.g2up$q.value, MAES=res$MAES.g2up$FDR, MARP=res$MARP.g2up$q.value, row.names=feature.names)
		res$MAQ.g2down <- data.frame(MAFM=res$MAFM.g2down$q.value, MAES=res$MAES.g2down$FDR, MARP=res$MARP.g2down$q.value, row.names=feature.names)
		
		#put all q-values together
		cat("colliding q-values...\n")
		res$ALLQ.g2up <- data.frame(cbind(res$MAQ.g2up, res$ENSQ.g2up), row.names=feature.names)
		res$ALLQ.g2down <- data.frame(cbind(res$MAQ.g2down, res$ENSQ.g2down), row.names=feature.names)
		
		#if user set flag, plot FDR information from every MA & ENS method
		if(plot.fdr){
			cat("plotting FDR figures...\n")
			ty.up <- factor(c(rep("m", ncol(res$MAQ.g2up)), rep("e", ncol(res$ENSQ.g2up))))
			tiff(paste(results.dir, "FDR_up.tif", sep="/"))
			plotFDR(res$ALLQ.g2up, types=ty.up, title="Significances (up)")
			dev.off()
			ty.down <- factor(c(rep("m", ncol(res$MAQ.g2down)), rep("e", ncol(res$ENSQ.g2down))))
			tiff(paste(results.dir, "FDR_down.tif", sep="/"))
			plotFDR(res$ALLQ.g2down, types=ty.down, title="Significances (down)")
			dev.off()
		}
		
		#if user set flag, plot FDR information from every MA & ENS method
		if(plot.ranks){
			cat("plotting rank heatmaps...\n")
			tiff(paste(results.dir, "ranks_up.tif", sep="/"))
			plotHeatMap(as.matrix(res$ranks.g2up), title="Ranks (up)")
			dev.off()
			tiff(paste(results.dir, "ranks_down.tif", sep="/"))
			plotHeatMap(as.matrix(res$ranks.g2down), title="Ranks (down)")
			dev.off()
			
		}
	}
	
	#put everything into result list
	res$feature.names <- feature.names
	res$results.dir <- results.dir
	
	#write everything that is in res to the results.dir
	if (write.all) {
		
		if (is.null(cluster)) {
			
			lapply(1:length(res), function(rr) {
						write.csv(res[[rr]], paste(results.dir, paste(names(res)[rr], 
												"csv", sep = "."), sep = "/"), row.names = TRUE)
						cat("writing ", paste(results.dir, names(res)[rr], 
										sep = "/"), ".csv...\n", sep = "")
					})
			
		} else{
			
			parLapply(cluster,1:length(res), function(rr) {
						write.csv(res[[rr]], paste(results.dir, paste(names(res)[rr], 
												"csv", sep = "."), sep = "/"), row.names = TRUE)
						cat("writing ", paste(results.dir, names(res)[rr], 
										sep = "/"), ".csv...\n", sep = "")
					})
			
		}
		
	}
	
	cat("Successfully finished Meta Analysis Ensemble Method...\n")
	return(res)
	
}
