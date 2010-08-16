## Class to wrap parameters for Effect Size Meta Analysis

## function hiding callof ES method by Choi
## A: list o ExpressionSets
## cl: list classes to be tested for first class 
## should be in the first sample of the first study)
## useREM: boolean indictaing wether to use REM ES model
## nperm: number of permutations to calculate FDR (default is 5000)

doES <- function(A, cl, cl.val= c(0,1), useREM=TRUE, nperm= 5000, cluster = NULL){
#  if (is.null(cluster)) {
#	  warning("cluster should not be NULL, computing time may be long")
#  }	
	if(!all(levels(factor(cl.val)) == levels(factor(unlist(cl))))){
		stop("factor levels not the same as cl.val")
	}
	if(length(levels(unlist(cl))) != 2){
		stop("ES wrapper works with two factors only")
	}
	if(length(A) == 0){
		stop("called ES with an empty list of ExpressionSets")
	}
	if(cl[[1]][1] != cl.val[1]){
		stop(paste("first sample in list of ExpressionSets has to be of class", cl.val[1]))
	}
	res <- list()
	res.temp <- zScoreFDRClust(A, cl, useREM=useREM, nperm=nperm, cluster = cluster)
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

####################################################
####### modified methods with optional clusters: ###
####################################################

## zScoreFDR-method with cluster-option
zScoreFDRClust = function (esets, classes, useREM = TRUE, nperm = 1000, CombineExp = 1:length(esets), cluster = NULL) 
{
	for (i in 1:length(classes)) {
		if (!is.factor(classes[[i]])) {
			classes[[i]] <- factor(classes[[i]])
		}
		if (nlevels(classes[[i]]) != 2) {
			stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
		}
		else {
			Ref <- levels(classes[[i]])[1] 
			## cluster
			if (is.null(cluster)) {
				classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == Ref, 0, 1))
				
			} else{
				classes[[i]] <- parSapply(cluster,classes[[i]], function(x) ifelse(x == Ref, 0, 1))
				
			}
		}
	}
	num.studies <- length(esets)
	num.genes <- nrow(exprs(esets[[1]]))
	zscoresAll <- zScoresClust(esets, classes, useREM = useREM, CombineExp = CombineExp, cluster = cluster)
	MuQu <- zscoresAll[, c("MUvals", "MUsds", "Qvals", "df", 
					"Qpvalues", "Chisq")]
	zscore <- zscoresAll[, c(paste("zSco_Ex_", 1:num.studies, 
							sep = ""), "zSco")]
	aperms <- .replicateClust(nperm, zScores(esets, lapply(classes, 
							sample), useREM, CombineExp = CombineExp)[, c(paste("zSco_Ex_", 
									1:num.studies, sep = ""), "zSco")], simplify = FALSE)
	k <- c("pos", "neg", "two.sided")
	All <- vector(mode = "list", length = 3)
	for (l in 1:length(k)) {
		theFDR <- .multExpFDRClust(zscore, aperms, type = k[l],cluster=cluster)
		n <- num.studies + 1
		i <- 1:n
		theResult <- matrix(NA, ncol = n * 2, nrow = nrow(zscore))
		rownames(theResult) <- rownames(zscore)
		stopifnot(rownames(theResult) == rownames(theFDR))
		theResult[, 2 * i - 1] <- zscore
		theResult[, 2 * i] <- theFDR
		colnames(theResult) <- 1:(2 * n)
		colnames(theResult)[2 * i - 1] <- paste("zSco_Ex_", i, 
				sep = "")
		colnames(theResult)[2 * i] <- paste("FDR_Ex_", i, sep = "")
		colnames(theResult)[(2 * n - 1):(2 * n)] <- c("zSco", 
				"FDR")
		All[[l]] <- cbind(theResult, MuQu)
	}
	names(All) <- k
	return(All)
}

## zScore-method with cluster-option
zScoresClust = function (esets, classes, useREM = TRUE, CombineExp = 1:length(esets), cluster = NULL) 
{
	num.studies <- length(esets)
	num.genes <- nrow(exprs(esets[[1]]))
	if (num.studies != length(classes)) 
		stop("Length of classes must be equal to length of esets.")
	for (i in 1:num.studies) {
		if (!is.factor(classes[[i]])) {
			classes[[i]] <- factor(classes[[i]])
		}
		if (nlevels(classes[[i]]) != 2) {
			stop("Error: Each list in the argument \"classes\" must contain only 2 levels.")
		}
		else {
			Ref <- levels(classes[[i]])[1]
			
			## cluster
			if (is.null(cluster)) {			
				classes[[i]] <- sapply(classes[[i]], function(x) ifelse(x == Ref, 0, 1))
				
			} else{
				
				classes[[i]] <- parSapply(cluster,classes[[i]], function(x) ifelse(x == Ref, 0, 1))				
			}
		}
	}
	tau2 <- function(Q, num.studies, my.weights) {
		vwts <- rowSums(my.weights)
		tmp2 <- rowSums(my.weights^2)
		tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
		return(tau2)
	}
	theNames <- featureNames(esets[[1]])
	for (i in 2:num.studies) stopifnot(identical(theNames, featureNames(esets[[i]])))
	ds <- matrix(NA, ncol = num.studies, nrow = num.genes)
	vars <- matrix(NA, ncol = num.studies, nrow = num.genes)
	for (i in 1:length(esets)) {
		my.d.adj <- dstar(getdF(esets[[i]], classes[[i]]), length(classes[[i]]))
		ds[, i] <- as.numeric(my.d.adj)
		vars[, i] <- as.numeric(sigmad(my.d.adj, sum(classes[[i]] == 
										0), sum(classes[[i]] == 1)))
	}
	sepZscores <- ds/sqrt(vars)
	effects <- ds
	effectsVar <- vars
	colnames(sepZscores) <- paste("zSco_Ex_", 1:num.studies, 
			sep = "")
	colnames(effects) <- paste("Effect_Ex_", 1:num.studies, sep = "")
	colnames(effectsVar) <- paste("EffectVar_Ex_", 1:num.studies, 
			sep = "")
	ds <- ds[, CombineExp]
	vars <- vars[, CombineExp]
	num.studies <- length(CombineExp)
	df <- num.studies - 1
	Qvals <- f.Q(ds, vars)
	if (useREM) 
		vars <- vars + tau2(Qvals, num.studies, my.weights = 1/vars)
	wt <- 1/vars
	MUvals <- rowSums(ds * wt)/rowSums(wt)
	MUsds <- sqrt(1/rowSums(wt))
	zSco <- MUvals/MUsds
	Qpvalues <- 1 - pchisq(Qvals, df)
	Chisq <- 1 - pchisq(zSco^2, 1)
	theResult <- cbind(sepZscores, zSco, MUvals, MUsds, Qvals, 
			df, Qpvalues, Chisq, effects, effectsVar)
	rownames(theResult) <- theNames
	return(theResult)
}

## replicate-method with cluster-option
.replicateClust = function (n, expr, simplify = TRUE, cluster = NULL) {
	if (is.null(cluster)) {		
		res = sapply(integer(n), eval.parent(substitute(function(...) expr)), 
				simplify = simplify)
		
	} else{
		
		res = parSapply(cluster, integer(n), eval.parent(substitute(function(...) expr)), 
				simplify = simplify)
		
	}
	return(res)
}

## multExpFDR-method with cluster-option
.multExpFDRClust = function (theScores, thePermScores, type = "pos", cluster = NULL) 
{
	numberOfPermutations <- length(thePermScores)
	if (!type %in% c("two.sided", "pos", "neg")) 
		stop("Wrong type!")
	ff <- function(x) -x
	biggerEq <- function(x, y) {
		y <- sort(y, decreasing = TRUE)
		a <- match(x, x)
		b <- x %in% y
		d <- match(x, sort(c(x, y), decreasing = TRUE))
		return(d - a + b)
	}
	theFDR <- matrix(NA, nrow = nrow(theScores), ncol = ncol(theScores))
	rownames(theFDR) <- rownames(theScores)
	if (type == "two.sided") {
		theScores <- abs(theScores)
		
		## cluster
		if (is.null(cluster)) {
			thePermScores <- lapply(thePermScores, abs)
			
		} else{
			
			thePermScores <- parLapply(cluster,thePermScores, abs)				
		}
		
	}
	if (type == "neg") {
		theScores <- -theScores
		
		## cluster
		if (is.null(cluster)) {			
			thePermScores <- lapply(thePermScores, ff)
			
		} else{
			
			thePermScores <- parLapply(cluster,thePermScores, ff)				
		}
	}
	for (i in 1:ncol(theScores)) {
		ord <- order(theScores[, i], decreasing = TRUE)
		z <- theScores[ord, i]
		
		## cluster
		if (is.null(cluster)) {
			
			randomZ <- as.vector(sapply(thePermScores, function(x) x[,i]))
			
		} else{
			
			randomZ <- as.vector(parSapply(cluster,thePermScores, function(x) x[,i]))				
		}
		
		randomZ <- sort(randomZ, decreasing = TRUE)
		numberisBigger <- biggerEq(z, randomZ)
		theFDR[ord, i] <- numberisBigger/((1:length(z)) * numberOfPermutations)
	}
	return(theFDR)
}