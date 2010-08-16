##data generation for meta analysis validation

## function to sample sample sizes in studies
## if no t provided
## size: number of studies, default= 1
.getNbSamples <- function(size=1){
	k1 <- sample(c(runif(6*10*size, 8, 15), runif(1*10*size, 5, 14), runif(3*10*size, 15, 30)), size)
	return(floor(k1))
}

## function to sample a number of studies, used for validation of meta analysis methods
## g: number of genes
## perc.sig: percentage of significant genes
## i: number of studies
## k: number of samples for each of the two groups/study
## err.span: span for var in error

## Value: a list of ExpressionSets, with a nassigned group label "group1" or "group2"
generateRandomMAData <- function(g=5000, perc.sig=0.1, i=2, k=0, err.span= c(0.1, 0.2)){
	if(k == 0 && i != 0){
		k = .getNbSamples(i*2)
	}
	
	#ids of significant genes
	id.sig <- sample(seq(1:g), g*perc.sig, replace=F)
	
	TAU <- matrix(0, nrow=g, ncol=i)
	DELTA <- runif(i, 0, 2)
	DELTA.sigma <- runif(i, 0.4, 0.5)
	
	B = matrix(0, nrow=g, ncol=i)
	MU <- sample(c(runif(2*g, 0.2, 1), runif(8*g, 1, 2.2)),g, replace=F)
	#determin wether MU is up or down 
	MU.negpos <- rbinom(length(MU), 1, 0.5)
	MU.sigma <- runif(g, 0.1, 0.5)
	
	
	for(ii in 1:ncol(TAU)){
		TAU[,ii] <- rnorm(g, DELTA[ii], DELTA.sigma[ii])
		B[,ii] <- abs(rnorm(g, MU[ii], MU.sigma[ii]))*(-1)^MU.negpos
		B[-id.sig,ii] <- 0
	}
	
	
	
	a = sample(c(runif(9*g, 4.5, 9), runif(1*g, 6, 12)),g, replace=F)
	A <- vector("list",i)
	ES <- vector("list", i)
	for(x in 1:i){
		k.study <- k[2*x-1] + k[2*x]
		A[[x]] <- matrix(rep(a, k.study), ncol=sum(k.study), nrow=g, byrow=F)
		A[[x]] <- A[[x]] + TAU[,x] + rnorm(ncol(A[[x]])*g, 0, runif(1, 0.1, 0.2))
		A[[x]][,-(1:k[[x]][1])] <- A[[x]][,-(1:k[[x]][1])] + B[,x]
		rownames(A[[x]]) <- paste("ART", seq(1:g), sep="_")
		appendix <- rep("D",g)
		appendix[MU.negpos == 0] <- "U"	
		appendix[-id.sig] <- ""
		#rownames(A[[x]])[id.sig] <- paste(rownames(A[[x]])[id.sig], "S", sep="_")
		rownames(A[[x]])[id.sig] <- paste(rownames(A[[x]])[id.sig], appendix[id.sig], sep="_")
		samplenames <- paste("Sample",1:sum(k[2*x-1],k[2*x]), sep="")
		adf.data <- data.frame(group=c(rep("group1", k[2*x-1]), rep("group2", k[2*x])))
		rownames(adf.data) <- samplenames
		adf.meta <- data.frame(labelDescription=c("class sample x belongs to"))
		adf <- new("AnnotatedDataFrame")
		pData(adf) <- adf.data
		varMetadata(adf) <- adf.meta
		ES[[x]] <- new("ExpressionSet", exprs=A[[x]], phenoData=adf)
	}
	
	return(ES)
}
