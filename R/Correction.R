## class for performing basic validation of studies in use

# function to correct if variance within a feature is 0,
# to correct for this a minimal random noise is added to the existing values
# es: list of expressionsets to correct
# mean: mean of normally distributed noise (default 0)
# sd: standard deviation of normally distributed noise (default is mean of group-wise feature-wise sd-> 0, otherwise give our own)

# Value: list of corrected Expressionsets

#default.replace.sd <- 1*10^(-4)

corrVar <- function(es, cl, cl.val=c(0,1), mean= 0, sd=0, na.rm=TRUE){
  cat("starting variance correction...")
  res <- lapply(1:length(es), function(i){
    group1 <- exprs(es[[i]])[,cl[[i]] == cl.val[1]]
    group2 <- exprs(es[[i]])[,cl[[i]] == cl.val[2]]

    if(sd == 0){
      sd1 <- mean(rowSds(group1))
      sd2 <- mean(rowSds(group2))
    } else {
      sd1 <- sd
      sd2 <- sd
    }
    
    group1[apply(group1, 1, var, na.rm=na.rm) == 0,] <- group1[apply(group1, 1, var) == 0,] + rnorm(ncol(group1), mean, sd1)    
    group2[apply(group2, 1, var, na.rm=na.rm) == 0,] <- group2[apply(group2, 1, var) == 0,] + rnorm(ncol(group2), mean, sd2)

    exprs(es[[i]])[,match(colnames(group1), sampleNames(es[[i]]))] <- group1
    exprs(es[[i]])[,match(colnames(group2), sampleNames(es[[i]]))] <- group2
    return(es[[i]])
  })
  return(res)
}

# function to correct simply correct for missing values within a study
# missing values are replaced by group mean if available
# if there is only one value left, it will be resued and variance correction will be performed
# es: list of expressionsets to correct
# cl: list of factors defining class of samples
# cl.cal: values of cl (max 2)
#conservative implies that at max a perctenage of na.thres NA values may be in a group, otherwise up to n(group)-2 values can be NA per group
# na.thres: max percentage of NAs per group
# exclude: whether or not to exclude features that do no fulfil criteria from all sets
# Value: list of corrected Expressionsets
corrNA <- function(es, method="madam", cl, cl.val=c(0,1), conservative=TRUE, na.thres=0.5, na.abs.thres=2, exclude=TRUE){
  if(!method %in% c("madam","knn"))
    stop(paste(method, "is not an allowed for method"))
  names <- sapply(1:length(es), function(i){
    identical(featureNames(es[[1]]), featureNames(es[[i]]))
  })
  if(!all(names)){
    stop("feature names have to be identical for all sets!")
  }
  cat("starting NA correction...")
  filt <- c()
  for(i in 1:length(es)){
    if(method == "madam"){
      group1 <- exprs(es[[i]])[,cl[[i]] == cl.val[1]]
      group2 <- exprs(es[[i]])[,cl[[i]] == cl.val[2]]

      #get number of NAs
      nas1 <- apply(group1, 1, function(g)sum(is.na(g)))
      nas2 <- apply(group2, 1, function(g)sum(is.na(g)))
      sds1 <- rowSds(group1, na.rm=TRUE)
      sds2 <- rowSds(group2, na.rm=TRUE)

      #conservative implies that at max a perctenage of na.thres NA values may be in a group
      if(conservative){
	for(nr in 1:nrow(group1)){
	  if(ncol(group1) - nas1[nr] >= ncol(group1)*na.thres){
	    group1[nr,is.na(group1[nr,])] <- rnorm(nas1[nr], mean(group1[nr,], na.rm=TRUE), sds1[nr])
	  } else {
	    filt <- c(filt, nr)
	  }
	  if(ncol(group2) - nas2[nr] >= ncol(group2)*na.thres){
	    group2[nr,is.na(group2[nr,])] <- rnorm(nas2[nr], mean(group2[nr,], na.rm=TRUE), sds2[nr])
	  } else {
	    filt <- c(filt, nr)
	  }
	}
      } else {
	for(nr in 1:nrow(group1)){
		#if at least two NA and not all values NA, replace with Mean Group1
		if(nas1[nr] < ncol(group1)){
		  #too many NAS, not na.abs.thres not NAs left GROUP 1
		  if(ncol(group1)-na.abs.thres < nas1[nr]){
		    group1[nr,is.na(group1[nr,])] <- rnorm(nas1[nr], mean(group1[nr], na.rm=TRUE), mean(sds1, na.rm=TRUE))
		  } else {
		    group1[nr,is.na(group1[nr,])] <- rnorm(nas1[nr], mean(group1[nr,], na.rm=TRUE), sds1[nr])
		  }
		} else {
		  filt <- c(filt, nr)
		}
		#if at least two NA and not all values NA, replace with Mean Group2
		if(nas2[nr] < ncol(group2)){
		  #too many NAS, not 2 not NAs left GROUP 1
		  if(ncol(group2)-na.abs.thres < nas2[nr]){
		    group2[nr,is.na(group2[nr,])] <- rnorm(nas2[nr], mean(group2[nr], na.rm=TRUE), mean(sds2, na.rm=TRUE))
		  } else {
		    group2[nr,is.na(group2[nr,])] <- rnorm(nas2[nr], mean(group2[nr,], na.rm=TRUE), sds2[nr])
		  }
		} else {
		  filt <- c(filt, nr)
		}
	}
      }
      exprs(es[[i]])[,match(colnames(group1), sampleNames(es[[i]]))] <- group1
      exprs(es[[i]])[,match(colnames(group2), sampleNames(es[[i]]))] <- group2
    }
    if(method == "knn"){
      temp <- impute.knn(exprs(es[[i]]))$data
      for(nr in 1:nrow(temp)){
	if(sum(is.na(temp[nr,])) > 0){
	  filt <- c(filt, nr)
	}
      }
      exprs(es[[i]]) <- temp
    }

  }
  filt.out <- unique(filt)

  if(exclude && !is.null(filt.out)){
    es <- lapply(es, function(esx){
      return(esx[-filt.out,])
    })
  }
  return(es)
}

