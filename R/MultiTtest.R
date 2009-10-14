## class for performing two group t-test on a list of ExpressionSets

##TODO: test wether rownames are in the same order for all sets
## best way by applying lapply -> building intersection 
## and if length of intersection == number of features(A[1])
## everyting is ok
## use match to reorder rownames

## function to perform t-tests on a list of ExpressionSets
## A: list of ExpressionSets
## cl: list of factors to discriminate between the two groups
## val: level values of cl (default: 0,1))
## side: either greater, less or two.sided (default)
## - greater: group 1 > group 2
## - less: group 1 < group 2
##
##Value: data.frame with ids as rownames and p-values for each set as columns
multiTtest <- function(es, cl, val=c(0,1), alternative= "two.sided"){
  if(!alternative %in% c("two.sided", "less", "greater"))
    stop(paste(alternative, "is not an allowed for alternative"))

  if(typeof(es) != "list"){
    stop("es has to be a list of ExpressionSets")
  }
  ##TODO: check length of es
  ##TODO: check type of es only ExpressionSet
  cat("performing multiple tests...\n")
  if(!all(levels(factor(val)) == levels(factor(unlist(cl))))){
    stop("factor levels not the same as val")
  }
  if(length(levels(unlist(cl))) != 2){
    stop("multiTTest works with two factors only")
  }
  
  #transform es to list of exprs-matrices
  A <- lapply(es, exprs)

  #test each matrix
  res <- sapply(1:length(A), function(i){
    print(i)
    f <- cl[[i]]
    #split groups, by using values group defintion are set to be correct
    #otherwise the order might differ!
    G1 <- A[[i]][,f==val[[1]]]
    G2 <- A[[i]][,f==val[[2]]]
    
    #test each gene
    temp <- sapply(1:nrow(G1), function(j){
      ltry <- try(t.test(x=G1[j,], y=G2[j,], alternative=alternative)$p.value, silent=TRUE)
      if(inherits(ltry,"try-error")) {
	return(NA)
      } else {
	return(ltry)
      }
    })
    cbind(temp)
  })

  #set names
  rownames(res) <- rownames(A[[1]])
  colnames(res) <- paste("p", 1:ncol(res), sep="")
  #exclude NAs
  res <- na.exclude(res)

  return(res)
}
