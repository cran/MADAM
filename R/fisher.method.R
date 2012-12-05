fisher.method <- function(pvals, method=c("fisher"), p.corr=c("bonferroni","BH","none"), zero.sub=0.00001, na.rm=FALSE, mc.cores=NULL){
  stopifnot(method %in% c("fisher"))
  stopifnot(p.corr %in% c("none","bonferroni","BH"))
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  if(is.null(dim(pvals)))
    stop("pvals must have a dim attribute")
  p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
  ##substitute p-values of 0
  pvals[pvals == 0] <- zero.sub
  if(is.null(mc.cores)){
    fisher.sums <- data.frame(do.call(rbind, apply(pvals, 1, fisher.sum, zero.sub=zero.sub, na.rm=na.rm)))
  } else {
    fisher.sums <- multicore::mclapply(1:nrow(pvals), function(i){
      fisher.sum(pvals[i,], zero.sub=zero.sub, na.rm=na.rm)
    }, mc.cores=mc.cores)
    fisher.sums <- data.frame(do.call(rbind, fisher.sums))
  }
    
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S, df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
                              bonferroni = p.adjust(fisher.sums$p.value, "bonferroni"),
                              BH = p.adjust(fisher.sums$p.value, "BH"),
                              none = fisher.sums$p.value)
  return(fisher.sums)
}
