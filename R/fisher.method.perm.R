fisher.method.perm <- function(pvals, p.corr=c("bonferroni","BH","none"), zero.sub=0.00001, B=10000, mc.cores=NULL, blinker=1000){
  stopifnot(is.na(blinker)||blinker>0)
  stopifnot(p.corr %in% c("none","bonferroni","BH"))
  stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
  stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
  if(is.null(dim(pvals)))
    stop("pvals must have a dim attribute")
  p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
  pvals[pvals==0] <- zero.sub
  
  res.perm <- lapply(1:nrow(pvals),function(i){
    if(!is.na(blinker) & i%%blinker==0)
      message("=", appendLF = FALSE)
    ##which studies contribute to S (don't have a NA in row i)
    good.p <- which(!is.na(pvals[i,]))
    S.obs= fisher.sum(pvals[i,good.p], na.rm=FALSE)
    if(is.null(mc.cores)){
      S.rand <- unlist(lapply(1:B, function(b){
        ##get non NA p-values from studies contributing to S
        myp <- sapply(good.p, function(pc){
          sample(na.exclude(pvals[,pc]),1)
        })
        fisher.sum(myp)$S
      }))
    } else {
      S.rand <- unlist(mclapply(1:B, function(b){
        ##get non NA p-values from studies contributing to S
        myp <- sapply(good.p, function(pc){
          sample(na.exclude(pvals[,pc]),1)
        })
        fisher.sum(myp)$S
      }, mc.cores=mc.cores))
    }
    p.value <- sum(S.rand>=S.obs$S)/B
    data.frame(S=S.obs$S, num.p=S.obs$num.p, p.value=p.value)
  })
  res.perm <- data.frame(do.call(rbind, res.perm))
  
  if(!is.na(blinker) && blinker>0)
    message()
  ## rownames(res.perm) <- rownames(pvals)
  res.perm$p.adj <- switch(p.corr,
                              bonferroni = p.adjust(res.perm$p.value, "bonferroni"),
                              BH = p.adjust(res.perm$p.value, "BH"),
                              none = res.perm$p.value)
  return(res.perm)
  
}
