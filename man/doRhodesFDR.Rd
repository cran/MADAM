\name{doRhodesFDR}
\Rdversion{1.1}
\alias{doRhodesFDR}
\title{
Function for performing classical Rhodes MA method by deriving significance from p-values and permutating labels.
}
\description{
Perform method as proposed by Rhodes.
}
\usage{
doRhodesFDR(data, B = 10000, zero.subst = default.zero.subst, cluster = NULL)
}
\arguments{
  \item{data}{
Data frame or matrix containing p-values from various methods.
}
  \item{B}{
Number of permutations.
}
  \item{zero.subst}{
Value to replace zeros with.
}
  \item{cluster}{
Snow cluster object to perform permutations with.
}
}
\value{
A data.frame containing significances and ranks
}
\references{
Rhodes DR, Barrette TR, Rubin MA, Ghosh D, Chinnaiyan AM.
Meta-analysis of microarrays: interstudy validation of gene expression profiles reveals pathway dysregulation in prostate cancer.
Cancer Res. 2002 Aug 1;62(15):4427-33.
}
\author{
Karl Kugler, UMIT
}
\examples{
\dontrun{
set.seed(123)
es1 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
rownames(exprs(es1)) <- LETTERS[1:10]
es2 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
rownames(exprs(es2)) <- LETTERS[1:10]
es3 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
rownames(exprs(es3)) <- LETTERS[1:10]
A <- list(es1, es2, es3)
cl=list(factor(c(rep(0,5), rep(1,5))),
factor(c(rep(0,5), rep(1,5))), factor(c(rep(0,5), rep(1,5))))
cl.val=c(0,1)
clust <- makeCluster(1, type="MPI")
pval <- multiTtest(A, cl, val=cl.val, alternative= "two.sided")
rownames(pval) <- rownames(exprs(A[[1]]))
res <- doRhodesFDR(data=pval, B = 10, zero.subst = default.zero.subst, cluster = clust)
stopCluster(clust)
res
}
}
\keyword{univar}
