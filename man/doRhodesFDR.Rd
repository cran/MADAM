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
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
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
Karl Kugler, karl@eigenlab.net
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
pval <- multiTtest(A, cl, alternative= "two.sided")
rownames(pval) <- rownames(exprs(A[[1]]))
res <- doRhodesFDR(data=pval, B = 10, zero.subst = default.zero.subst, cluster = NULL)
res
}
\keyword{univar}
