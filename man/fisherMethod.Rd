\name{fisherMethod}
\Rdversion{1.1}
\alias{fisherMethod}
\title{
Perform classical Fisher Method for combining p-values in a parallel manner.
}
\description{
Function for performing classical Fisher method by deriving significance from chi square distribution.
}
\usage{
fisherMethod(data, zero.subst = default.zero.subst)
}
\arguments{
  \item{data}{
Data.frame with features in rows and studies in columns.
}
  \item{zero.subst}{
Since p-values of 0 cause problems, these are substituted by a very small number (default is \code{default.zero.subst}).
}
}
\value{
list of significances and ranks
}
\references{
Fisher RA
Statistical Methods for Research Workers. 
Oliver & Boyd, Edinburgh 1925.
}
\author{
Karl Kugler, UMIT
}
\examples{
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
pval <- multiTtest(A, cl, val=cl.val, alternative= "two.sided")
rownames(pval) <- rownames(exprs(A[[1]]))
res <- fisherMethod(data=pval)
res
}
\keyword{univar}

