\name{multiTtest}
\Rdversion{1.1}
\alias{multiTtest}
\title{
Function to perform t-tests on a list of ExpressionSets.
}
\description{
This function helps in executing a t-test (as implemented in the \pkg{stats} package) on a list of ExpressionSet objects, which is the main data type used in allmost all MADAM based functions.
}
\usage{
multiTtest(es, cl, val = c(0, 1), alternative = "two.sided")
}
\arguments{
  \item{es}{
List of ExpressionSets.
}
  \item{cl}{
List of factors to discriminate between the two groups.
}
  \item{val}{
Level values of cl (default: 0,1)).
}
  \item{alternative}{
Alternative hypothesis, either greater, less or two.sided (default)
  - greater: 	group 1 > group 2
  - less: 	group 1 < group 2
}
}
\value{
Data.frame with id's as rownames and p-values for each set as columns.
}
\author{
Karl Kugler, UMIT
}
\seealso{
stats::t.test
}
\examples{
set.seed(123)
library(genefilter)
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
\keyword{ univar }

