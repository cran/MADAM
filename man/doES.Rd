\name{doES}
\Rdversion{1.1}
\alias{doES}
\title{
Method to wrap parameters for Effect Size Meta Analysis.
}
\description{
Function hiding call of ES method by Choi implemented in GeneMeta package.
}
\usage{
doES(A, cl, val = c(0, 1), useREM = TRUE, nperm = 5000)
}
\arguments{
  \item{A}{
List o ExpressionSets.
}
  \item{cl}{
List classes to be tested for first class. (should be in the first sample of the first study)
}
  \item{val}{
Level values of cl (default: 0,1)).
}
  \item{useREM}{
Boolean indicating whether to use REM ES model.
}
  \item{nperm}{
Number of permutations to calculate FDR (default is 5000).
}
}
\value{
list with information on up- and down regulation of features
}
\references{
Choi et al, Combining multiple microarray studies and modeling
     interstudy variation. Bioinformatics, 2003, i84-i90.
}
\author{
Karl Kugler, UMIT
}
\note{
Since the current implementation of MADAM only wraps to the method within the \pkg{GeneMeta} package, it cannot be run in a parallel manner.
}
\seealso{
\code{\link{zScoreFDR}}
}
\examples{
library(GeneMeta)
set.seed(123)
es1 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
es2 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
es3 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
A <- list(es1, es2, es3)
cl=list(factor(c(rep(0,5), rep(1,5))),
factor(c(rep(0,5), rep(1,5))), factor(c(rep(0,5), rep(1,5))))
cl.val=c(0,1)
res <- doES(A, cl, val = cl.val, useREM = TRUE, nperm = 10)
res
}
\keyword{univar}
