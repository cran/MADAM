\name{doRP}
\Rdversion{1.1}
\alias{doRP}
\title{
Function for calling RankProduct meta-analysis method by Breitling.
}
\description{
Class to wrap parameters for Rank Product Meta-Analysis by Breilting.
}
\usage{
doRP(A, cl, val = c(0, 1), nperm = 1000, gene.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
List of ExpressionSets.
}
  \item{cl}{
List classes to be tested for first class. 
}
  \item{val}{
Level values of cl (default: 0,1)). First sample of each study has to be of first type.
}
  \item{nperm}{
Number of permutations to calculate FDR (default is 1000)
}
  \item{gene.names}{
If not provided within the list of ExpressionSets, feature names can be set here.
}
}
\value{
List of significances and ranks for both up- and downregulation.
}
\references{
 Breitling, R., Armengaud, P., Amtmann, A., and Herzyk,  P.
     Rank Products:A simple, yet powerful, new method to  detect
     differentially regulated genes in replicated microarray
     experiments, _FEBS Letter_, 57383-92, 2004

}
\author{
Karl Kugler, UMIT
}
\note{
Since the current implementation of MADAM only wraps to the method within the \pkg{RankProd} package, it cannot be run in a parallel manner.
}

\seealso{
\code{\link{calculateRankProduct}}, \code{\link{calculateRankSum}}, \code{\link{RPadvance}}
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
res <- doRP(A, cl=cl, val = cl.val, nperm = 10)
res
}
\keyword{univar}

