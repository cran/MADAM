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
doRP(A, cl, cl.val = c(0, 1), nperm = 1000, gene.names = NULL, cluster = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
List of ExpressionSets.
}
  \item{cl}{
List classes to be tested for first class. 
}
  \item{cl.val}{
Level values of cl (default: 0,1)). First sample of each study has to be of first type.
}
  \item{nperm}{
Number of permutations to calculate FDR (default is 1000)
}
  \item{gene.names}{
If not provided within the list of ExpressionSets, feature names can be set here.
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
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
Karl Kugler, karl@eigenlab.net
}
\note{
Since the current implementation of MADAM only wraps to the method within the \pkg{RankProd} package, it cannot be run in a parallel manner.
}

\seealso{
\code{\link{calculateRankProduct}}, \code{\link{calculateRankSum}}, \code{\link{RPadvance}}
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
res <- doRP(A, cl=cl, nperm = 10, cluster = NULL)
res
}
\keyword{univar}

