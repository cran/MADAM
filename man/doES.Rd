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
doES(A, cl, cl.val = c(0, 1), useREM = TRUE, nperm = 5000, cluster = NULL)
}
\arguments{
  \item{A}{
List o ExpressionSets.
}
  \item{cl}{
List classes to be tested for first class. (should be in the first sample of the first study)
}
  \item{cl.val}{
Level values of cl (default: 0,1)).
}
  \item{useREM}{
Boolean indicating whether to use REM ES model.
}
  \item{nperm}{
Number of permutations to calculate FDR (default is 5000).
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
}
}
\value{
list with information on up- and down regulation of features
}   
\note{
This function is based on the implementation of the GeneMeta package Version 1.20.0 by Choi et al. Any modification made by MADAM developers are only for enabling cluster computing.
}
\references{
Choi et al, Combining multiple microarray studies and modeling
     interstudy variation. Bioinformatics, 2003, i84-i90.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\note{
Since the current implementation of MADAM only wraps to the method within the \pkg{GeneMeta} package, it cannot be run in a parallel manner.
}
\seealso{
\code{\link{zScoreFDR}}
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
res <- doES(A, cl, useREM = TRUE, nperm = 10, cluster = NULL)
res
}
\keyword{univar}
