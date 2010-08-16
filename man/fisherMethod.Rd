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
fisherMethod(data, zero.subst = default.zero.subst, cluster = NULL)
}
\arguments{
  \item{data}{
Data.frame with features in rows and studies in columns.
}
  \item{zero.subst}{
Since p-values of 0 cause problems, these are substituted by a very small number (default is \code{default.zero.subst}).
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
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
Karl Kugler, karl@eigenlab.net
}
\examples{

set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
pval <- multiTtest(A, cl, alternative= "two.sided")
rownames(pval) <- rownames(exprs(A[[1]]))
res <- fisherMethod(data=pval)
res
}
\keyword{univar}

