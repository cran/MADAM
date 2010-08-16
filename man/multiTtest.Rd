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
multiTtest(es, cl, cl.val = c(0, 1), alternative = "two.sided", cluster = NULL)
}
\arguments{
  \item{es}{
List of ExpressionSets.
}
  \item{cl}{
List of factors to discriminate between the two groups.
}
  \item{cl.val}{
Level values of cl (default: 0,1)).
}
  \item{alternative}{
Alternative hypothesis, either greater, less or two.sided (default)
  - greater: 	group 1 > group 2
  - less: 	group 1 < group 2
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
}
}
\value{
Data.frame with id's as rownames and p-values for each set as columns.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\seealso{
stats::t.test
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
pval <- multiTtest(A, cl, alternative= "two.sided", cluster = NULL)
rownames(pval) <- rownames(exprs(A[[1]]))
res <- fisherMethod(data=pval)
res
}
\keyword{ univar }

