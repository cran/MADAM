\name{corrVar}
\Rdversion{1.1}
\alias{corrVar}
\title{
Function to replace variances of zero
}
\description{
In very rare cases measurements may show a variance of zero for some features. This might be induced by either technical or data manipulation means. Since such non-variance may cause problems in the later analysis, we integrated a method for dealing with this. First, we calculate the variance for all features within one group, and then use the mean m of it for adding some noise n to the features with zero variance: n ~ N(0,m) for each group separately.
}
\usage{
corrVar(es, cl, cl.val = c(0, 1), mean = 0, sd = 0, 
        na.rm = TRUE)
}
\arguments{
  \item{es}{
A list of ExpressionSet objects, all sharing the same feature names.
}
  \item{cl}{
A list of factor variables, labeling the class for each of the samples in es.
}
  \item{cl.val}{
Level values of cl (default: 0,1))..
}
  \item{mean}{
The mean used for N (default is set to zero)
}
  \item{sd}{
The standard deviation of n. It will be either calculated for each group, when setting it to zero (default), or can be given to the function.
}
  \item{na.rm}{
If set to 'TRUE' (default), missing values are excluded from calculating the variances.
}
}
\value{
The method returns the list of ExpressionSets that was provided, but now with variances of zero replaced by the new ones.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\seealso{
\code{\link{corrVar}}
}
\examples{
set.seed(666)
A <- generateRandomMAData(g = 20, perc.sig = 0.1, i = 3)
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
#generate some variance of 0
exprs(A[[1]])[1,] <- 5
min(rowSds(exprs(A[[1]])))
es <- corrVar(A, cl)
min(rowSds(exprs(es[[1]])))
}
\keyword{manip}

