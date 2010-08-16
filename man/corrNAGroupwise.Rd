\name{corrNAGroupwise}
\alias{corrNAGroupwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to filter features with missing values in a groupwise manner
}
\description{
Filtering for features with too many missing values is an important task. If feature information is not uniformly distributed between the experimental groups it may lead to problems with later analysis. An overall threshold for feature information may be reached for the whole study but not within each group.

This method checks for a certain ratio of numbers to present in every experimental group. By default bad features are filtered out (drop=TRUE), optionally the information about which features are bad can be returned.
}
\usage{
corrNAGroupwise(es = NULL, cl = NA, Th = 0.5, drop=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
A single ExpressionSet object or a list of ExpressionSet objects
}
  \item{cl}{
In case \code{es} is a single ExpressionSet object a factor with as many entries as \code{es} has samples, or in case \code{es} is a list, a list of the same length containing factors
}
  \item{Th}{
A percentage of at least \code{Th} measurements has to be present in each feature and group for a feature to be valid
}
  \item{drop}{
 By default bad features will be removed from the corresponding ExpressionSet, if \code{drop} is set to \code{FALSE} they will be kept and a list indicating which features are bad will be returned instead.
}

}
\value{
By default a list with ExpressionSet object will be returned. They are identical to the input \code{es}, but bad features are removed. If \code{drop} is set to \code{FALSE} a list will be returned providing bad feature information for each ExpressionSet object. This means that the indiced of those feature that do not fullfill the quality criterion are being listed.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\seealso{
\code{\link{corrVar}},\code{\link{corrNA}}
}
\examples{
set.seed(666)
es <- generateRandomMAData(g=10, i=3)

#generate some missingness
exprs(es[[1]])[1,1:5] <- NA

#prepare group information
es.g <- lapply(es, function(x){x$group})

#standard call with list
corrNAGroupwise(es=es, cl=es.g)[[1]]

#standard call with single object
corrNAGroupwise(es=es[[1]], cl=es.g[[1]])[[1]]

#do not drop bad features
corrNAGroupwise(es=es[[1]], cl=es.g[[1]], drop=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ NA }

