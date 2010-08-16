\name{corrNA}
\Rdversion{1.1}
\alias{corrNA}
\title{
Functions to impute missing values
}
\description{
Missing values might cause problems in the analysis pipeline. Within MADAM a simple imputation method was integrated. This function bases on calculating mean values and sampling variance robust values around it.
%
A more sophisticated method using k nearest neighbors as provided by the \pkg{impute} package, which is wrapped to by this function as well.
}
\usage{
corrNA(es, method = "madam", cl, cl.val = c(0, 1), conservative = TRUE, 
na.thres = 0.5, na.abs.thres = 2, exclude = TRUE)
}
\arguments{
  \item{es}{
A list of ExpressionSet objects, all sharing the same feature names.
}
  \item{method}{
Selection of the used imputation method. Current options are "madam" and "knn". If using "madam" the MADAM based method is used, by selecting "knn" the knn method from the impute package is used.
}
  \item{cl}{
A list of factor variables, labelling the class for each of the samples in es.
}
  \item{cl.val}{
Level values of cl (default: 0,1)).
}
  \item{conservative}{~
A conservative substitution at least a percentage of na.thres features needs to be present in each group to impute missing values with the 
}
  \item{na.thres}{
The percentage threshold of values that have at least to be present in order to impute missing values using the MADAM method.
}
  \item{na.abs.thres}{
An absolute number of entries needed to impute missing values, if mean should not be calculated the conservative (percentage) way.
}
  \item{exclude}{
If set to 'TRUE', features still containing missing values after imputation are removed from all studies.
}
}
\value{
The method returns the list of ExpressionSets that was provided, but now with missing values inserted and failed features excluded, if this option was set.
}
\references{
     Hastie, T., Tibshirani, R., Sherlock, G., Eisen, M., Brown, P. and
     Botstein, D., Imputing Missing Data for Gene Expression Arrays,
     Stanford University Statistics Department Technical report (1999),
     http://www-stat.stanford.edu/~hastie/Papers/missing.pdf 

     Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown,
     Trevor Hastie, Robert Tibshirani, David Botstein and Russ B.
     Altman, Missing value estimation methods for DNA microarrays
     BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525

}
\author{
Karl Kugler, karl@eigenlab.net
}
\seealso{
\code{\link{corrNA}}
}
\examples{
set.seed(666)
A <- generateRandomMAData(g = 20, perc.sig = 0.1, i = 3)
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
#generate some bad measurements
exprs(A[[1]])[1,1] <- NA
exprs(A[[1]])[4,2:4] <- NA
exprs(A[[2]])[3,4:6] <- NA
sapply(A, function(a){sum(is.na(exprs(a)))})
A <- corrNA(A, method="madam", cl)
sapply(A, function(a){sum(is.na(exprs(a)))})
}
\keyword{manip}
