\name{doEnsemble}
\Rdversion{1.1}
\alias{doEnsemble}
\title{
Method to perform an ensemble approach for meta-analysis.
}
\description{
This function allows performing a complete meta-analysis for a set of studies. Including an effect size, a significance and a ranking approach. Additionally an exploratory ensemble approach may be performed. The user may choose to plot some figures and report the complete results into a directory.
}
\usage{
doEnsemble(A, cl, cl.val = c(0, 1), feature.names = NULL, do.FM = TRUE, 
perm.ES = 1000, perm.RP = 1000, perm.ENS.RP = 1000, perm.ENS.RS = 1000, 
do.ENS.FM = TRUE, cluster = NULL, useREM = TRUE, results.dir = getwd(), 
write.all = FALSE, plot.fdr = FALSE, plot.ranks = FALSE)
}
\arguments{
  \item{A}{
A list of ExpressionSet objects, all sharing the same feature names.
}
  \item{cl}{
A list of factor variables, labeling the class for each of the samples in A.
}
  \item{cl.val}{
Level values of cl (default: 0,1)).
}
  \item{feature.names}{
names of features common to all sets. If NULL, featureNames(A[[1]]) is used.
}
  \item{do.FM}{
Doing meta analysis based on significance (Fisher method).
}
  \item{perm.ES}{
Number of permutations for calculating FDR for ES Choi.
}
  \item{perm.RP}{
Number of permutations for rank product approach.
}
  \item{perm.ENS.RP}{
Number of permutations for Ensemble rank product approach.
}
  \item{perm.ENS.RS}{
Number of permutations for Ensemble rank sum approach.
}
  \item{do.ENS.FM}{
Doing Ensemble Fisher-Method.
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
}
  \item{useREM}{
boolean parameter indicating whether to use ES REM model (default is TRUE)
}
  \item{results.dir}{
Directory where results are written as csv. 
The script will create a subfolde with date and time, where everything will be stored.
}
  \item{write.all}{
If TRUE all results are written as a csv in a subfolder of results.dir.
}
  \item{plot.fdr}{
Plot FDR figures for MA & Ensemble methods.
}
  \item{plot.ranks}{
Plot heatmaps for single MA methods.
}
}
\value{
The method returns the list of objects with all the results for both up- and downregulation. Wording for upregulation in the second group is g2up, and for downregulation in the second group is referred to as g2down. Additionally the target directory for all resulting output files is reported.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\note{
When performing an ensemble approach with meta-analysis data coming from the same set of studies, the result has to be considered to be biased towards the alternative hypothesis due to correlations within the data. The result has then to be considered for exploratory use only!
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <-  lapply(A, function(a){factor(as.numeric(a$group)-1)})
res <- doEnsemble(A, cl, feature.names=NULL, do.FM= TRUE, 
	   perm.ES=10, perm.RP=10, perm.ENS.RP= 10, perm.ENS.RS= 10, 
	   do.ENS.FM= TRUE, cluster= NULL, useREM= TRUE, 
	   results.dir=getwd(), write.all=FALSE, plot.fdr= FALSE, plot.ranks=FALSE)
names(res)
}
\keyword{ manip }

