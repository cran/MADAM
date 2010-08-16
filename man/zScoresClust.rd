\name{zScoresClust}
\Rdversion{1.1}
\alias{zScoresClust}
\title{
Tool for Meta-analysis of gene expression data.
}
\description{
A meta-analysis function for computing zScores with optional cluster.
}
\usage{
zScoresClust(esets, classes, useREM=TRUE, 
             CombineExp=1:length(esets), cluster = NULL)
}
\arguments{
  \item{esets}{
A 'list' of 'ExpressionSet's, one expression set per
          experiment.  All experiments must have the same
          variables(genes).
}
  \item{classes}{
A 'list' of class memberships, one per experiment. Each
          'list' can only contain 2 levels.}
  \item{useREM}{
A 'logical' value indicating whether or not to use a REM,
          'TRUE', or a FEM, 'FALSE', for combining the z scores.
}
 \item{CombineExp}{
'vector' of integer- which experiments should be
          combined-default:all experiments
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
}
}
\details{
     The function 'zScores' implements the approach of Choi et al. for
     for a set of 'ExpressionSet's. The function 'zScoreFDR' computes a FDR for each gene, both for each
     single experiment and for the combined experiment. The FDR is
     calculated as described in Choi et al. Up to now ties in the
     zscores are not taken into account in the calculation. The
     function might produce incorrect results in that case. The
     function also computes zScores, both for the combines experiment
     and for each single experiment.
     
     'Clust'-functions are the same functions, extended for making cluster-computation possible and optional
}
\value{
 A 'matrix' with one row for each probe(set) and the following
     columns:

zSco_Ex_: For each single experiment the standardized mean difference,
          'Effect_Ex_', divided by the estimated standard deviation,
          the square root of the 'EffectVar_Ex_' column.

  MUvals: The combined standardized mean difference (using a FEM or
          REM)

   MUsds: The standard deviation of the 'MUvals'.

    zSco: The z statistic - the 'MUvals' divided by their standard
          deviations, 'MUsds'.

   Qvals: Cochran's Q statistic for each gene.

      df: The degree of freedom for the Chi-square distribution.  This
          is equal to the number of combined experiments minus one.

Qpvalues: The probability that a Chi-square random variable, with 'df'
          degrees of freedom) has a higher value than the value from
          the Q statistic.

   Chisq: The probability that a Chi-square random variate (with 1
          degree of freedom) has a higher value than the value of
          zSco^2.

Effect_Ex_: The standardized mean difference for each single
          experiment.

EffectVar_Ex_: The variance of the standardized mean difference for
          each single experiment.
     Note that the three column names that end in an underscore are
     replicated, once for each experiment that is being analyzed.

}
\references{
Choi et al, Combining multiple microarray studies and modeling
     interstudy variation. Bioinformatics, 2003, i84-i90.
}
\author{
 M. Ruschhaupt
 Clustering option added by Matthias Wieser, UMIT
}
\seealso{
\code{\link{zScores}}
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
res <- zScoresClust(A, cl, useREM = TRUE, 
                    CombineExp=1:length(A), cluster = NULL)
res
}
\keyword{univar}
