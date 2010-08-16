\name{generateRandomMAData}
\Rdversion{1.1}
\alias{generateRandomMAData}
\title{
Calculate the rank product and asses the significance
}
\description{
With this function the rank product for each feature (row) in data is calculated first. Then, the significance is calculated by permutating the ranks B times and comparing the derived random rank products with the original one. Since this might turn out to be computationally demanding, a snow cluster object has to be provided. 
}
\usage{
generateRandomMAData(g = 5000, perc.sig = 0.1, i = 2, 
k = 0, err.span = c(0.1, 0.2))
}
\arguments{
  \item{g}{
Number of features per study.
}
  \item{perc.sig}{
Percentage of overall significant features.
}
  \item{i}{
Number of studies.
}
  \item{k}{
Number of samples for each of the two groups/study. If set to 0, numbers are sampled randomly
}
  \item{err.span}{
Span for var in error.
}
}
\value{
A list of ExpressionSet Objets with features names "ART". Artificially down- or upregulated features are marked with "_D" or "_U"
}
\details{
In order to validate our ensemble approach, two data models were merged and adapted for our
purposes. This model driven approach enabled the testing and validating of our work on three
levels: First, on the single study level, second on the information derived from the meta-analysis, and third
we were able to inspect our results on the ensemble approach level. The basic assumption in this model is,
that in a set of G features, every feature has a baseline concentration, being the same for all studies. 
This baseline concentration may differ from study to study and from sample to sample.

If no study size information is provided they are sampled as: between 8 and 16 with probability 0.6, with probability 0.1 between 5 and 14 and with probability 0.3 between 15 and 30.


}
\references{
Hong F, Breitling R: A comparison of meta-analysis methods for detecting differentially expressed genes in microarray experiments. Bioinformatics 2008, 24(3):374-382

Kong X, Mas V, Archer KJ: A non-parametric meta-analysis approach for combining independent
microarray datasets: application using two microarray datasets pertaining to chronic allograft
nephropathy. BMC Genomics 2008, 9:98
}
\author{
Karl Kugler, karl@eigenlab.net
}
\examples{
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = 0)
A[[1]]
exprs(A[[1]])
}
\keyword{datagen}
