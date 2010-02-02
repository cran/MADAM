\name{calculateRankSum}
\Rdversion{1.1}
\alias{calculateRankSum}
\title{
Calculate the rank sum and assess the significance
}
\description{
With this function the rank sum for each feature (row) in data is calculated first. Then, the significance is calculated by permutating the ranks B times and comparing the derived random rank products with the original one (similar to the rank product). Since this might turn out to be computationally demanding, a snow cluster object has to be provided. 

The basic idea is similar to the one of the rank product.
}
\usage{
calculateRankSum(data, B = 10000, cluster = NULL)
}
\arguments{
  \item{data}{
A matrix containing rank information from several studies. Features are stored in the rows, and the single studies in the columns.
}
  \item{B}{
The number of permutations used to assess the significance of the rank sum.
}
  \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, a single core snow cluster has to be established nevertheless. Leaving it to NULL (default) will result in an error.
}
}
\value{
A data.frame object with 4 columns, where the rownames are the rownames from the input matrix data.
\item{RS}{The rank sum}
\item{p.value}{The significances of RS}
\item{q.value}{The p-values corrected for the FDR}
\item{rank}{The rank of a feature, derived by ranking the q.values. Ties are ranked randomly.}
}

\references{
RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis.
Hong F, Breitling R, McEntee CW, Wittner BS, Nemhauser JL, Chory J.
Bioinformatics. 2006 Nov 15;22(22):2825-7. Epub 2006 Sep 18.
}
\author{
Karl Kugler, UMIT
}
\seealso{
\code{\link{calculateRankProduct}}, \code{\link{doRP}}
}
\examples{
#examples using Rmpi were set inactive for automated testing
#for increased compability with  Windows
\dontrun{
library(snow)
library(Rmpi)
set.seed(123)
M <- sapply(1:5, sample, x=1:10, size=10)
clust <- makeCluster(1, type="MPI")
res <- calculateRankSum(data=M, B = 100, cluster = clust)
stopCluster(clust)
res
}
}
\keyword{ univar }
