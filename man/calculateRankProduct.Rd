\name{calculateRankProduct}
\Rdversion{1.1}
\alias{calculateRankProduct}
\title{
Calculate the rank product and asses the significance
}
\description{
With this function the rank product for each feature (row) in data is calculated first. Then, the significance is calculated by permutating the ranks B times and comparing the derived random rank products with the original one. Since this might turn out to be computationally demanding, a snow cluster object has to be provided. 
}
\usage{
calculateRankProduct(data, B = 1000, cluster = NULL)
}
\arguments{
  \item{data}{
A matrix containing rank information from several studies. Features are stored in rows, and the single studies in columns.
}
  \item{B}{
The number of permutations used to assess the significance of the rank product.
}
 \item{cluster}{
A snow cluster object. If this package is used without parallel computing facilities, computing time may be much higher.
}
}
\value{
A data.frame object with 4 columns, where the rownames are the rownames from the input matrix data.
\item{RP}{The rank product}
\item{p.value}{The significances of RP}
\item{q.value}{The p-values corrected for the FDR}
\item{rank}{The rank of a feature, derived by ranking the q.values. Ties are ranked randomly.}
}
\note{
This function is based on the implementation of the RankProd package Version 2.20.0 by Hong F et al. Any modification made by MADAM developers are only for enabling cluster computing.
}
\references{
RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis.
Hong F, Breitling R, McEntee CW, Wittner BS, Nemhauser JL, Chory J.
Bioinformatics. 2006 Nov 15;22(22):2825-7. Epub 2006 Sep 18.

}
\author{
Karl Kugler, karl@eigenlab.net
}
\seealso{
\code{\link{calculateRankSum}}
}
\examples{
set.seed(123)
M <- sapply(1:5, sample, x=1:10, size=10)
res <- calculateRankProduct(data=M, B = 100, cluster = NULL)
res
}
\keyword{ univar }
