\name{plotHeatMap}
\Rdversion{1.1}
\alias{plotHeatMap}
\title{
Function for plotting ranklist heatmap.
}
\description{
This helps in plotting a heatmap of ranks from a set of studies. Typically the ranks are derived from the significances or a test statistic, but may of course represent effect sizes as well.
}
\usage{
plotHeatMap(rankList, title = "", top.col = NULL, side.col = NULL, 
col = rainbow(nrow(rankList)))
}
\arguments{
  \item{rankList}{
Data.frame containing ranks, rows representing the features and  columns the studies or methods...
}
  \item{title}{
Image title
}
  \item{top.col}{
Vector of colors used for plotting column headers (default= NULL)
}
  \item{side.col}{
Vector indicating how to print row side colors (default= NULL)
}
  \item{col}{
Colors to be used for plotting ranks (default: rainbow(nrow(rankList))
}
}
\value{
a plot
}
\details{
A ranklist heatmap is plotted using the function \code{heatmap.2} from the \code{gplots} package. Clustering is done using the Euclidian distance.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\seealso{
gplots::heatmap.2
}
\examples{
set.seed(123)
M <- sapply(1:5, sample, x=1:100, size=100)
plotHeatMap(M)
plotHeatMap(M, title = "", top.col = NULL, side.col = NULL, 
			col = rainbow(nrow(M)))
}
\keyword{ hplot }
