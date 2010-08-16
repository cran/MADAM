\name{plotMAVolcano}
\Rdversion{1.1}
\alias{plotMAVolcano}
\title{
Function for MA-volcano plot
}
\description{
Function to plot a volcano plot using information from Effect Size based MA and significance based MA.
}
\usage{
plotMAVolcano(x, title = "", col = "black", pch = ".", points = NULL,
points.col = "red", points.pch = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix with two columns (first column containing effect size, second the significance).
}
  \item{title}{
Title of plot.
}
  \item{col}{
Color for volcano plot.
}
  \item{pch}{
Character type for volcano plot.
}
  \item{points}{
A matrix with two columns: specifing if certain points are to be drawn separately.
}
  \item{points.col}{
Character type for plotting points.
}
  \item{points.pch}{
color for plotting points
}
}
\value{
Plots a MA-Volcano plot.
}
\author{
Karl Kugler, karl@eigenlab.net
}
\examples{
set.seed(123)
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = rep(5, 6))
cl <- lapply(A, function(a){factor(as.numeric(a$group)-1)})
res <- doES(A, cl, useREM = TRUE, nperm = 10, cluster = NULL)
plotMAVolcano(data.frame(res$g2up$MUvals, res$g2up$FDR))
}
\keyword{ hplot }
