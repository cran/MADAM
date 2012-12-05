\name{plotPvsRank}
\alias{plotPvsRank}
\title{
Plot list of p-values against their ranks
}
\description{
This function plot a list of adjusted p-values against the index of the
entries ranked ascendingly by their p-value. This allows assesing how the distributions of p-values
compare between different studies or against a meta-analysis.
}
\usage{
plotPvsRank(pvals, s.legend = NULL, colos = NULL, ylim = NULL, log = "y", plot.title = NULL, zero.sub = 1e-05, legend.pos = "topright", lty = 1, lwd = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvals}{
A (named) list containing two data.frames or arrays with two
columns. The first columns contains the raw p-values and the second one
the adjusted ones.
}
  \item{s.legend}{
The names of the studies (used in the plot legend). Either this or names
of \code{pvals} have to be set.
}
  \item{colos}{
Colors to be used for the plotting. If less than 12 studies are given in
\code{pvals} and \code{colos} is not set, colors will be inferred from \code{RColorBrewer}.
}
  \item{ylim}{
Limits for plotting the y-axis (Must be between 0 and 1).
}
  \item{log}{
Which axis should be printed log scaled. Either ('y' (default),'x','xy' or '')
}
  \item{plot.title}{
    The plot title.
}
  \item{zero.sub}{
Replacemetn for 0s.
}
  \item{legend.pos}{
Position of the legend.
}
  \item{lty}{
Line with passed to \code{plot} and \code{legend}.
}
  \item{lwd}{
Line type passed to \code{plot} and \code{legend}.
}
}
\details{
  A graphical representation of studie p-values and meta-analysis
  results as given in Rhodes et al. 2002.

  If you have a named list and \code{s.legend} set, the later parameter
  will be used in the legend plot.
  
  Having \code{log} set to \code{y} or \code{xy} and the y-axis including 0, this
  will result in an error.

  In case you want to print the study raw p-values against their
  rank, just have the second column contain the values from the first one.
}
\value{
A plot
}
\references{
Rhodes, D. R., (2002). Meta-analysis of microarrays: interstudy validation of gene expression profiles reveals pathway dysregulation in prostate cancer. Cancer research, 62(15), 4427-33.
}
\author{
Karl Kugler <karl@eigenlab.net>
}

\examples{
pvals <- lapply(1:11, function(i){runif(1000)})

set.seed(12345)
pvals <- lapply(1:11, function(i){pnorm(rnorm(1000, sd=sample(c(1,1.2,1.4),1)), lower=TRUE)})
pvals <- lapply(pvals, function(x){cbind(x, p.adjust(x,method="BH"))})
names(pvals) <- paste("study",1:11)
plotPvsRank(pvals, lwd=2, lty=1)

##make line for last study look different
plotPvsRank(pvals, lwd=c(rep(1,10),2), lty=c(rep(1,10),2))

#log scale only x
plotPvsRank(pvals, log="x", lwd=2, lty=1)

#inverse y,  place legend at the bottom, and add title
plotPvsRank(pvals, log="y", ylim=c(0.0001,1),lwd=2, lty=1, legend.pos="bottomright", plot.title="my title")

}
