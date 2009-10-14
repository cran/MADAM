\name{plotFDR}
\Rdversion{1.1}
\alias{plotFDR}
\title{
Plot results significance against the rank of the significance.
}
\description{
For plotting the significances versus the ranks, the significances are first ranked. If ties occur, the first observered significance will get the lowest rank.

This helps in assesing the distribution of significances for meta analysis. Three types of studies can be annotated: single studies, meta-analysis of studies, and ensemble methods that combine meta-analysis.
}
\usage{
plotFDR(x, title = "", types = factor(rep("s", ncol(x))))
}
\arguments{
  \item{x}{
A matrix with the rownames being features and colnames being studies/methods. Entries are significances.
}
  \item{title}{
Title for image.
}
  \item{types}{
Three level factor indicating if a x comes from a normal method or a ensemble approach, difference will be in color and line type allowed levels are "s" single study, "m" meta analysis, and "e" ensemble approach.
}
}
\value{
a plot
}
\author{
Karl Kugler, UMIT
}
\examples{
set.seed(123)
es1 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
rownames(exprs(es1)) <- LETTERS[1:10]
es2 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
rownames(exprs(es2)) <- LETTERS[1:10]
es3 <- new("ExpressionSet", exprs=matrix(runif(100,1,10), nrow=10))
rownames(exprs(es3)) <- LETTERS[1:10]
A <- list(es1, es2, es3)
cl=list(factor(c(rep(0,5), rep(1,5))),
factor(c(rep(0,5), rep(1,5))), factor(c(rep(0,5), rep(1,5))))
cl.val=c(0,1)
pval <- multiTtest(A, cl, val=cl.val, alternative= "two.sided")
rownames(pval) <- rownames(exprs(A[[1]]))
res <- fisherMethod(data=pval)
PP <- data.frame(pval, res$p.value)
colnames(PP) <- c(paste("study", 1:3), "fisher")
plotFDR(PP, types=c(rep("s", 3), "m"))
}
\keyword{ hplot }

