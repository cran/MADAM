\name{fishersum}
\Rdversion{1.1}
\alias{fishersum}
\title{
Function to calculate Fisher's sum.
}
\description{
This function calculates the sum, that is used for Fisher's method of combining p-values.
}
\usage{
fishersum(P)
}
\arguments{
  \item{P}{
Vector of p-values.
}
}
\value{
Returns the sum of -2*log of the vector of p-values.
}
\references{
Fisher RA
Statistical Methods for Research Workers. 
Oliver & Boyd, Edinburgh 1925.
}
\author{
Karl Kugler, UMIT
}
\examples{
set.seed(123)
P1 <- c(1,1,1,1)
P2 <- c(0.05, 0.05, 0.05, 0.05)
S1 <- fishersum(P1)
S2 <- fishersum(P2)
1-pchisq(S1, df=4*2)
1-pchisq(S2, df=4*2)
}
\keyword{univar}

