\name{getNbSamples}
\Rdversion{1.1}
\alias{getNbSamples}
\title{
Generate random number of samples for studies.
}
\description{
Function to draw random sample sizes in studies if not provided.
}
\usage{
getNbSamples(size = 1)
}
\arguments{
  \item{size}{
Should be two times the number of simulated studies
}
}
\details{
Data generation for meta analysis validation.
}
\value{
A vector of length size, representing a paired number of samples per study.
}
\author{
Karl Kugler, UMIT
}
\examples{
set.seed(1234)
#implicit call of getNbSamples
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = 0)
#explicit call
k.expl <- getNbSamples(3*2)
k.expl
A <- generateRandomMAData(g = 100, perc.sig = 0.1, i = 3, k = k.expl)
}
\keyword{utilities}


