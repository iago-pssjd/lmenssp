\name{data.sim.ibm}
\alias{data.sim.ibm}
\docType{data}
\title{
A simulated data set under a mixed model with random intercept, integrated Brownian motion and 
multivariate normal response distribution 
}
\description{
A longitudinal data set containing measurements on 500 subjects, each with their own sequence of follow-up time points
collected at irregularly spaced follow-up time points.
}
\usage{data("data.sim.ibm")}
\format{
  A data frame with 8462 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{a numeric vector: subject identification}
    \item{\code{sex}}{a numeric vector: 0 = male, 1 = female}
    \item{\code{bage}}{a numeric vector: baseline age in years}
    \item{\code{fu}}{a numeric vector: follow-up time in years}
    \item{\code{pwl}}{a numeric vector: piecewise linear term, 
                     calculated as, max(0, age-56.5)}
    \item{\code{log.egfr}}{a numeric vector: response variable representing 
                           natural logarithm of estimated glomerular filtration rate}
  }
}
\details{
For the details of the true parameter settings, see Diggle, Sousa and Asar (2015). 
}
\references{
Diggle PJ, Sousa I, Asar O (2015) Real time monitoring of progression towards renal failure in primary care patients.
\emph{Biostatistics}, 16(3), 522-536.
}
\examples{
data(data.sim.ibm)
str(data.sim.ibm)
}
\keyword{datasets}
