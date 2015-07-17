\name{data.sim.ibm.heavy}
\alias{data.sim.ibm.heavy}
\docType{data}
\title{
A simulated data set under a mixed model with random intercept, integrated Brownian motion 
and multivariate t response distribution
}
\description{
A longitudinal data set containing measurements on 100 subjects, each with their own sequence of follow-up time points
collected at irregularly spaced follow-up time points.
}
\usage{data("data.sim.ibm.heavy")}
\format{
  A data frame with 1636 observations on the following 6 variables.
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
For the true parameter settings of the fixed effects and variance parameters, see Diggle, Sousa and Asar (2015). 
True degrees-of-freedom parameter is set to 5.
}
%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Asar O, Ritchie J, Kalra P, Diggle PJ (2015) Acute kidney injury amongst chronic kidney disease 
patients: a case-study in statistical modelling. To be submitted.

Diggle PJ, Sousa I, Asar O (2015) Real time monitoring of progression towards renal failure in primary care patients.
\emph{Biostatistics}, 16(3), 522-536. 

}
\examples{
data(data.sim.ibm.heavy)
str(data.sim.ibm.heavy)
}
\keyword{datasets}
