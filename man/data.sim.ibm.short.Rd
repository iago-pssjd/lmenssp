\name{data.sim.ibm.short}
\alias{data.sim.ibm.short}
\docType{data}
\title{
Shorter version of \code{data.sim.ibm}
}
\description{
Includes data for the first 50 patients of \code{data.sim.ibm}.
}
\usage{data("data.sim.ibm.short")}
\format{
  A data frame with 820 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{a numeric vector: subject identification}
    \item{\code{sex}}{a numeric vector: 0 = male , 1 = female}
    \item{\code{bage}}{a numeric vector: baseline age in years}
    \item{\code{fu}}{a numeric vector: follow-up time in years}
    \item{\code{pwl}}{a numeric vector: piecewise linear term, 
                     calculated as, max(0, age-56.5)}
    \item{\code{log.egfr}}{a numeric vector: response variable representing 
                           natural logarithm of estimated glomerular filtration rate}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(data.sim.ibm.short)
str(data.sim.ibm.short)
}
\keyword{datasets}
