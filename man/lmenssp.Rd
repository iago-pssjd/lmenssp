\name{lmenssp}
\alias{lmenssp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to obtain the maximum likelihood estimates of the parameters for linear mixed effects models with
random intercept and a non-stationary stochastic process component
}
\description{
Obtains the maximum likelihood estimates of the parameters for linear mixed effects models with
random intercept and a non-stationary stochastic process component
}
\usage{
lmenssp(formula, data = NULL, id, process = "bm", timeVar, init = NULL, tol = 1e-05, 
maxiter = 100, silent = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a typical \code{R} formula for the fixed effects component of the model
}
  \item{data}{
a data frame from which the variables are to be extracted
}
  \item{id}{
a numerical vector for subject identification
}
  \item{process}{
a character string for the non-stationary process: \code{"bm"} for Brownian motion,
\code{"ibm"} for integrated Brownian motion and \code{"iou"} for integrated 
Ornstein-Uhlenbeck process
}
  \item{timeVar}{
a numerical vector for the time variable
}
  \item{init}{
a numerical vector of initial values to start the Fisher-Scoring algorithm; if the user does not 
provide their own values \code{lmenssp} obtains these by fitting a random intercept and random 
slope model using the \code{lme} function of the \pkg{nlme} package
}
  \item{tol}{
a numerical value for the maximum tolerance to assess the convergence
}
  \item{maxiter}{
a numerical value for the number of iterations for the Fisher-Scoring algorithm
}
  \item{silent}{
a character string, if set to \code{FALSE} the details of the Fisher-Scoring steps are printed when 
the algorithm is running
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns the results as lists
}
\references{
Diggle PJ, Sousa I, Asar O (2014) Real time monitoring of progression towards renal failure in primary care patients.
Tentatively accepted by \emph{Biostatistics}.

Taylor JMG, Cumberland WG, Sy JP (1994) A Stochastic Model for Analysis of Longitudinal AIDS Data. 
\emph{Journal of the American Statistical Association}, \bold{89}, 727-736.
}
\author{
Ozgur Asar
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# loading the package and the data set
library("lmenssp")
data("data.sim.ibm.short")

fit.ibm <- lmenssp(log.egfr ~ sex + bage + fu + pwl, data = data.sim.ibm.short,
  id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maxmimum likelihood estimation}
\keyword{stochastic processes}% __ONLY ONE__ keyword per line
