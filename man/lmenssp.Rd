\name{lmenssp}
\alias{lmenssp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to obtain the maximum likelihood estimates of the parameters for linear mixed effects models with
random intercept and a stationary or non-stationary stochastic process component, under multivariate normal 
response distribution
}
\description{
Obtains the maximum likelihood estimates of the parameters for linear mixed effects models with
random intercept and a stationary or non-stationary stochastic process component, under multivariate 
normal response distribution
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
a character string for the stochastic process: \code{"bm"} for Brownian motion,
\code{"ibm"} for integrated Brownian motion, \code{"iou"} for integrated 
Ornstein-Uhlenbeck process, \code{"sgp-powered-power-method"} for stationary process 
with powered correlation function (for the usage, see the details below) and 
\code{"sgp-matern-kappa"} for stationary process with Matern correlation 
function (for the usage, see the details below)
}
  \item{timeVar}{
a numerical vector for the time variable
}
  \item{init}{
a numerical vector of initial values for the variance parameters to 
start the Fisher-Scoring or Nelder-Mead algorithms; 
if the user does not provide their own values \code{lmenssp} obtains these by fitting a 
random intercept and random slope model using the \code{lme} function of the \pkg{nlme} 
package; for details see the details below
}
  \item{tol}{
a numerical value for the maximum tolerance to assess the convergence
}
  \item{maxiter}{
a numerical value for the number of iterations for the Fisher-Scoring or Nelder-Mead algorithms
}
  \item{silent}{
a character string, if set to \code{FALSE} the details of the Fisher-Scoring steps are printed when 
the algorithm is running
}
}
\details{
For \code{"process"}, \code{"sgp-powered-power-method"} is a general form for stationary process 
with powered correlation function. \code{"power"} is the shapre parameter and 
corresponds to \code{"c"} in \deqn{exp(-|t-s|^c/\phi),}  
where \eqn{t} and \eqn{s} are two time points and \eqn{\phi} is the scale parameter,  
and \code{"method"} might be \code{"fs"} for Fisher-Scoring or \code{"nm"} for Nelder-Mead. 
Some examples are: \code{"sgp-powered-1-fs"} for stationary process with exponential correlation 
function with Fisher-Scoring algorithm and 
\code{"sgp-powered-2-nm"} for stationary process with Gaussian correlation 
function with Nelder-Mead algorithm. Similarly, \code{"sgp-matern-kappa"} is a general form for 
stationary process with Matern correlation 
function. \code{"kappa"} is the shape parameter and corresponds to \eqn{\kappa} in 
\deqn{\left\{ 2^{\kappa-1} \Gamma \left( \kappa \right) \right \}^{-1} \left(|t - s|/\nu \right)^{\kappa} K_{\kappa} \left(|t - s|/\nu \right),}
\eqn{t} and \eqn{s} are two time points and \eqn{\nu} is the scale parameter. 
An example is \code{"sgp-matern-0.5"} for 
stationary process with exponential correlation function. Nelder-Mead algorithm is automatically specified 
for the choice of Matern, i.e. Fisher-Scoring is not available. 

\code{"init"} assumes the following: 

- 3 element vectors for \code{"process = bm"} and \code{"process = ibm"}, with initials for the 
variances of random intercept, stochastic process and measurement error, respectively 
   
- 4 element vector for \code{"process = iou"}, with initials for the variance of 
random intercept variance, (two) parameters of the stochastic process, variance of measurement error

- 3 element vector for \code{"process = sgp-powered-power-fs"}, with initials for 
log of the fraction of the variance of random intercept and variance of the process, 
log(\eqn{\phi}) and log of the fraction of the variance of measurement error and 
variance of the process.

- \code{NULL} for any of the specification of \code{"process"}, in which case \code{lmenssp} 
finds the initals internally using the \code{lme} function of the \pkg{nlme} package.
}
\value{
Returns the results as lists
}
\references{
Diggle PJ (1988) An approach to the analysis of repeated measurements. \emph{Biometrics}, \bold{44}, 959-971.

Diggle PJ, Sousa I, Asar O (2015) Real time monitoring of progression towards renal failure in primary care patients.
\emph{Biostatistics}, 16(3), 522-536. 

Taylor JMG, Cumberland WG, Sy JP (1994) A Stochastic Model for Analysis of Longitudinal AIDS Data. 
\emph{Journal of the American Statistical Association}, \bold{89}, 727-736.

}
\author{
Ozgur Asar, Peter J. Diggle
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# loading the data set and subsetting it for the first 20 patients 
# for the sake illustration of the usage of the functions
data(data.sim.ibm)
data.sim.ibm.short <- data.sim.ibm[data.sim.ibm$id <= 20, ]

# fitting the model with integrated Brownian motion 
fit.ibm <- lmenssp(log.egfr ~ sex + bage + fu + pwl, data = data.sim.ibm.short,
  id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu, 
  silent = FALSE)
fit.ibm
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum likelihood estimation}
\keyword{stochastic processes}% __ONLY ONE__ keyword per line
