\name{lmenssp.heavy}
\alias{lmenssp.heavy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to obtain the maximum likelihood estimates of the parameters for linear mixed effects models with
random intercept and a stationary or non-stationary stochastic process component, under multivariate t 
response distribution
}
\description{
Obtains the maximum likelihood estimates of the parameters by expectatio-maximisation (E-M) algorithm 
for linear mixed effects models with
random intercept and a stationary or non-stationary stochastic process component, under multivariate 
normal t distribution
}
\usage{
lmenssp.heavy(formula, data, id, timeVar, init.em = NULL, maxiter.em = 100, 
tol.em = 0.001, process, silent = TRUE, dof.est = c(0.1, 10, 1e-04), 
tol.cd = 0.001, tol.lmenssp = 10^-5, init.lmenssp = NULL, maxiter.lmenssp = 100, 
silent.lmenssp = TRUE)
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
  \item{timeVar}{
a numerical vector for the time variable
}
  \item{init.em}{
a vector of inital values for the E-M algorithm
}
  \item{maxiter.em}{
a numerical value for the maximum number of iterations for the E-M algorithm
}
  \item{tol.em}{
a numerical value for the maximum tolerance to assess the convergence of the E-M algorithm
}
  \item{process}{
a character string for the stochastic process: \code{"bm"} for Brownian motion,
\code{"ibm"} for integrated Brownian motion, \code{"iou"} for integrated 
Ornstein-Uhlenbeck process, \code{"sgp-powered-power-method"} for stationary process 
with powered correlation function 
and \code{"sgp-matern-kappa"} for stationary 
process with Matern correlation function
}
  \item{silent}{
a character string, if set to \code{FALSE} the details of the E-M steps are printed when 
the algorithm is running
}
  \item{dof.est}{
a vector of three elements, to be passed to \code{optimize}
}
  \item{tol.cd}{
a numerical value for the tolerance of central-difference approximation
}
  \item{tol.lmenssp}{
a numerical value for the maximum tolerance to assess the convergence, to be passed to \code{lmenssp}
}
  \item{init.lmenssp}{
a vector of inital values, to be passed to \code{lmenssp}
}
  \item{maxiter.lmenssp}{
a numerical value of the number of iterations for the Fisher-Scoring or Nelder-Mead algorithms, to be passed to \code{lmenssp}
}
  \item{silent.lmenssp}{
a character string, if set to \code{FALSE} the details of the Fisher-Scoring steps are printed when 
the algorithm is running, to be passed to \code{lmenssp}
}
}
\details{
\code{lmenssp.heavy} calls \code{lmenssp} inside. 

Whilst theoretical standard errors are calculated and reported only for the fixed effects, 
central-difference based standard errors are calculated and reported for all the parameter estimates.

There are more than one way of specifying \code{init.em}, it can be set to: 

1) fixed effects, random effects parameters and degrees-of freedom, 

2) only the degrees-of-freedom, and 

3) NULL; for this specification, \code{lmenssp.heavy} finds the intials internally.

For the details of \code{"process"}, see \code{lmenssp}.

In \code{dof.est}, first and second elements are the minimum and maximum values of the search and the 
third element is the tolerance. It is passed to \code{optimize}.

}

\value{
Returns a list of results.
}
\references{
Asar O, Ritchie J, Kalra P, Diggle PJ (2015) Acute kidney injury amongst chronic kidney disease 
patients: a case-study in statistical modelling. To be submitted.

Pinheiro JC, Liu C, Wu YN. (2001) Efficient algorithms for robust estimation in linear mixed-effects 
models using the multivariate t distribution. \emph{Journal of Computational and Graphical Statistics} \bold{10}, 249-276.

}
\author{
Ozgur Asar, Peter J. Diggle
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# loading the data set and subsetting it for the first 20 patients 
# for the sake illustration of the usage of the functions
data(data.sim.ibm.heavy)
data.sim.ibm.heavy.short <- data.sim.ibm.heavy[data.sim.ibm.heavy$id <= 20, ]

# estimating the parameters
# tol.em is set to 10^-1 and tol.lmenssp to 10^-2 only for illustration, 
# decrease these values in your applications
fit.heavy <- lmenssp.heavy(formula = log.egfr ~ sex + bage + fu + pwl, 
       data = data.sim.ibm.heavy.short, id = data.sim.ibm.heavy.short$id, 
       timeVar = data.sim.ibm.heavy.short$fu, init.em = 5, maxiter.em = 1000, 
       tol.em = 10^-1, process = "ibm", silent = FALSE, 
       dof.est = c(0.1, 10, 0.0001), tol.cd = 0.001, tol.lmenssp = 10^-2, 
       silent.lmenssp = FALSE)
fit.heavy
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{maximum likelihood estimation}
\keyword{E-M algorithm}% __ONLY ONE__ keyword per line
