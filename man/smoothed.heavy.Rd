\name{smoothed.heavy}
\alias{smoothed.heavy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for smoothing under multivariate t response distribution
}
\description{
Smooths random components of the mixed model with a stationary or non-stationary stochastic process 
component, under multivariate t response distribution
}
\usage{
smoothed.heavy(formula, data, id, process, timeVar, estimate, subj.id = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a typical R formula for the fixed effects component of the model
}
  \item{data}{
a data frame from which the variables are to be extracted
}
  \item{id}{
a vector for subject identification
}
  \item{process}{
a character string, \code{"bm"} for Brownian motion, \code{"ibm"} for 
integrated Brownian motion, \code{"iou"} for integrated 
Ornstein-Uhlenbeck process, \code{"sgp-powered-power-method"} for 
stationary process with powered correlation function, and 
\code{"sgp-matern-kappa"} for stationary process with Matern correlation function
}
  \item{timeVar}{
a vector for the time variable
}
  \item{estimate}{
a vector for the maximum likelihood estimates
}
  \item{subj.id}{
a vector of IDs of the subject for whom smoothing is to be carried out
}
}
\details{
For details of \code{"process"} see \code{lmenssp}
}
\value{
Returns the results as lists for the random intercept and 
stochastic process 
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
%%  ~~further notes~~
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

# a formula to be used 
formula <- log.egfr ~ sex + bage + fu + pwl

# estimating the parameters
# tol.em is set to 10^-1 and tol.lmenssp to 10^-2 only for illustration, 
# decrease these values in your applications
fit.heavy <- lmenssp.heavy(formula = formula, data = data.sim.ibm.heavy.short, 
       id = data.sim.ibm.heavy.short$id, timeVar = data.sim.ibm.heavy.short$fu, init.em = 5, 
       maxiter.em = 1000, tol.em = 10^-1, 
       process = "ibm", silent = FALSE, dof.est = c(0.1, 10, 0.0001), tol.cd = 0.001,
       tol.lmenssp = 10^-2, silent.lmenssp = FALSE)
fit.heavy

# smoothing for the patients with ID = 1, 2, 3, 4
smo.heavy <- smoothed.heavy(formula = formula, data = data.sim.ibm.heavy.short, 
      id = data.sim.ibm.heavy.short$id, process = "ibm", timeVar = data.sim.ibm.heavy.short$fu, 
      estimate = fit.heavy$est, subj.id = c(1, 2, 3, 4))
smo.heavy 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smoothing}
\keyword{multivariate t distribution}
