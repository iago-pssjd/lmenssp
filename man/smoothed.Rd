\name{smoothed}
\alias{smoothed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for smoothing under multivariate normal response distribution
}
\description{
Smooths random components of the mixed model with a stationary or non-stationary stochastic process 
component, under multivariate normal response distribution
}
\usage{
smoothed(formula, data = NULL, id, process = "bm", timeVar, estimate, 
subj.id = NULL, fine = NULL, eq.forec = NULL, uneq.forec = NULL)
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
  \item{fine}{
a numerical value for smoothing at fine intervals within the follow-up period 
}
  \item{subj.id}{
a vector of IDs of the subject for whom smoothing is to be carried out
}
  \item{eq.forec}{
a two element vector for equally spaced forecasting
}
  \item{uneq.forec}{
a two-column data frame or matrix for forecasting at desired time points
}
}
\details{
For details of \code{"process"} see \code{lmenssp}.
}
\value{
Returns the results as lists for the random intercept and 
stochastic process 
}
\references{
Asar O, Ritchie J, Kalra P, Diggle PJ (2015) Acute kidney injury amongst chronic kidney disease 
patients: a case-study in statistical modelling. To be submitted.

Diggle PJ (1988) An approach to the analysis of repeated measurements. \emph{Biometrics}, \bold{44}, 959-971.


Diggle PJ, Sousa I, Asar O (2015) Real time monitoring of progression towards renal failure in primary care patients.
\emph{Biostatistics}, 16(3), 522-536. 
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

# model formula to be used below
formula <- log.egfr ~ sex + bage + fu + pwl

# obtaining the maximum likelihood estimates of the model
# parameters for the model with integrated Brownian motion
fit.ibm <- lmenssp(formula = formula, data = data.sim.ibm.short,
  id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu, silent = FALSE)
fit.ibm

# smoothing for subject with ID=1 and 2
subj.id <- c(1, 2)
smo.res <- smoothed(formula = formula, data = data.sim.ibm.short, 
    id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu, 
    estimate = fit.ibm$estimate[, 1], subj.id = subj.id)
smo.res

# smoothing with fine interval of 0.01 within the follow-up period
smo.within <- smoothed(formula = formula, data = data.sim.ibm.short, 
    id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu,  
    estimate = fit.ibm$estimate[, 1], subj.id = subj.id, fine = 0.01)
smo.within

# one, two and three month forecasting for patients with IDs = 1 and 2
eq.forecast <- smoothed(formula = formula, data = data.sim.ibm.short, 
    id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu, 
    estimate = fit.ibm$estimate[, 1], subj.id = subj.id, 
    eq.forec = c(1/12, 3))
eq.forecast

# forecasting at arbitrary time points for patients with IDs = 1 and 2
uneq.forec <- data.frame(c(1, 1, 1, 2, 2), c(1/12, 2/12, 6/12, 1/12, 3/12))
uneq.forecast <- smoothed(formula = formula, data = data.sim.ibm.short, 
    id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu, 
    estimate = fit.ibm$estimate[, 1], uneq.forec = uneq.forec)
uneq.forecast

## smoothing for a new (hypothetical) patient 
data.501   <- data.frame(id = c(501, 501, 501), sex = c(0, 0, 0),
  bage = c(50, 50, 50), fu = c(0, 0.2, 0.4),
  pwl = c(0, 0, 0), log.egfr = c(4.3, 2.1, 4.1))
new.id     <-  501

# at observed time points
smo.501 <- smoothed(formula = formula, data = data.501, 
    id = data.501$id, process = "ibm", timeVar = data.501$fu, 
    estimate = fit.ibm$estimate[, 1], subj.id = new.id)
smo.501

# at fine interval of 0.01 within the follow-up period
smo.within.501 <- smoothed(formula = formula, data = data.501, 
    id = data.501$id, process = "ibm", timeVar = data.501$fu,  
    estimate = fit.ibm$estimate[, 1], subj.id = new.id, fine = 0.01)
smo.within.501

# one, two and three month forecasting 
eq.forecast.501 <- smoothed(formula = formula, data = data.501, 
    id = data.501$id, process = "ibm", timeVar = data.501$fu, 
    estimate = fit.ibm$estimate[, 1], subj.id = new.id, 
    eq.forec = c(1/12, 3))
eq.forecast.501

# forecasting at arbitrary time points 
uneq.forec.501 <- data.frame(c(501, 501, 501), c(1/12, 2/12, 4/12))
uneq.forecast.501 <- smoothed(formula = formula, data = data.501, 
    id = data.501$id, process = "ibm", timeVar = data.501$fu, 
    estimate = fit.ibm$estimate[, 1], uneq.forec = uneq.forec.501)
uneq.forecast.501
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smoothing}
\keyword{multivariate normal distribution}% __ONLY ONE__ keyword per line
