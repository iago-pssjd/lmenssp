\name{boot.nm}
\alias{boot.nm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate bootstrap standard errors
}
\description{
Calculates bootstrap standard errors for the parameter estimates obtained 
by \code{lmenssp} when Nelder-Mead algorithm is used
}
\usage{
boot.nm(formula, data, id, timeVar, result, matern = TRUE, kappa.or.power, 
nboot = 100, tol.lmenssp = 1e-08, maxiter.lmenssp = 500)
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
  \item{result}{
a matrix of results obtained by \code{lmenssp}, see the example below
}
  \item{matern}{
a logical variable, \code{TRUE} corresponds to Matern correlation function, \code{FALSE} corresponds to powered 
correlation function 
}
  \item{kappa.or.power}{
a numerical value for the shape parameter, it corresponds to \eqn{\kappa} if \code{matern = TRUE} 
and \eqn{\phi} if \code{matern = FALSE}
}
  \item{nboot}{
a numerical value for number of bootstrap sample
}
  \item{tol.lmenssp}{
a numerical value for the tolerance, to be passed to \code{lmenssp}
}
  \item{maxiter.lmenssp}{
a numerical value for the maximum number of iterations, to be passed to \code{lmenssp}
}
}
\details{
This function consider parametric bootstrap based on the fitted model. 
The recommended number of bootstrap replications is at least 100. 
For the details of \eqn{\kappa} and \eqn{\phi} in \code{kappa.or.power}, 
see the details section of \code{lmenssp} function.
}
\value{
Returns a list of results
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
# loading the data set and subsetting it for the first 5 patients 
# for the sake illustration of the usage of the functions
data(data.sim.ibm)
data.sim.ibm.short <- data.sim.ibm[data.sim.ibm$id <= 5, ]

# model formula to be used below
formula <- log.egfr ~ sex + bage + fu + pwl

# fitting the mixed model with Matern, kappa = 0.5
fit.matern <- lmenssp(formula = formula, data = data.sim.ibm.short,
  id = data.sim.ibm.short$id, process = "sgp-matern-0.5", timeVar = data.sim.ibm.short$fu, 
  init = c(-13, 1, -1), silent = FALSE)
fit.matern

# bootstrapping the standard errors, nboot is set to 2 for illustration
# set nboot to at least 100 in your applications
fit.matern.boot <- boot.nm(formula = formula, data = data.sim.ibm.short, 
                           id = data.sim.ibm.short$id, timeVar = data.sim.ibm.short$fu, 
                           result = fit.matern$est, matern = TRUE, kappa.or.power = 0.5, 
                           nboot = 2)
fit.matern.boot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bootstrap}
\keyword{standard error calculation}% __ONLY ONE__ keyword per line
