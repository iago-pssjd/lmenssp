\name{filtered}
\alias{filtered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for filtering
}
\description{
Filters random components of the mixed model with a non-stationary stochastic process component
}
\usage{
filtered(formula, data = NULL, id, process = "bm", timeVar, estimate, subj.id)
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
integrated Brownian motion, and \code{"iou"} for integrated 
Ornstein-Uhlenbeck process
}
  \item{timeVar}{
a vector for the time variable at which filtered values are required
}
  \item{estimate}{
a vector for the maximum likelihood estimates
}
  \item{subj.id}{
a vector for IDs of the subjects for whom filtering will be carried out
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
}
\author{
Ozgur Asar
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

# loading the package and the data set
library(lmenssp)
data(data.sim.ibm.short)

# model formula to be used below
formula <- log.egfr ~ sex + bage + fu + pwl

# obtaining the maximum likelihood estimates of the model
# parameters for the model with integrated Brownian motion
fit.ibm <- lmenssp(formula = formula, data = data.sim.ibm.short,
  id = data.sim.ibm.short$id, process = "ibm", timeVar = data.sim.ibm.short$fu, silent = FALSE)
fit.ibm

# filtering for subjects with ID=1 and 2
subj.id <- c(1, 2)
fil.res <- filtered(formula = formula, data = data.sim.ibm.short, 
    id = data.sim.ibm.short$id, process = "ibm",  timeVar = data.sim.ibm.short$fu, 
    estimate = fit.ibm$estimate[, 1], subj.id = subj.id)
fil.res

# filtering for a new (hypothetical) subject
data.501   <- data.frame(id = c(501, 501, 501), sex = c(0, 0, 0),
   bage = c(50, 50, 50), fu = c(0, 0.2, 0.4), 
   pwl = c(0, 0, 0), log.egfr = c(4.3, 2.1, 4.1))
data.501
fil.501 <- filtered(formula = formula, data = data.501, 
   id = data.501$id, process = "ibm", timeVar = data.501$fu,  
   estimate = fit.ibm$estimate[, 1], subj.id = 501)
fil.501

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{filtering}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
