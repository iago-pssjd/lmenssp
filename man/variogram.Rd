\name{variogram}
\alias{variogram}
\title{
A function for calculating the empirical variogram of a data set with irregularly spaced follow-up time points
}
\description{
Calculates empirical variogram for data sets with irregularly spaced time points, and 
plots the result
}
\usage{
variogram(resid, timeVar, id, binwidth, numElems = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resid}{
a vector of empirical residuals
}
  \item{timeVar}{
a vector for the time variable
}
  \item{id}{
a vector for subject identification
}
  \item{binwidth}{
a numerical value for the bin length
}
  \item{numElems}{
a numerical value for the elimination of the bins with less than that number 
of elements
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns mid values and means of the bins, and numbers of elements falling into the bins
}
%\references{
%% ~put references to the literature/web site here ~
%}
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

# obtaining empirical residuals by a linear model
# and calculating the empirical variogram
lm.fit <- lm(log.egfr ~ sex + bage + fu + pwl, data = data.sim.ibm.short)
variogram(resid = resid(lm.fit), timeVar = data.sim.ibm.short$fu, id = data.sim.ibm.short$id,
  binwidth = 0.1, numElems = 20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variogram}
