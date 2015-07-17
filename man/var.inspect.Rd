\name{var.inspect}
\alias{var.inspect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for calculating empirical variances with respect to time for data sets with regularly or irregularly spaced follow-up time points
}
\description{
Calculates empirical variances for data sets with regularly or irregularly spaced time points, and plots the result
}
\usage{
var.inspect(resid, timeVar, binwidth, numElems = 0, irregular = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resid}{
a vector of empirical residuals
}
  \item{timeVar}{
a vector for the time variable
}
  \item{binwidth}{
a numerical value for the bin length, to be used for irregularly spaced data
}
  \item{numElems}{
a numerical value for the elimination of the bins with less than that number of elements
}
  \item{irregular}{
a character string, \code{FALSE} indicates the data are collected at regular time points
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns mid values and variances of the bins, and numbers of elements falling into the bins for \code{irregular = TRUE}, 
and unique time points and variances, and number of the elements for the time points for \code{irregular = FALSE}. 
}
\references{
Asar O, Ritchie J, Kalra P, Diggle PJ (2015) Acute kidney injury amongst chronic kidney disease 
patients: a case-study in statistical modelling. To be submitted.

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

# obtaining empirical residuals by a linear model
# and calculating the empirical variances
lm.fit <- lm(log.egfr ~ sex + bage + fu + pwl, data = data.sim.ibm.short)
var.inspect(resid = resid(lm.fit), timeVar = data.sim.ibm.short$fu, binwidth = 0.1,
   numElems = 20, irregular = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variances}
