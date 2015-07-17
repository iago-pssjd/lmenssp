\name{qqplot.t}
\alias{qqplot.t}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile-quantile plot for univariate t distribution
}
\description{
Calculates empirical quantiles of univariate data and theoretical quantiles 
of a t distribution with a given degrees-of-freedom
}
\usage{
qqplot.t(x, dof, print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector or column matrix for data  
}
  \item{dof}{
a numerical value for degrees-of-freedom
}
  \item{print}{
a logical variable, if \code{"print = TRUE"} empirical and theoretical quantiles 
are reported and the user can use them to draw a new plot
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a list of output is returned if \code{"print = TRUE"}
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

\seealso{
\code{\link{qqplot}}
}
\examples{
set.seed(1)
rs <- rt(500, 5)
qqplot.t(rs, 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{quantile-quantile plot}
\keyword{univariate t distribution}% __ONLY ONE__ keyword per line
