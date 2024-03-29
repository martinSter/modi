% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotIT.R
\name{plotIT}
\alias{plotIT}
\title{Plot of  infection times of the EA algorithm}
\usage{
plotIT(infection.time, weights, cutpoint)
}
\arguments{
\item{infection.time}{vector of infection.times of the observations}

\item{weights}{vector of (survey) weights of the observations}

\item{cutpoint}{a cutpoint to for declaring outliers}
}
\description{
The (weighted) cdf of infection times is plotted. The infection times
jumps of the cdf are shown by the points with the
same infection times stacked vertically and respecting the weights.
}
\details{
The infection times of \code{EAdet} are the main input. In addition the weights
may be needed. The default cutpoint from \code{EAdet} may be used for the cutpoint.
Points that are never infected have a missing infection time. These missing infection times
are (temporarily) imputed by 1.2 times the maximum infection time
to show them on the plot marked with an x.
}
\examples{
it <- c(rep(NA, 3), rep(1:7, times=c(1, 4, 10, 8, 5, 3, 2)))
wt <- rep(c(1,2,5), times=12)
plotIT(it, wt, 6)
}
\author{
Beat Hulliger
}
