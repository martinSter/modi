% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modi-internal.R
\name{EA.dist}
\alias{EA.dist}
\title{Utility function for EAdet and EAimp}
\usage{
EA.dist(
  data,
  n,
  p,
  weights,
  reach,
  transmission.function,
  power,
  distance.type,
  maxl
)
}
\arguments{
\item{data}{a data frame or matrix with data.}

\item{n}{number of rows.}

\item{p}{number of columns.}

\item{weights}{a vector of positive sampling weights.}

\item{reach}{if \code{reach = "max"} the maximal nearest
neighbor distance is used as the basis for the transmission
function, otherwise the weighted}

\item{transmission.function}{form of the transmission function
of distance d: \code{"step"} is a heaviside function which jumps
to \code{1} at \code{d0}, \code{"linear"} is linear between
\code{0} and \code{d0}, \code{"power"} is \code{(beta*d+1)^(-p)}
for \code{p = ncol(data)} as default, \code{"root"} is the function
\code{1-(1-d/d0)^(1/maxl)}.}

\item{power}{sets \code{p = power}.}

\item{distance.type}{distance type in function \code{dist()}.}

\item{maxl}{maximum number of steps without infection.}
}
\description{
Calculation of distances for EPIDEMIC Algorithm for
multivariate outlier detection and imputation
}
\author{
Beat Hulliger
}
\keyword{internal}
