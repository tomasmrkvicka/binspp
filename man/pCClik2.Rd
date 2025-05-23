% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pCClik2}
\alias{pCClik2}
\title{Evaluate unnormalized likelihood for auxiliary variable}
\usage{
pCClik2(thetaprop, CC)
}
\arguments{
\item{thetaprop}{A vector of theta1 and theta2.}

\item{CC}{A coordinate matrix of points.}
}
\value{
A list of the computed likelihood and a row vector of summed interaction between points.
}
\description{
Calculates the unnormalized likelihood for an auxiliary variable by evaluating pairwise interaction between points. The interaction thresholds are derived from the input theta vector.
}
