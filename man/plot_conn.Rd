% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estgtp.R
\name{plot_conn}
\alias{plot_conn}
\title{plot_conn}
\usage{
plot_conn(X, C)
}
\arguments{
\item{X}{The input set from the estgtp function}

\item{C}{Prepared parameter from the estgtp function}
}
\description{
Auxiliary function to plot partial results during evaluation of estgtp.
}
\details{
Auxiliary function which plots next step of partial results during
calculation of the estgtp function.
}
\examples{
library(spatstat)
kappa = 10
omega = .1
lambda= .5
theta = 10

X = rgtp(kappa, omega, lambda, theta, win = owin(c(0, 1), c(0, 1)))
plot_conn(X$X, X$C)

}
