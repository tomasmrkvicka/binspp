#' @title Bayesian inference for Neyman-Scott point processes
#'
#' @description The Bayesian MCMC estimation of parameters for Thomas-type cluster
#'  point process with various inhomogeneities. It allows for inhomogeneity in
#'  (i) distribution of parent points, (ii) mean number of points in a cluster,
#' (iii) cluster spread. The package also allows for the Bayesian MCMC
#' algorithm for the homogeneous generalized Thomas process.
#' The cluster size is allowed to have a variance that is greater or less
#' than the expected value (cluster sizes are over or under dispersed).
#' Details are described in Dvořák, Remeš, Beránek & Mrkvička (2022)
#' (\doi{10.48550/arXiv.2205.07946}).
#'
#' @docType package
#' @name binspp
#' @author Tomas Mrkvicka <mrkvicka.toma@gmail.com> (author),
#' Jiri Dvorak <dvorak@karlin.mff.cuni.cz> (author),
#' Ladislav Beranek <beranek@jcu.cz> (author),
#' Radim Remes <inrem@jcu.cz> (author, creator)
#'
#' @useDynLib binspp, .registration=TRUE
#' @import spatstat
#' @import Rcpp
#' @importFrom spatstat.model ppm
#' @importFrom spatstat.geom area owin ppp im
#' @importFrom spatstat.random rpoint rpoispp
#' @importFrom mvtnorm rmvnorm
#' @importFrom cluster pam
#' @importFrom stats runif median pnorm quantile rnorm rpois vcov
#' @importFrom graphics abline axis barplot hist lines par points
## #' @importFrom FNN
## #' @importFrom Rcpp
#' @importFrom VGAM dgenpois0
## #' @importFrom nlme
## #' @importFrom spatstat
#' @note
#'     License: GPL-3
## #'     Encoding: UTF-8
#'
#' @references
#'     Anderson, C. Mrkvička T. (2020). Inference for cluster point
#'     processes with over- or under-dispersed cluster sizes,
#'     \emph{Statistics and computing} \strong{30}, 1573–1590,
#'     \doi{10.1007/s11222-020-09960-8}.
#'
#'     Kopecký J., Mrkvička T. (2016). On the Bayesian estimation
#'     for the stationary Neyman-Scott point processes,
#'     Applications of Mathematics \strong{61}/\strong{4}, 503-514.
#'     Available from: \url{https://am.math.cas.cz/am61-4/9.html}.
#'
#'     Dvořák, J., Remeš, R., Beránek, L., Mrkvička, T. (2022). binspp:
#'     An R Package for Bayesian Inference for Neyman-Scott Point Processes
#'     with Complex Inhomogeneity Structure. \emph{arXiv}.
#'     \doi{10.48550/ARXIV.2205.07946}.
NULL
