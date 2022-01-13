#' binspp Bayesian inference for Neyman-Scott point processes
#'
#' Simulation of generalized Thomas process, Bayesian MCMC estimation of parameters for generalized Thomas process.
#'
#' @docType package
#' @name binspp
#' @description This package introduces the Bayesian MCMC algorithm for generalized Thomas process,
#'    allowing cluster size to have a variance that is greater or less than the expected value
#'    (cluster sizes are over- or under-dispersed).
#' @author Tomas Mrkvicka <mrkvicka.toma@gmail.com> (author), Ladislav Beranek <beranek@jcu.cz> (maintainer), Radim Remes <inrem@jcu.cz> (maintainer)
#'
#' @useDynLib binspp
#' @import Rcpp VGAM spatstat FNN cluster mvtnorm
#' @note
#'     LinkingTo: Rcpp, RcppArmadillo, RcppEigen
#'
#'     License: GPL-3
#'
#'     Encoding: UTF-8
#'
NULL
