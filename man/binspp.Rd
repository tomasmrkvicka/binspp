% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binspp-package.R
\name{binspp}
\alias{binspp}
\title{Bayesian inference for Neyman-Scott point processes}
\description{
The Bayesian MCMC estimation of parameters for Thomas-type cluster
 point process with various inhomogeneities. It allows for inhomogeneity in
 (i) distribution of parent points, (ii) mean number of points in a cluster,
(iii) cluster spread. The package also allows for the Bayesian MCMC
algorithm for the homogeneous generalized Thomas process.
The cluster size is allowed to have a variance that is greater or less
than the expected value (cluster sizes are over or under dispersed).
Details are described in Dvořák, Remeš, Beránek & Mrkvička (2022)
(\doi{10.48550/arXiv.2205.07946}).
}
\note{
License: GPL-3
}
\references{
Anderson, C. Mrkvička T. (2020). Inference for cluster point
    processes with over- or under-dispersed cluster sizes,
    \emph{Statistics and computing} \strong{30}, 1573–1590,
    \doi{10.1007/s11222-020-09960-8}.

    Kopecký J., Mrkvička T. (2016). On the Bayesian estimation
    for the stationary Neyman-Scott point processes,
    Applications of Mathematics \strong{61}/\strong{4}, 503-514.
    Available from: \doi{10.1007/s10492-016-0144-8}.

    Dvořák, J., Remeš, R., Beránek, L., Mrkvička, T. (2022). binspp:
    An R Package for Bayesian Inference for Neyman-Scott Point Processes
    with Complex Inhomogeneity Structure. \emph{arXiv}.
    \doi{10.48550/ARXIV.2205.07946}.
}
\author{
Tomas Mrkvicka <mrkvicka.toma@gmail.com> (author),
Jiri Dvorak <dvorak@karlin.mff.cuni.cz> (author),
Ladislav Beranek <beranek@jcu.cz> (author),
Radim Remes <inrem@jcu.cz> (author, creator),
Jaewoo Park <jwpark88@yonsei.ac.kr> (contributor),
Sujeong Lee <dltnwjd2304@gmail.com> (contributor)
}
