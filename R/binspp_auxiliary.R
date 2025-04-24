utils::globalVariables(c("W"))

#' Simulate a realization of Thomas-type cluster point process with complex inhomogeneities
#'
#' @description The means to simulate realizations
#'       from the Thomas-type cluster point process with complex
#'       inhomogeneities are provided.
#'
#' @details A realization of a Thomas-type cluster
#'          point process model with possible inhomogeneity
#'          (described by covariates) are produced in any or all of the following model
#'          components: intensity function of the parent process, mean number
#'          of points in a cluster, scale of the clusters.
#'          Model parametrization is described in the documentation to the
#'          function [estintp()]. The parent process is generated in the dilated
#'          observation window \emph{W_dil} to avoid edge-effects,
#'          the resulting point pattern is eventually truncated to the smaller
#'          observation window \emph{W}.
#'
#' @param kappa intensity or intensity function of the parent process, scalar or pixel image object of class [spatstat.geom::im()] from the \pkg{spatstat} package.
#' @param alpha scalar, influences the mean number of points in individual clusters, see Details.
#' @param omega scalar, influences the spread of individual clusters, see Details.
#' @param W the observation window where the realization is to be generated, in the \cr [spatstat.geom::owin()] format of the \pkg{spatstat} package.
#' @param W_dil the observation window dilated by the assumed maximal cluster radius, as a binary mask with the same resolution as the covariates.
#' @param betavec vector of parameters describing the dependence of the intensity function of the parent process on covariates in the list \emph{z_beta}.
#' @param alphavec vector of parameters describing the dependence of the mean number of points in a cluster on covariates in the list \emph{z_alpha}.
#' @param omegavec vector of parameters describing the dependence of the spread of the clusters on covariates in the list \emph{z_omega}.
#' @param z_beta list of covariates describing the intensity function of the parent process, each covariate being a pixel image as used in the \pkg{spatstat} package.
#' @param z_alpha list of covariates describing the location-dependent mean number of points in a cluster, each covariate being a pixel image as used in the \pkg{spatstat} package.
#' @param z_omega list of covariates describing the location-dependent scale of a cluster, each covariate being a pixel image as used in the \pkg{spatstat} package.
#'
#' @return A planar point pattern, object of the type [spatstat.geom::ppp()]
#'         used in the \pkg{spatstat} package.
#'
#' @md
#' @examples
#'
#' library(spatstat)
#' # Unit square observation window:
#' W <- owin()
#'
#' # Dilation of the observation window:
#' W_dil <- dilation(W, 0.1)
#' W_dil <- as.mask(W_dil)
#'
#' # Define covariates:
#' f1 <- function(x, y) { x }
#' f2 <- function(x, y) { y }
#' f3 <- function(x, y) { 1 - (y - 0.5) ^ 2 }
#' cov1 <- as.im(f1, W = W_dil)
#' cov2 <- as.im(f2, W = W_dil)
#' cov3 <- as.im(f3, W = W_dil)
#'
#'
#' # Stationary Thomas process:
#' X <- rThomasInhom(kappa = 50, alpha = log(10), omega = log(0.01),
#'        W = W, W_dil = W_dil)
#' plot(X)
#'
#'
#' # Thomas-type cluster process with inhomogeneity in all model components:
#' X <- rThomasInhom(kappa = 10, betavec = c(1), z_beta = list(cov1),
#'             alpha = log(10), alphavec = c(1), z_alpha = list(cov2),
#'             omega = log(0.01), omegavec = c(1), z_omega = list(cov3),
#'             W = W, W_dil = W_dil)
#' plot(X)
#'
#' @export
rThomasInhom <- function(kappa, alpha, omega, W, W_dil,
                         betavec = NULL, alphavec = NULL, omegavec = NULL,
                         z_beta = NULL, z_alpha = NULL, z_omega = NULL){

  if (length(z_beta)  != length(betavec)){ stop("Error: length of z_beta does not equal the length of betavec!")}
  if (length(z_alpha) != length(alphavec)){ stop("Error: length of z_alpha does not equal the length of alphavec!")}
  if (length(z_omega) != length(omegavec)){ stop("Error: length of z_omega does not equal the length of omegavec!")}

  x.out <- NULL
  y.out <- NULL

  if (!is.null(betavec)){
    B = betavec[1] * z_beta[[1]]
    if (length(betavec) > 1){
      for (i in 2:length(betavec)){B = B + betavec[i] * z_beta[[i]]}
    }
    intensity_parents <- kappa * exp(B)
    parents <- rpoispp(intensity_parents)
  } else {
    parents <- rpoispp(kappa, win = W_dil)
  }

  if (parents$n > 0){
    if (!is.null(alphavec)){
      B = alphavec[1] * z_alpha[[1]]
      if (length(alphavec) > 1){
        for (i in 2:length(alphavec)){B = B + alphavec[i] * z_alpha[[i]]}
      }
      alpha_surface <- exp(alpha + B)
      alphas <- alpha_surface[parents]
    } else {
      alphas <- rep(exp(alpha), times=parents$n)
    }
    if (!is.null(omegavec)){
      B = omegavec[1] * z_omega[[1]]
      if (length(omegavec)>1){
        for (i in 2:length(omegavec)){B = B + omegavec[i] * z_omega[[i]]}
      }
      omega_surface <- exp(omega + B)
      omegas <- omega_surface[parents]
    } else {
      omegas <- rep(exp(omega), times = parents$n)
    }

    n.offsprings <- rpois(n = parents$n, lambda = alphas)

    for (i in 1:parents$n){
      if (n.offsprings[i] == 0){ next }
      aux.dis <- rnorm(n = 2 * n.offsprings[i], mean = 0, sd = omegas[i])

      x.aux <- parents$x[i] + aux.dis[1:n.offsprings[i]]
      y.aux <- parents$y[i] + aux.dis[(n.offsprings[i] + 1):(2 * n.offsprings[i])]
      x.out <- c(x.out, x.aux)
      y.out <- c(y.out, y.aux)
    }

  }
  else return(rpoispp(0))

  ok <- spatstat.geom::inside.owin(x.out, y.out, W)
  out <- ppp(x = x.out[ok], y = y.out[ok], window = W)
  return(out)
}


# @TODO
intrho = function(kappa, bet, z_beta){
  s=0*z_beta[[1]]
  for (i in 1:length(bet)){
    s = s + bet[i] * z_beta[[i]]
  }
  s2 = 1 * exp(s)
  return(sum(s2$v, na.rm = TRUE) * z_beta[[1]]$xstep * z_beta[[1]]$ystep)
}

#' @title Estimate the first-order inhomogeneity
#'
#' @description For exploratory purposes it may be useful to perform the first
#'       step of the analysis only, to investigate the dependence of
#'       the intensity function of the parent process on given covariates,
#'       without running the MCMC chain.
#'
#' @details The calling the [spatstat.model::ppm()] function from the \pkg{spatstat}
#'          package, with some additional computations useful when preparing
#'          the run of the MCMC chain, is mainly performed in this function.
#'          The function also contains a simple way to plot the estimated
#'          intensity function of the parent process.
#'
#' @param X observed point pattern in the [spatstat.geom::ppp()] format of the \pkg{spatstat} package.
#' @param z_beta list of covariates describing the intensity function of the parent process, each covariate being a pixel image as used in the \pkg{spatstat} package.
#' @param W_dil the observation window dilated by the assumed maximal cluster radius.
#' @param plot logical, should the estimates intensity function of the parent process be plotted?
#'
#' @return List containing the output of the [spatstat.model::ppm()]
#'         function from the \pkg{spatstat} package, along with some auxiliary
#'         objects useful for running the MCMC chain.
#'
#' @md
#' @examples
#'
#' library(spatstat)
#' # Prepare the dataset:
#' X = trees_N4
#' x_left = x_left_N4
#' x_right = x_right_N4
#' y_bottom = y_bottom_N4
#' y_top = y_top_N4
#'
#' z_beta = list(refor = cov_refor, slope = cov_slope)
#'
#' # Determine the union of rectangles:
#' W = owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))
#' if (length(x_left) >= 2){
#'   for (i in 2:length(x_left)){
#'     W2 = owin(c(x_left[i], x_right[i]), c(y_bottom[i], y_top[i]))
#'     W = union.owin(W, W2)
#'   }
#' }
#'
#' # Dilated observation window:
#' W_dil = dilation.owin(W, 100)
#'
#'
#' # Estimating the intensity function of the parent process:
#' aux = first_step(X, z_beta, W_dil, plot = TRUE)
#'
#' @export
first_step <- function(X, z_beta, W_dil, plot = TRUE){
  if (length(z_beta) > 0){
    # Estimate the log-linear model for intensity function depending on covariates
    FS = ppm(X ~ ., covariates = z_beta)
    bet=FS$coef[-1]

    # Estimated intensity function
    B = bet[1] * z_beta[[1]]
    if (length(bet) > 1){
      for (i in 2:length(bet)){B = B + bet[i] * z_beta[[i]]}
    }
    B = exp(B) # Covariate part of intensity
    if (plot){plot(exp(FS$coef[1]) * B,
                   main = "Plot of estimated 1st-order intensity function")}

    # Estimated integral of the intensity function over W_dil
    integralrho = intrho(1, bet, z_beta)
  } else {
    # No covariates provided, estimating homogeneous model
    FS = spatstat.model::ppm(X ~ 1)

    # Estimated intensity function
    B = spatstat.geom::as.im(1, W = W_dil)
    if (plot) {plot(exp(FS$coef) * B,
                    main = "Plot of estimated 1st-order intensity function")}

    # Estimated integral of the intensity function over W_dil
    integralrho = area(W_dil)
  }

  return (list(FS = FS, B = B, integralrho = integralrho))
}


NewCenterPoint = function(W_dil, B){
  NK = TRUE
  while(NK){
    P = rpoint(n=1,f=B,win=W_dil)
    if(spatstat.geom::inside.owin(P$x[1],P$y[1],W_dil)) NK=FALSE
  }
  c(P$x[1],P$y[1])
}


pbetXC = function(Y, z_beta, z_alpha, z_omega, NStep, SamplingFreq, x_left,
                  x_right,y_bottom, y_top, output.first_step, control,
                  AreaW, AreaMRW, W_dil, Wpix, verbose = TRUE){
  # Set parameters of prior distributions, if not specified in the "control" list
  if (is.null(control$Prior_alpha_mean)){Prior_alpha_mean = 3.0} else {Prior_alpha_mean = control$Prior_alpha_mean}
  if (is.null(control$Prior_alpha_SD)){Prior_alpha_SD = 2.0} else {Prior_alpha_SD = control$Prior_alpha_SD}
  if (is.null(control$Prior_omega_mean)){Prior_omega_mean = log(sqrt(AreaW)/20)} else {Prior_omega_mean = control$Prior_omega_mean}
  if (is.null(control$Prior_omega_SD)){Prior_omega_SD = log(3+sqrt(AreaW)/40)} else {Prior_omega_SD = control$Prior_omega_SD}

  if (is.null(control$Prior_alphavec_SD)){
    if (length(z_alpha)>0){
      Prior_alphavec_SD = rep(0, times=length(z_alpha))
      for (i in 1:length(z_alpha)){
        Prior_alphavec_SD[i] = 2/max(z_alpha[[i]]$v, na.rm=TRUE)
      }
    } else {Prior_alphavec_SD = NULL}
  } else {Prior_alphavec_SD = control$Prior_alphavec_SD}

  if (is.null(control$Prior_omegavec_SD)){
    if (length(z_omega)>0){
      Prior_omegavec_SD = rep(0, times=length(z_omega))
      for (i in 1:length(z_omega)){
        Prior_omegavec_SD[i] = 2/max(z_omega[[i]]$v, na.rm=TRUE)*log(3+sqrt(AreaW)/20)
      }
    } else {Prior_omegavec_SD = NULL}
  } else {Prior_omegavec_SD = control$Prior_omegavec_SD}

  aux.podil <- 30

  salpha = Prior_alpha_SD/aux.podil;       # standard deviation for update of alpha
  somega = Prior_omega_SD/aux.podil;       # standard deviation for update of omega
  salphabet = Prior_alphavec_SD/aux.podil # standard deviation for update of alphabet
  somegabet = Prior_omegavec_SD/aux.podil # standard deviation for update of omegabet

  if (length(z_alpha) != length(Prior_alphavec_SD)){stop("Error: length of z_alpha does not equal the length of Prior_alphavec_SD!")}
  if (length(z_omega) != length(Prior_omegavec_SD)){stop("Error: length of z_omega does not equal the length of Prior_omegavec_SD!")}

  while((alpha <- rnorm(1,Prior_alpha_mean,Prior_alpha_SD)) <Prior_alpha_mean || alpha > (Prior_alpha_mean+Prior_alpha_SD)){};
  while((omega <- rnorm(1,Prior_omega_mean,Prior_omega_SD)) <Prior_omega_mean || omega > (Prior_omega_mean+Prior_omega_SD)){};
  if(length(salphabet)>0){alphabet <- rmvnorm(1,mean = rep(0,times=length(salphabet)), sigma = salphabet*diag(length(salphabet))/100)}else {alphabet <- NULL};
  if(length(somegabet)>0){omegabet <- rmvnorm(1,mean = rep(0,times=length(somegabet)), sigma = somegabet*diag(length(somegabet))/100)}else {omegabet <- NULL};
  lambda = exp(output.first_step$FS$coef[1]);
  kappa = max(lambda/exp(alpha),0.00001);
  CC = rpoispp(kappa*output.first_step$B);
  CC = t(rbind(CC$x,CC$y))

  integralrho = output.first_step$integralrho

  integral = KumulaVsechC(CC,z_alpha,z_omega,alpha,alphabet,omega,omegabet,x_left,x_right,y_bottom,y_top)
  logP = logpXCbetC(Y,CC, z_alpha, z_omega, alpha,alphabet,omega,omegabet,AreaW, integral)

  if (length(z_beta)>0){
    pBX_0 = c(kappa, alpha, alphabet, omega, omegabet, logP, integral, rep(1,length(output.first_step$FS$coef[-1])), 0, 0, 0, 0, 0, 0);
  } else {
    pBX_0 = c(kappa, alpha, alphabet, omega, omegabet, logP, integral, NULL, 0, 0, 0, 0, 0, 0);
  }
  pBX = matrix(NA, nrow=1+floor(NStep/SamplingFreq), ncol=length(pBX_0))
  pBX[1,] = pBX_0

  alphaAccepts = rep(0, times=1000)
  omegaAccepts = rep(0, times=1000)
  parentAccepts = rep(0, times=1000)

  start.time = Sys.time()
  if (verbose) {
    cat("Start of run: ")
    cat(as.character(as.POSIXct(start.time, origin="1970-01-01"),usetz=F),sep="")
    cat("\n")
  }

  for(step in 1:NStep){

    if (verbose) {
      if ((step %% SamplingFreq) == 0){cat(".")}
      if ((step %% 1000) == 0){
        cat(" Iteration no. ")
        cat(step)
        cat("; estimated end of run: ")
        cat(as.character(as.POSIXct((NStep-step)*(Sys.time()-start.time)/step + Sys.time(), origin="1970-01-01"),usetz=F),sep="")
        cat("\n")
      }
    }
    S = StepbetC(kappa, z_alpha, z_omega, alpha, salpha, alphabet, salphabet, omega, somega, omegabet, somegabet,
                 Y, CC, logP, integral, integralrho, x_left, x_right, y_bottom, y_top, Wpix, AreaW, output.first_step$FS$coef[1],
                 Prior_alpha_mean, Prior_alpha_SD, Prior_omega_mean, Prior_omega_SD, Prior_alphavec_SD, Prior_omegavec_SD)
    kappa = S$kappa;
    alpha = S$alpha;
    alphabet = S$alphabet;
    omega = S$omega;
    omegabet = S$omegabet;
    logP = S$logP;
    integral = S$integral;
    alphaAccepts[1+(step %% 1000)] = S$alphaAccept;
    omegaAccepts[1+(step %% 1000)] = S$omegaAccept;

    NewCenter = NewCenterPoint(W_dil,output.first_step$B)
    U = StepMovePointC(kappa, z_alpha, z_omega, alpha, alphabet, omega, omegabet, Y, CC, logP, integral, integralrho, x_left, x_right, y_bottom, y_top, AreaW, NewCenter)
    CC = U$CC;
    logP = U$logP;
    integral = U$integral;
    parentAccepts[1+(step %% 1000)] = U$parentAccept;

    if( (step %% SamplingFreq) == 0){
      if (length(z_beta)>0){
        XX = suppressWarnings(ppp(x=CC[,1],y=CC[,2],window=W_dil))
        SS = suppressWarnings(spatstat.model::ppm(XX ~ ., covariates = z_beta))
        pv = rep(0,length(output.first_step$FS$coef[-1]))
        for(i in 2:length(SS$coef)){
          pv[i-1]=2*(1-pnorm(abs(SS$coef[i])/sqrt(vcov(SS)[i,i]),0,1))
        }
      } else {
        pv = NULL
      }
      pBX[1+floor(step/SamplingFreq),] = c(S$kappa,S$alpha,S$alphabet,S$omega,S$omegabet,U$logP,S$integral,pv,dim(CC)[1],min(1,exp(S$logProbAcceptAlpha)), min(1,exp(S$logProbAcceptOmega)), mean(alphaAccepts), mean(omegaAccepts), mean(parentAccepts))

      # plot(Y)
      # points(CC, col="red", pch=20)
    }
  }
  return(list(pBX=pBX, priorParameters=list(Prior_alpha_mean=Prior_alpha_mean,Prior_alpha_SD=Prior_alpha_SD,Prior_omega_mean=Prior_omega_mean,Prior_omega_SD=Prior_omega_SD,Prior_alphavec_SD=Prior_alphavec_SD,Prior_omegavec_SD=Prior_omegavec_SD)))
}


#' @title Estimation of Thomas-type cluster point process with complex inhomogeneities
#'
#' @description The Bayesian MCMC estimation of
#'     parameters for Thomas-type cluster point process with inhomogeneity
#'     is performed in any of the following parts: (i) distribution of parent
#'     points, (ii) mean number of points in a cluster, (iii) cluster spread.
#'     The process is observed in the observation window \emph{W} which
#'     is a union of aligned rectangles, aligned with the coordinate axes.
#'     The inhomogeneities are described through a parametric model
#'     depending on covariates. The estimation algorithm is described in
#'     Dvořák, Remeš, Beránek & Mrkvička (2022)
#'     (\doi{10.48550/arXiv.2205.07946}).
#'
#' @details
#' # Details
#' ## Parametric model
#' The model for the intensity function of the parent process is the following:
#' \eqn{f(u) = kappa * exp(beta_1 * z\_beta_1(u) + … + beta_k * z\_beta_k(u))},
#' where \eqn{(kappa, beta_1, …, beta_k)} is the vector of parameters
#' and \eqn{z\_beta = (z\_beta_1, …, z\_beta_k)} is the list of covariates.
#' Note that choosing \eqn{k = 0} is acceptable, resulting in a homogeneous
#' distribution of parents. In such a case \emph{z_beta} must be an empty list
#' or NULL. Furthermore, the list z_beta must contain named covariates in order
#' to properly function with the function [spatstat.model::ppm()] from the
#' \pkg{spatstat} package
#' which is used in the first step to estimate the parameters
#' \eqn{(kappa, beta_1, …, beta_k)}. Note that due to identifiability issues
#' the covariate lists \emph{z_beta} and \emph{z_alpha} must be disjoint.
#'
#' The model for the mean number of points in a cluster corresponding
#' to the parent at location \eqn{u} is the following:
#' \eqn{g(u) = exp(alpha + alpha_1 * z\_alpha_1(u) + … + alpha_l * z\_alpha_l(u))},
#' where \eqn{(alpha, alpha_1, …, alpha_l)} is the vector of parameters and
#' \eqn{z\_alpha = (z\_alpha_1, …, z\_alpha_l)} is the list of covariates.
#' Note that choosing \eqn{l = 0} is acceptable, resulting in a constant model.
#' In such a case \emph{z_alpha} must be an empty list or NULL. Note that
#' due to identifiability issues the covariate lists \emph{z_beta} and
#' \emph{z_alpha} must be disjoint.
#'
#' The model for the scale of a cluster corresponding to the parent at
#' location \eqn{u} is the following:
#' \eqn{h(u) = exp(omega + omega_1 * z\_omega_1(u) + … + omega_m * z\_omega_m(u))},
#' where \cr
#' \eqn{(omega, omega_1, …, omega_m)} is the vector of parameters and \cr
#' \eqn{z\_omega = (z\_omega_1, …, z\_omega_m)} is the list of covariates.
#' Note that choosing \eqn{m = 0} is acceptable, resulting in a constant model.
#' In such a case \emph{z_omega} must be an empty list or NULL.
#'
#' ## Observation window and its dilation
#' The observation window must be provided as the union of aligned rectangles,
#' aligned with the coordinate axes. This, however, allows the analysis
#' of point patterns observed in rather irregular regions by approximating
#' the region by a union of aligned rectangles. The structure of the vectors
#' \emph{x_left, x_right, y_bottom} and \emph{y_top} is such that the first
#' rectangle is constructed using the function [spatstat.geom::owin()] from the
#' \pkg{spatstat} package as
#' \code{owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))},
#' and similarly for the other rectangles. Naturally, a rectangular window
#' can be used and in such a case the vectors \emph{x_left} to \emph{y_top}
#' each contain a single element.
#'
#' ## Covariates
#' The covariates must be provided as pixel images of the class
#' [spatstat.geom::im()] used in the \pkg{spatstat} package. It is recommended
#' that all the covariates have the same pixel resolution. However, it is
#' necessary that all the covariates in the list \emph{z_beta} have the same
#' resolution, all the covariates in the list \emph{z_alpha} have the same
#' resolution and all the covariates in the list \emph{z_omega} have
#' the same distribution. The covariates must be provided in the dilated
#' observation window \emph{W_dil}, with NA values at pixels lying outside
#' \emph{W_dil}.
#'
#' ## Control
#' The control list must contain the following elements: \emph{NStep}
#' (the required number of MCMC iterations to be performed), \emph{BurnIn}
#' (burn-in, how many iterations at the beginning of the chain will be
#' disregarded when computing the resulting estimates – note that this choice
#' can be updated after the computation without re-running the chain,
#' see the function [re_estimate()]), \emph{SamplingFreq} (sampling frequency
#' for estimating the posterior distributions). Additionally,
#' the hyperparameters for the prior distributions should be given, see below.
#' Note that some default values for the hyperparameters are provided but it is
#' \strong{strongly encouraged} that the hyperparameter values are given
#' by the user, based on the actual knowledge of the problem at hand.
#'
#' ## Prior distributions and hyperparameters
#' The prior distribution for \emph{alpha} is normal with
#' \code{mean = Prior_alpha_mean} and \cr \code{SD = Prior_alpha_SD}.
#'
#' The prior distribution for the vector \eqn{(alpha_1, …, alpha_l)} is
#' centered normal with diagonal variance matrix and the vector of
#' \code{SDs = Prior_alphavec_SD}.
#'
#' The prior distribution for \emph{omega} is normal with
#' \code{mean = Prior_omega_mean} and \cr \code{SD = Prior_omega_SD}.
#'
#' The prior distribution for the vector \eqn{(omega_1, …, omega_m)}
#' is centered normal with diagonal variance matrix and the vector of
#' \code{SDs = Prior_omegavec_SD}.
#'
#' The hyperparameters should be provided in the control list. However,
#' the following default choices are applied if the hyperparameter values
#' are not provided by user or are given as NULL: \cr
#' \code{Prior_alpha_mean = 3},
#' \code{Prior_alpha_SD = 2}, \code{Prior_omega_mean = log(sqrt(area(W) / 20))},
#' \code{Prior_omega_SD = log(3 + sqrt(area(W) / 40))},
#' \code{Prior_alphavec_SD[i] = 2 / max(z_alpha_i)},
#' \code{Prior_omegavec_SD[i] = 2 / max(z_omega_i) * log(3 + sqrt(area(W) / 20))}.
#'
#' ## Output
#' The output of the function is given by the list containing the parameter
#' estimates along with the 2.5% and 97.5% quantiles of the posterior
#' distributions. Also, several auxiliary objects are included in the list
#' which are needed for the [print_outputs()] and [plot_outputs()] functions.
#'
#'
#' @param X observed point pattern in the [spatstat.geom::ppp()] format of
#'          the \pkg{spatstat} package.
#' @param control list specifying various tuning constants for the MCMC
#'                estimation. See also Details.
#' @param x_left vector describing the observation window, contains
#'               the lower x-coordinate of the corners of each rectangle.
#' @param x_right vector describing the observation window, contains
#'                the higher x-coordinate of the corners of each rectangle.
#' @param y_bottom vector describing the observation window, contains the smaller y-coordinate of the corners of each rectangle.
#' @param y_top vector describing the observation window, contains the higher y-coordinate of the corners of each rectangle.
#' @param W_dil the observation window dilated by the assumed maximal cluster radius.
#' @param z_beta list of covariates describing the intensity function of the parent process, each covariate being a pixel image as used in the \pkg{spatstat} package.
#' @param z_alpha list of covariates describing the location-dependent mean number of points in a cluster, each covariate being a pixel image as used in the \pkg{spatstat} package.
#' @param z_omega list of covariates describing the location-dependent scale of a cluster, each covariate being a pixel image as used in the \pkg{spatstat} package.
#' @param verbose logical (TRUE or FALSE). For suppressing information messages to the console set value to FALSE. Defaults to TRUE.
#'
#' @return The output of the function is given by the list containing
#'     the parameter estimates along with the 2.5% and 97.5% quantiles
#'     of the posterior distributions. Also, several auxiliary objects are
#'     included in the list which are needed for the [print_outputs()] and
#'     [plot_outputs()] functions.
#'
#' @md
#' @examples
#'
#' library(spatstat)
#' # Prepare the dataset:
#' X = trees_N4
#' x_left = x_left_N4
#' x_right = x_right_N4
#' y_bottom = y_bottom_N4
#' y_top = y_top_N4
#'
#' z_beta = list(refor = cov_refor, slope = cov_slope)
#' z_alpha = list(tmi = cov_tmi, tdensity = cov_tdensity)
#' z_omega = list(slope = cov_slope, reserv = cov_reserv)
#'
#' # Determine the union of rectangles:
#' W = NULL
#' W = owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))
#' if (length(x_left) >= 2) {
#'   for (i in 2:length(x_left)) {
#'     W2 = owin(c(x_left[i], x_right[i]), c(y_bottom[i], y_top[i]))
#'     W = union.owin(W, W2)
#'   }
#' }
#'
#' # Dilated observation window:
#' W_dil = dilation.owin(W, 100)
#'
#'
#' # User-specified hyperparameters for prior distributions:
#' control = list(NStep = 100, BurnIn = 50, SamplingFreq = 5,
#'     Prior_alpha_mean = 3, Prior_alpha_SD = 2, Prior_omega_mean = 5.5,
#'     Prior_omega_SD = 5, Prior_alphavec_SD = c(4.25, 0.012),
#'     Prior_omegavec_SD = c(0.18,0.009))
#'
#' # MCMC estimation:
#' Output = estintp(X, control, x_left, x_right, y_bottom, y_top,
#'     W_dil, z_beta, z_alpha, z_omega, verbose = FALSE)
#'
#' # Text output + series of figures:
#' print_outputs(Output)
#' plot_outputs(Output)
#'
#' @export
#'
estintp <- function(X, control, x_left, x_right, y_bottom, y_top, W_dil,
                    z_beta=NULL, z_alpha=NULL, z_omega=NULL, verbose = TRUE){
  Y = t(rbind(X$x,X$y)) # Y just contains the points of point process and nothing else.

  NStep = control$NStep
  SamplingFreq = control$SamplingFreq
  BurnIn = control$BurnIn

  # Wpix is a binary image of dimensions of W_dil, pixels with "1" belong to W, pixels with "0" do not
  fW <- function(x,y){as.numeric(spatstat.geom::inside.owin(x,y,W))}
  Wpix <- spatstat.geom::as.im(fW, W=W_dil)
  Wpix$v[is.na(Wpix$v)] <- 0

  # WMRpix is a binary image of dimensions of W_dil, pixels with "1" belong to W, pixels with "NA" do not
  fWMR <- function(x,y){as.numeric(spatstat.geom::inside.owin(x,y,W_dil))}
  WMRpix <- spatstat.geom::as.im(fWMR, W=W_dil)

  AreaW = spatstat.geom::area.owin(W)
  AreaMRW = spatstat.geom::area.owin(W_dil)

  if (length(z_alpha)>0){
    for (i in 1:length(z_alpha)){
      z_alpha[[i]]$v[is.na(z_alpha[[i]]$v)] <- 0
    }
  }
  if (length(z_omega)>0){
    for (i in 1:length(z_omega)){
      z_omega[[i]]$v[is.na(z_omega[[i]]$v)] <- 0
    }
  }

  output.first_step = first_step(X,z_beta,W_dil,plot=FALSE)

  pBX_out = pbetXC(Y, z_beta, z_alpha, z_omega, NStep, SamplingFreq, x_left, x_right, y_bottom, y_top, output.first_step, control, AreaW, AreaMRW, W_dil, Wpix, verbose = verbose)
  pBX = pBX_out$pBX

  Postkappa = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),1]
  Postalpha = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),2]
  if(length(z_alpha)>0) Postalphabet = as.matrix(pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),3:(3+length(z_alpha)-1)])
  Postomega = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),3+length(z_alpha)]
  if(length(z_omega)>0) Postomegabet = as.matrix(pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),(4+length(z_alpha)):(4+length(z_alpha)+length(z_omega)-1)])
  PostlogP = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),(4+length(z_alpha)+length(z_omega))]
  if(length(z_beta)>0) Postbeta = as.matrix(pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),(6+length(z_alpha)+length(z_omega)):(dim(pBX)[2]-6)])


  Kappahat = median(Postkappa)
  Alphahat = median(Postalpha)
  if(length(z_alpha)>0) Alphabethat = apply(Postalphabet,2,median) else Alphabethat = NULL
  Omegahat = median(Postomega)
  if(length(z_omega)>0) Omegabethat = apply(Postomegabet,2,median) else Omegabethat = NULL
  if(length(z_beta)>0) Betahat = apply(Postbeta,2,median) else Betahat = NULL

  KappaCI = quantile(Postkappa, probs = c(0.025,0.975))
  AlphaCI = quantile(Postalpha, probs = c(0.025,0.975))
  Q = function(x){quantile(x, probs = c(0.025,0.975))}
  if (length(z_alpha) > 0)
    AlphabetCI = apply(Postalphabet, 2, Q)
  else
    AlphabetCI = NULL
  OmegaCI = quantile(Postomega, probs = c(0.025,0.975))
  if (length(z_omega) > 0)
    OmegabetCI = apply(Postomegabet, 2, Q)
  else
    OmegabetCI = NULL
  if (length(z_beta) > 0 )
    BetaCI = apply(Postbeta, 2, Q)
  else
    BetaCI = NULL

  res = list(Kappahat=Kappahat,
             Alphahat=Alphahat,
             Alphabethat=Alphabethat,
             Omegahat=Omegahat,
             Omegabethat=Omegabethat,
             KappaCI=KappaCI,
             AlphaCI=AlphaCI,
             AlphabetCI=AlphabetCI,
             OmegaCI=OmegaCI,
             OmegabetCI=OmegabetCI,
             Betahat=Betahat,
             BetaCI=BetaCI,
             bet=output.first_step$FS$coef[-1],
             first_step=output.first_step,
             pBX=pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),],
             pBX_backup=pBX,
             priorParameters=pBX_out$priorParameters,
             control=control,
             WMRpix=WMRpix,
             W_dil=W_dil,
             noAlpha=length(z_alpha),
             noOmega=length(z_omega),
             z_beta=z_beta,
             z_alpha=z_alpha,
             z_omega=z_omega
  )

  return(res)
}


#' Re-estimate the posterior distributions with different burn-in
#'
#' @description After running the MCMC chain for the given number of steps, the trace plots may indicate that too small value of burn-in was used in the first place. This function enables re-estimating the posterior distributions with a different value of burn-in, without the need to run the MCMC chain again.
#'
#' @details The output of the main function binspp contains all
#'          the intermediate states of the chain (sampled with the required
#'          frequency) no matter what the original value of burn-in was.
#'          This enables simple and quick re-estimation of the posterior
#'          distributions with either higher or lower value of burn-in
#'          than the one used originally. The output of this function has
#'          the same structure as the output of the main function [estintp()].
#'
#' @param Output list, output of the main function estintp.
#' @param BurnIn new value of burn-in.
#'
#' @return List containing the parameter estimates along
#'         with the 2.5% and 97.5% quantiles of the posterior distributions,
#'         along with auxiliary objects needed for printing and plotting
#'         the outputs.
#'
#' @md
#' @examples
#'
#' library(spatstat)
#' # Prepare the dataset:
#' X = trees_N4
#' x_left = x_left_N4
#' x_right = x_right_N4
#' y_bottom = y_bottom_N4
#' y_top = y_top_N4
#'
#' z_beta = list(refor = cov_refor, slope = cov_slope)
#' z_alpha = list(tmi = cov_tmi, tdensity = cov_tdensity)
#' z_omega = list(slope = cov_slope, reserv = cov_reserv)
#'
#' # Determine the union of rectangles:
#' W = owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))
#' if (length(x_left) >= 2) {
#'   for (i in 2:length(x_left)) {
#'     W2 = owin(c(x_left[i], x_right[i]), c(y_bottom[i], y_top[i]))
#'     W = union.owin(W, W2)
#'   }
#' }
#'
#' # Dilated observation window:
#' W_dil = dilation.owin(W, 100)
#'
#'
#' # Default parameters for prior distributions:
#' control = list(NStep = 100, BurnIn = 50, SamplingFreq = 5)
#'
#' # MCMC estimation:
#' Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil,
#'                  z_beta, z_alpha, z_omega, verbose = FALSE)
#'
#' # Text output + series of figures:
#' print_outputs(Output)
#' plot_outputs(Output)
#'
#' # Recompute the outputs when another value of burn-in is desired,
#' # without running the chain again:
#' Out2 <- re_estimate(Output, BurnIn = 80)
#' print_outputs(Out2)
#' plot_outputs(Out2)
#'
#' @export
re_estimate <- function(Output, BurnIn = 0) {

  pBX = Output$pBX_backup
  SamplingFreq = Output$control$SamplingFreq

  Postkappa = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),1]
  Postalpha = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),2]
  if(length(Output$z_alpha)>0) Postalphabet = as.matrix(pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),3:(3+length(Output$z_alpha)-1)])
  Postomega = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),3+length(Output$z_alpha)]
  if(length(Output$z_omega)>0) Postomegabet = as.matrix(pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),(4+length(Output$z_alpha)):(4+length(Output$z_alpha)+length(Output$z_omega)-1)])
  PostlogP = pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),(4+length(Output$z_alpha)+length(Output$z_omega))]
  if(length(Output$z_beta)>0) Postbeta = as.matrix(pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),(6+length(Output$z_alpha)+length(Output$z_omega)):(dim(pBX)[2]-6)])

  Kappahat = median(Postkappa)
  Alphahat = median(Postalpha)
  if(length(Output$z_alpha)>0) Alphabethat = apply(Postalphabet,2,median) else Alphabethat = NULL
  Omegahat = median(Postomega)
  if(length(Output$z_omega)>0) Omegabethat = apply(Postomegabet,2,median) else Omegabethat = NULL
  if(length(Output$z_beta)>0) Betahat = apply(Postbeta,2,median) else Betahat = NULL

  KappaCI = quantile(Postkappa, probs = c(0.025,0.975))
  AlphaCI = quantile(Postalpha, probs = c(0.025,0.975))
  Q = function(x){quantile(x, probs = c(0.025,0.975))}
  if(length(Output$z_alpha)>0) AlphabetCI = apply(Postalphabet,2,Q) else AlphabetCI = NULL
  OmegaCI = quantile(Postomega, probs = c(0.025,0.975))
  if(length(Output$z_omega)>0) OmegabetCI = apply(Postomegabet,2,Q) else OmegabetCI = NULL
  if(length(Output$z_beta)>0) BetaCI = apply(Postbeta,2,Q) else BetaCI = NULL

  control2 = Output$control
  control2$BurnIn = BurnIn

  res = list(Kappahat=Kappahat,
             Alphahat=Alphahat,
             Alphabethat=Alphabethat,
             Omegahat=Omegahat,
             Omegabethat=Omegabethat,
             KappaCI=KappaCI,
             AlphaCI=AlphaCI,
             AlphabetCI=AlphabetCI,
             OmegaCI=OmegaCI,
             OmegabetCI=OmegabetCI,
             Betahat=Betahat,
             BetaCI=BetaCI,
             bet=Output$first_step$FS$coef[-1],
             first_step=Output$first_step,
             pBX=pBX[(round(BurnIn/SamplingFreq)+1):(dim(pBX)[1]),],
             pBX_backup=pBX,
             priorParameters=Output$priorParameters,
             control=control2,
             WMRpix=Output$WMRpix,
             W_dil=Output$W_dil,
             noAlpha=Output$noAlpha,
             noOmega=Output$noOmega,
             z_beta=Output$z_beta,
             z_alpha=Output$z_alpha,
             z_omega=Output$z_omega
  )

  return(res)
}


#' Text output describing the posterior distributions
#'
#' @description The summaries of the posterior distributions in the text form
#'          are provided.
#'
#' @details The parameter estimates (sample medians
#'          from the empirical posterior distributions) and the 2.5%
#'          and 97.5% quantiles from the empirical posterior
#'          distributions are printed. \cr
#'          Additionally, during the run of the MCMC chain the significance
#'          of the covariates in the list \emph{z_beta} with respect to the
#'          current population of parent points is repeatedly tested.
#'          This function prints the medians of the series of p-values
#'          obtained in this way, together with the corresponding
#'          2.5% and 97.5% sample quantiles of the p-values for each covariate.
#'
#' @param Output list, output of the main function [estintp()].
#'
#' @return Text output summarizing the posterior distributions.
#'
#' @examples
#'
#' library(spatstat)
#' # Prepare the dataset:
#' X = trees_N4
#' x_left = x_left_N4
#' x_right = x_right_N4
#' y_bottom = y_bottom_N4
#' y_top = y_top_N4
#'
#' z_beta = list(refor = cov_refor, slope = cov_slope)
#' z_alpha = list(tmi = cov_tmi, tdensity = cov_tdensity)
#' z_omega = list(slope = cov_slope, reserv = cov_reserv)
#'
#' # Determine the union of rectangles:
#' W = owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))
#' if (length(x_left) >= 2) {
#'   for (i in 2:length(x_left)) {
#'     W2 = owin(c(x_left[i], x_right[i]), c(y_bottom[i], y_top[i]))
#'     W = union.owin(W, W2)
#'   }
#' }
#'
#' # Dilated observation window:
#' W_dil = dilation.owin(W, 100)
#'
#'
#' # Default parameters for prior distributions:
#' control = list(NStep = 100, BurnIn = 20, SamplingFreq = 5)
#'
#'
#' # MCMC estimation:
#' Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil,
#'                  z_beta, z_alpha, z_omega, verbose = FALSE)
#'
#'
#' # Text output
#' print_outputs(Output)
#'
#' @md
#' @export
print_outputs <- function(Output){

  cat("Estimate of kappa (intensity of parent process): ")
  cat(Output$Kappahat)
  cat("\n")

  cat("Estimate of alpha (related to mean number of points in a cluster): ")
  cat(Output$Alphahat)
  cat("\n")

  cat("Estimate of exponential of alpha (related to mean number of points in a cluster): ")
  cat(exp(Output$Alphahat))
  cat("\n")

  if (length(Output$z_alpha)>0){
    cat("Estimate of regression parameters influencing alpha: ")
    cat(Output$Alphabethat)
    cat("\n")
  }

  cat("Estimate of omega (related to standard deviation of cluster size): ")
  cat(Output$Omegahat)
  cat("\n")

  cat("Estimate of exponential of omega (related to standard deviation of cluster size): ")
  cat(exp(Output$Omegahat))
  cat("\n")

  if (length(Output$z_omega)>0){
    cat("Estimate of regression parameters influencing omega: ")
    cat(Output$Omegabethat)
    cat("\n")
  }

    cat("Estimate of regression parameters from the first step, intercept: ")
    cat(Output$first_step$FS$coef[1])
    cat("\n")

  if (length(Output$z_beta)>0){
    cat("Estimate of regression parameters from the first step, covariates: ")
    cat(Output$bet)
    cat("\n")
  }

  if (length(Output$z_beta)>0){
    cat("Estimate of p-values for significance of covariates influencing the intensity function, obtained from the second step: ")
    cat(Output$Betahat)
    cat("\n")
    cat("\n")
  }


  cat("Posterior confidence interval for kappa: \n")
  print(Output$KappaCI)
  cat("\n")

  cat("Posterior confidence interval for alpha: \n")
  print(Output$AlphaCI)
  cat("\n")

  if (length(Output$z_alpha)>0){
    cat("Posterior confidence intervals for regression parameters influencing alpha: \n")
    print(Output$AlphabetCI)
    cat("\n")
  }

  cat("Posterior confidence interval for omega: \n")
  print(Output$OmegaCI)
  cat("\n")

  if (length(Output$z_omega)>0){
    cat("Posterior confidence intervals for regression parameters influencing omega: \n")
    print(Output$OmegabetCI)
    cat("\n")
  }

  if (length(Output$z_beta)>0){
    cat("Posterior confidence intervals for p-value for significance of covariates influencing the intensity function: \n")
    print(Output$BetaCI)
    cat("\n")
  }
}


#' Graphical output describing the posterior distributions
#'
#' @description A graphical representation of the posterior distributions
#'          in terms of histograms and trace plots.
#'
#' @details If the covariate list \emph{z_beta} was non-empty, the estimated
#'          intensity function of the parent process is plotted. Then,
#'          the estimated surface representing the location dependent
#'          mean number of points in a cluster is plotted, and similarly,
#'          the estimated surface representing the location dependent scale
#'          of clusters is plotted. \cr
#'          After that, histograms of the sample posterior distributions
#'          of the individual parameters are plotted, together with
#'          the histograms of p-values giving significance of the individual
#'          covariates in \emph{z_beta} with respect to the population
#'          of parent points. \cr
#'          Then, the trace plots for individual model parameters are plotted,
#'          with highlighted sample median (full red line) and sample 2.5%
#'          and 97.5% quantiles (dashed red lines), and similarly for
#'          the p-values giving significance of the individual covariates
#'          in \emph{z_beta} with respect to the population of parent points.
#'          \cr
#'          Additionally, the following graphs are also plotted:
#' * trace plot for the log-likelihood of the model,
#' * trace plot for the number of parent points,
#' * trace plot for the probability of accepting proposed updates of \eqn{(alpha, alpha_1, …, alpha_l)}, \emph{}
#' * trace plot for the fraction of accepted updates of \eqn{alpha, alpha_1, …, alpha_l} in the last 1000 iterations,
#' * trace plot for the probability of accepting proposed updates of \eqn{omega, omega_1, …, omega_m}, \emph{}
#' * trace plot for the fraction of accepted updates of \eqn{omega, omega_1, …, omega_m} in the last 1000 iterations,
#' * trace plot for the fraction of accepted updates of parent points in the last 1000 iterations.
#'
#' @param Output list, output of the main function [estintp()].
#'
#' @return Series of plots providing a graphical representation
#'         of the posterior distributions in terms of histograms and
#'         trace plots.
#'
#' @examples
#'
#' library(spatstat)
#' # Prepare the dataset:
#' X = trees_N4
#' x_left = x_left_N4
#' x_right = x_right_N4
#' y_bottom = y_bottom_N4
#' y_top = y_top_N4
#'
#' z_beta = list(refor = cov_refor, slope = cov_slope)
#' z_alpha = list(tmi = cov_tmi, tdensity = cov_tdensity)
#' z_omega = list(slope = cov_slope, reserv = cov_reserv)
#'
#' # Determine the union of rectangles:
#' W = owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))
#' if (length(x_left) >= 2) {
#'   for (i in 2:length(x_left)) {
#'     W2 = owin(c(x_left[i], x_right[i]), c(y_bottom[i], y_top[i]))
#'     W = union.owin(W, W2)
#'   }
#' }
#'
#' # Dilated observation window:
#' W_dil = dilation.owin(W, 100)
#'
#'
#' # Default parameters for prior distributions:
#' control = list(NStep = 100, BurnIn = 20, SamplingFreq = 5)
#'
#'
#' # MCMC estimation:
#' Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil,
#'                  z_beta, z_alpha, z_omega, verbose = FALSE)
#'
#' # Text output + series of figures:
#' print_outputs(Output)
#' plot_outputs(Output)
#'
#' @md
#' @export
plot_outputs <- function(Output){

  # Estimated intensity function
  if (length(Output$z_beta)>0){
    bet=Output$first_step$FS$coef[-1]
    B = bet[1]*Output$z_beta[[1]]
    if(length(bet)>1){
      for(i in 2:length(bet)){B=B+bet[i]*Output$z_beta[[i]]}
    }
    B = exp(B) #Covariate part of intensity
    plot(exp(Output$first_step$FS$coef[1])*B, main="Plot of estimated first-order intensity function")
  }

  # Plotting alpha(u):
  if (Output$noAlpha>0){
    fWMR <- function(x,y){as.numeric(spatstat.geom::inside.owin(x,y,Output$W_dil))}
    WMRpix2 <- spatstat.geom::as.im(fWMR, W=Output$W_dil, dimyx=Output$z_alpha[[1]]$dim)

    Balpha = Output$Alphabethat[1]*Output$z_alpha[[1]]
    if(Output$noAlpha>1){
      for(i in 2:Output$noAlpha){Balpha=Balpha+Output$Alphabethat[i]*Output$z_alpha[[i]]}
    }
    Balpha = exp(Output$Alphahat)*exp(Balpha)*WMRpix2
    plot(Balpha, main="Surface of alpha(c)")
  }

  # Plotting omega(u):
  if (Output$noOmega>0){
    fWMR <- function(x,y){as.numeric(spatstat.geom::inside.owin(x,y,Output$W_dil))}
    WMRpix2 <- spatstat.geom::as.im(fWMR, W=Output$W_dil, dimyx=Output$z_omega[[1]]$dim)

    Bomega = Output$Omegabethat[1]*Output$z_omega[[1]]
    if(Output$noOmega>1){
      for(i in 2:Output$noOmega){Bomega=Bomega+Output$Omegabethat[i]*Output$z_omega[[i]]}
    }
    Bomega = exp(Output$Omegahat)*exp(Bomega)*WMRpix2
    plot(Bomega, main="Surface of omega(c)")
  }


  # Plotting posterior distributions:
  Postkappa = Output$pBX[,1]
  Postalpha = Output$pBX[,2]
  if(Output$noAlpha>0) Postalphabet = as.matrix(Output$pBX[,3:(2+Output$noAlpha)]) else Postalphabet = NULL
  Postomega = Output$pBX[,3+Output$noAlpha]
  if(Output$noOmega>0) Postomegabet = as.matrix(Output$pBX[,(4+Output$noAlpha):(3+Output$noAlpha+Output$noOmega)]) else Postomegabet = NULL
  PostlogP = Output$pBX[,4+Output$noAlpha+Output$noOmega]
  if (length(Output$z_beta)>0) Postbeta = as.matrix(Output$pBX[,(6+Output$noAlpha+Output$noOmega):(dim(Output$pBX)[2]-6)]) else Postbeta = NULL
  noParents = Output$pBX[,dim(Output$pBX)[2]-5]
  probsAlpha = Output$pBX[,dim(Output$pBX)[2]-4]
  probsOmega = Output$pBX[,dim(Output$pBX)[2]-3]
  acceptsAlpha = Output$pBX[,dim(Output$pBX)[2]-2]
  acceptsOmega = Output$pBX[,dim(Output$pBX)[2]-1]
  acceptsParents = Output$pBX[,dim(Output$pBX)[2]]

  # Posterior histogram for kappa (Intensity of parent process)
  hist(Postkappa, xlab="kappa", main="Posterior distribution of kappa")
  # Posterior histogram for alpha (Mean number of points in a cluster)
  hist(Postalpha, xlab="alpha", main="Posterior distribution of alpha")
  # Posterior histogram for regression parameters influencing alpha (Mean number of points in a cluster)
  if (Output$noAlpha>0){
    for (i in 1:Output$noAlpha){
      hist(Postalphabet[,i], xlab=paste("alphavec(",i,")", sep=""), main=paste("Posterior distribution of alphavec(",i,")", sep=""))
    }
  }
  # Posterior histogram for omega (standard deviation of cluster size)
  hist(Postomega, xlab="omega", main="Posterior distribution of omega")
  # Posterior histogram for regression parameters influencing omega (standard deviation of cluster size)
  if (Output$noOmega>0){
    for (i in 1:Output$noOmega){
      hist(Postomegabet[,i], xlab=paste("omegavec(",i,")", sep=""), main=paste("Posterior distribution of omegavec(",i,")", sep=""))
    }
  }

  # Posterior histogram for p-value (regression parameters of covariates)
  if (length(Output$z_beta)>0){
    for (i in 1:length(Output$z_beta)){
      hist(Postbeta[,i], xlab="p-value", main=paste("Posterior distribution of p-value for beta(",i,")", sep=""))
    }
  }

  # For plotting nice axis labels
  DS = Output$control$BurnIn
  NS = Output$control$NStep
  labs = c(DS, DS+(NS-DS)*1/4, DS+(NS-DS)*2/4, DS+(NS-DS)*3/4, DS+(NS-DS)*4/4)


  # MCMC trace for kappa
  plot(Postkappa,type = 'l',xaxt = "n",xlab="MCMC step", ylab="kappa", main="MCMC trace for kappa")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)
  abline(h=Output$Kappahat, col="red")
  abline(h=Output$KappaCI[1], col="red", lty=2)
  abline(h=Output$KappaCI[2], col="red", lty=2)

  # MCMC trace for alpha
  plot(Postalpha,type = 'l',xaxt = "n",xlab="MCMC step", ylab="alpha", main="MCMC trace for alpha")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)
  abline(h=Output$Alphahat, col="red")
  abline(h=Output$AlphaCI[1], col="red", lty=2)
  abline(h=Output$AlphaCI[2], col="red", lty=2)

  # MCMC trace for regression parameters influencing alpha
  if(Output$noAlpha>0){
    for(i in 1:Output$noAlpha){
      plot(Postalphabet[,i],type = 'l',xaxt = "n",xlab="MCMC step", ylab=paste("alphavec(",i,")", sep=""), main=paste("MCMC trace for alphavec(",i,")", sep=""))
      axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)
      abline(h=Output$Alphabethat[i], col="red")
      abline(h=Output$AlphabetCI[1,i], col="red", lty=2)
      abline(h=Output$AlphabetCI[2,i], col="red", lty=2)
    }
  }
  # MCMC trace for omega
  plot(Postomega,type = 'l',xaxt = "n", xlab="MCMC step", ylab="omega", main="MCMC trace for omega")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)
  abline(h=Output$Omegahat, col="red")
  abline(h=Output$OmegaCI[1], col="red", lty=2)
  abline(h=Output$OmegaCI[2], col="red", lty=2)

  # MCMC trace for regression parameters influencing omega
  if(Output$noOmega>0){
    for(i in 1:Output$noOmega){
      plot(Postomegabet[,i],type = 'l',xaxt = "n",xlab="MCMC step", ylab=paste("omegavec(",i,")", sep=""), main=paste("MCMC trace for omegavec(",i,")", sep=""))
      axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)
      abline(h=Output$Omegabethat[i], col="red")
      abline(h=Output$OmegabetCI[1,i], col="red", lty=2)
      abline(h=Output$OmegabetCI[2,i], col="red", lty=2)
    }
  }


  # MCMC trace for p-values for beta
  if (length(Output$z_beta)>0){
    for (i in 1:length(Output$z_beta)){
      plot(Postbeta[,i],type = 'l',xaxt = "n",xlab="MCMC step", ylab="p-value", main=paste("MCMC trace for p-value for beta(",i,")", sep=""))
      axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)
      abline(h=Output$Betahat[i], col="red")
      abline(h=Output$BetaCI[1,i], col="red", lty=2)
      abline(h=Output$BetaCI[2,i], col="red", lty=2)
    }
  }


  # MCMC trace for log likelihood
  plot(PostlogP,type = 'l',xaxt = "n", xlab="MCMC step", ylab="log-likelihood", main="MCMC trace for log-likelihood")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)


  # MCMC trace for number of parents
  plot(noParents,type = 'l',xaxt = "n", xlab="MCMC step", ylab="parents", main="MCMC trace for number of parents")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)


  # MCMC trace for probability of accepting proposed updates of alpha and alphabet
  plot(probsAlpha,type = 'l',xaxt = "n", xlab="MCMC step", ylab="probability", main="Probability of accepting update of alpha, alphavec")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)


  # Fraction of accepted updates of alpha, alphabet in the past 1000 steps
  plot(acceptsAlpha,type = 'l',xaxt = "n", xlab="MCMC step", ylab="fraction", main="Fraction of accepted updates of alpha, alphavec in the past 1000 steps")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)


  # MCMC trace for probability of accepting proposed updates of omega and omegabet
  plot(probsOmega,type = 'l',xaxt = "n", xlab="MCMC step", ylab="probability", main="Probability of accepting update of omega, omegavec")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)


  # Fraction of accepted updates of omega, omegabet in the past 1000 steps
  plot(acceptsOmega,type = 'l',xaxt = "n", xlab="MCMC step", ylab="fraction", main="Fraction of accepted updates of omega, omegavec in the past 1000 steps")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)


  # Fraction of accepted updates of parents in the past 1000 steps
  plot(acceptsParents,type = 'l',xaxt = "n", xlab="MCMC step", ylab="fraction", main="Fraction of accepted updates of parents in the past 1000 steps")
  axis(1, at=(0:4)*(length(Postkappa)/4), labels=labs)

}
