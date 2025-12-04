##############################################################################################

#' @noRd

KumulaVsech <- function(CC, omega, x_left, x_right, y_bottom, y_top){
##  require(mvtnorm)
  S = 0;
  # Loop through each region
  for(i in 1:length(x_left)){
    # Function to calculate the probability mass in the current region for a given point
    .dist = function(C){pmvnorm(lower=c(x_left[i],y_bottom[i]), upper=c(x_right[i],y_top[i]), mean=C,sigma=diag(2)*omega^2)}
    S = S + sum(apply(CC, 1, .dist))
  };
  S
};

##############################################################################################

#' @noRd

GenerateCenterPoint <- function(W_dil){P = runifpoint(1,W_dil);c(P$x[1],P$y[1])}

##############################################################################################

#' @noRd

StepMovePoint <- function(kappa, alpha, omega, theta, r, likelihoodprev, rho0sum, X, CC, logP, integral, AreaMRW, W_dil,
                          x_left, x_right, y_bottom, y_top,
                          AreaW){

  if (runif(1) < 1/3) {

    Discard = ceiling(runif(1,min = 0, max = dim(CC)[1])); # Randomly select a point to discard
    int2 = integral - KumulaVsech(t(as.matrix(CC[Discard,])), omega, x_left, x_right, y_bottom, y_top); # Update integral
    Interaction = DeathInteractionLik(CC[Discard,], Discard, matrix(CC,ncol=2), rho0sum, theta, r) # Compute new interaction after death
    likelihood = Interaction$likelihood
    rhosum = Interaction$rhosum
    CC1 = CC[-Discard,]; # Remove the selected point

    NewCenter = GenerateCenterPoint(W_dil);
    Interaction = BirthInteractionLik(NewCenter, matrix(CC1,ncol=2), rhosum, theta, r) # Compute new interaction after birth
    likelihood = Interaction$likelihood
    rhosum = Interaction$rhosum
    CC1 = rbind(CC1, NewCenter); # Add the new point
    int2 = int2 + KumulaVsech(t(as.matrix(NewCenter)), omega, x_left, x_right, y_bottom, y_top);
    logP2 = logpXCbeta(X, CC1, alpha, omega, AreaW, int2);

    if(log(runif(1)) < (logP2 - logP + likelihood - likelihoodprev)) {
      Vystup = list(CC1, logP2, int2, likelihood, rhosum)} else {Vystup = list(CC, logP, integral, likelihoodprev, rho0sum)};

  } else {

    if(runif(1) < 1/2 || dim(CC)[1] < 3){ # Add a new point if random value is less than 1/2 or if there are fewer than 3 points
      NewCenter = GenerateCenterPoint(W_dil);
      Interaction = BirthInteractionLik(NewCenter, matrix(CC,ncol=2), rho0sum, theta, r)
      likelihood = Interaction$likelihood
      rhosum = Interaction$rhosum
      CC1 = rbind(CC, NewCenter);
      int2 = integral + KumulaVsech(t(as.matrix(NewCenter)), omega, x_left, x_right, y_bottom, y_top);

      logP2 = logpXCbeta(X, CC1, alpha, omega, AreaW, int2);
      if(log(runif(1)) < (logP2 - logP + likelihood - likelihoodprev + log(kappa*(AreaMRW)) - log(dim(CC1)[1]))) {
        Vystup = list(CC1, logP2, int2, likelihood, rhosum)} else {Vystup = list(CC, logP, integral, likelihoodprev, rho0sum)};

    } else { # Remove a point if there are more than 2 points

      Discard = ceiling(runif(1,min = 0, max = dim(CC)[1]));
      Interaction = DeathInteractionLik(CC[Discard,], Discard, matrix(CC,ncol=2), rho0sum, theta, r)
      likelihood = Interaction$likelihood
      rhosum = Interaction$rhosum
      int2 = integral - KumulaVsech(t(as.matrix(CC[Discard,])), omega, x_left, x_right, y_bottom, y_top);
      CC1 = CC[-Discard,];
      logP2 = logpXCbeta(X, CC1, alpha, omega, AreaW, int2);

      if(log(runif(1)) < (logP2 - logP + likelihood - likelihoodprev - log(kappa*(AreaMRW)) + log(dim(CC1)[1]))) {
        Vystup = list(CC1, logP2, int2, likelihood, rhosum)} else {Vystup = list(CC, logP, integral, likelihoodprev, rho0sum)};
    }

  }
  return(Vystup)
}

##############################################################################################

#' Generate auxiliary variable for given proposed parameters.
#' @export
#'
#' @description Generate auxiliary variable for given proposed parameters.
#'
#' @param kappa parameter to generate auxiliary variable. \eqn{\kappa} represents the -.
#' @param theta parameter vector
#' \ifelse{html}{
#' \out{&Theta; = (&theta;<sub>1</sub>, &theta;<sub>2</sub>). &theta;<sub>1</sub> represents -. &theta;<sub>2</sub> represents-.}}{
#' \eqn{\Theta = (\theta_1, \theta_2). \theta_1 represents -. \theta_2 represents -.}
#' }.
#' @param likelihoodprev previous likelihood value.
#' @param rho0sum Initial sum of interaction strengths.
#' @param CC cluster centers.
#' @param AreaMRW area of dilated window.
#' @param W_dil observation window dilated by the assumed maximal cluster radius.
#' @param niter number of iterations of MCMC.
#'
#' @return The output is a list of \code{CC}, \code{likelihoodprev}, \code{rho0sum}.
#' @export

AuxVarGen <- function(kappa, theta, likelihoodprev, rho0sum, CC, AreaMRW, W_dil, niter){

  r = coeff(theta)

  for(i in 1:niter){

    if(runif(1) < 1/3){
      Discard = ceiling(runif(1, min = 0, max = dim(CC)[1]));
      Interaction = DeathInteractionLik(CC[Discard,], Discard, matrix(CC,ncol=2), rho0sum, theta, r)
      likelihood = Interaction$likelihood
      rhosum = Interaction$rhosum
      CC1 = CC[-Discard,];

      NewCenter = GenerateCenterPoint(W_dil);
      Interaction = BirthInteractionLik(NewCenter, matrix(CC1,ncol=2), rhosum, theta, r)
      likelihood = Interaction$likelihood
      rhosum = Interaction$rhosum
      CC1 = rbind(CC1, NewCenter);
      if(log(runif(1)) < (likelihood - likelihoodprev)) {
        CC = CC1; likelihoodprev = likelihood; rho0sum = rhosum
      }

    } else {

      if(runif(1) < 1/2 || dim(CC)[1] < 3){
        NewCenter = GenerateCenterPoint(W_dil);
        Interaction = BirthInteractionLik(NewCenter, matrix(CC,ncol=2), rho0sum, theta, r)
        likelihood = Interaction$likelihood
        rhosum = Interaction$rhosum
        CC1 = rbind(CC, NewCenter);
        if(log(runif(1)) < (likelihood - likelihoodprev + log(kappa*(AreaMRW)) - log(dim(CC1)[1]))) {
          CC = CC1; likelihoodprev = likelihood; rho0sum = rhosum
        }

      } else {
        Discard = ceiling(runif(1,min = 0, max = dim(CC)[1]));
        Interaction = DeathInteractionLik(CC[Discard,], Discard, matrix(CC,ncol=2), rho0sum, theta, r)
        likelihood = Interaction$likelihood
        rhosum = Interaction$rhosum
        CC1 = CC[-Discard,];
        if(log(runif(1)) < (likelihood - likelihoodprev - log(kappa*(AreaMRW)) + log(dim(CC1)[1]))) {
          CC = CC1; likelihoodprev = likelihood; rho0sum = rhosum
        }
      }
    }
  }
  Vystup = list(CC, likelihoodprev, rho0sum)
  return(Vystup)
}

##############################################################################################

#' Estimation of interaction Neyman-Scott point process using auxiliary variable algorithm into Markov chain Monte Carlo.
#' @export
#'
#' @import spatstat
#' @import fields
#'
#' @description The Bayesian MCMC estimation of parameters for interaction Neyman-Scott point process using auxiliary variable algorithm into Markov chain Monte Carlo.
#'
#' @param X observed point pattern in the \code{\link[spatstat.geom]{ppp}} format of the \strong{spatstat} package.
#' @param control list specifying various tuning constants for the MCMC estimation. See also Details.
#' @param x_left vector describing the observation window, contains the lower x-coordinate of the corners of each rectangle.
#' @param x_right vector describing the observation window, contains the higher x-coordinate of the corners of each rectangle.
#' @param y_bottom vector describing the observation window, contains the smaller y-coordinate of the corners of each rectangle.
#' @param y_top vector describing the observation window, contains the higher y-coordinate of the corners of each rectangle.
#' @param W_dil observation window dilated by the assumed maximal cluster radius.
#' @param AreaW area of a window.
#' @param AreaMRW area of a dilated window.
#' @param radius The radius of dilation, passed as an argument to \code{\link[spatstat.geom]{dilation.owin}} to create a dilated window.
#'
#' @return The output of the function is given by the list containing the parameter estimates along with the 2.5\% and 97.5\% confidence interval of the posterior distributions, cluster centers, \emph{control}, \code{W_dil} and parameters.
#'
#' @details
#' \strong{\emph{Observation window and its dilation}}
#'
#' The observation window must be provided as the union of aligned rectangles, aligned with the coordinate axes.
#' This, however, allows the analysis of point patterns observed in rather irregular regions by approximating the region by a union of aligned rectangles.
#' The structure of the vectors \emph{x_left, x_right, y_bottom} and \emph{y_top} is such that the first rectangle is constructed using the function \code{\link[spatstat.geom]{owin}} from the spatstat package as \code{owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))}, and similarly for the other rectangles.
#' Naturally, a rectangular window can be used and in such a case the vectors \emph{x_left} to \emph{y_top} each contain a single element.
#'
#' \strong{\emph{control}}
#'
#' The control list must contain the following elements: \strong{\emph{NStep}} (the required number of MCMC iterations to be performed),
#' \strong{\emph{BurnIn}} (burn-in, how many iterations at the beginning of the chain will be disregarded when computing the resulting estimates),
#' \strong{\emph{SamplingFreq}} (sampling frequency for estimating the posterior distributions). Additionally, the hyperparameters for the prior distributions should be given, see below.
#'
#' \strong{\emph{Prior distributions and hyperparameters}}
#'
#' The prior distribution for each parameter (\emph{alpha}, \emph{omega}, \emph{kappa}, \emph{theta1}, \emph{theta2}) is normal with respective \code{mean=Prior_<parameter>_mean} and \code{SD=Prior_<parameter>_SD}.
#' During update, the Lower Bound for each parameter is \code{<parameter>_LB} and the Upper Bound is \code{<parameter>_UB}.
#' \emph{alpha} controls the expected number of occurrences or events around specified focal points or regions, \emph{omega} controls the width of cluster events activity and \emph{kappa} controls the overall intensity of the parent process.
#' \emph{theta1} controls the shape of the interaction function and \emph{theta2} controls the shape of the interaction function and provides the location of the peak value.
#'
#' @references Park, J., Chang, W., & Choi, B. (2022). An interaction Neyman-Scott point process model for coronavirus disease-19. \emph{Spatial Statistics}, 47, 100561.
#' @seealso \code{\link{print.output_estinternsp}}
#'
#' @examples
#' library(spatstat)
#' # library(spatstat.geom)
#' library(fields)
#' library(mvtnorm)
#' library(binspp)
#'
#' # Generate example window
#' W <- owin(xrange = c(0, 100), yrange = c(0, 50)) # example window (rectangle)
#' radius = 2
#' W_dil = dilation.owin(W, radius)
#' AreaW <- area(W)
#' AreaMRW <- area(W_dil)
#'
#' x_left = W$xrange[1]
#' x_right = W$xrange[2]
#' y_bottom = W$yrange[1]
#' y_top = W$yrange[2]
#'
#' # True parameters
#' alpha <- 15
#' omega <- 1.5
#' kappa <- 0.001
#' theta1 <- 3
#' theta2 <- 10
#'
#' # Generate parents process
#' CCC <- rpoispp(kappa, win = W_dil)
#' CC <- t(rbind(CCC$x,CCC$y))
#' CCdist <- rdist(CC)
#' r <- binspp::coeff(c(theta1,theta2))
#' r1 <- r[1]; r2 <- r[2]; t1 <- theta1; t2 <- theta2; t3 <- 0.5; R <- 0;
#' rho0sum <- rep(0,dim(CC)[1])                # row sum of rho matrix
#' res <- binspp::pCClik2(c(theta1,theta2), CC)
#' rho0sum <- res$rhosum
#' likelihoodprev <- res$likelihood
#'
#' # this is just for example, for a real estimate
#' # use value of at least niter = 10000
#' CC.true <- AuxVarGen(kappa, c(theta1,theta2), likelihoodprev, rho0sum, CC,
#'     AreaMRW, W_dil, 100)
#'
#' # Generate offsprings given parents
#' gaus <- function(n, omega) {
#'   matrix(rnorm(2 * n, mean=0, sd=omega), ncol=2)
#' }
#' parents <- ppp(CC.true[[1]][,1],CC.true[[1]][,2],window=W)
#' np <- npoints(parents)
#'
#' csize <- qpois(runif(np, min = dpois(0, alpha)), alpha)
#' noff <- sum(csize)
#' xparent <- parents$x
#' yparent <- parents$y
#' x0 <- rep.int(xparent, csize)
#' y0 <- rep.int(yparent, csize)
#' dd <- gaus(noff, omega)
#' xy <- xy.coords(dd)
#' dx <- xy$x; dy <- xy$y
#'
#' xoff <- x0 + dx; yoff <- y0 + dy
#' result <- ppp(xoff, yoff, window = W_dil, check = FALSE, marks = NULL)
#' X <- cbind(result$x,result$y) # Generate example data
#'
#' # this is just for example, for a real estimate
#' # use values of at least NStep = 20000 and BurnIn = 5000
#' control = list(NStep = 50, BurnIn = 25, SamplingFreq = 10,
#'    Prior_alpha_mean = 15, Prior_alpha_SD = 2, alpha_LB = 10, alpha_UB = 20,
#'    Prior_omega_mean = 1.5, Prior_omega_SD = 0.2, omega_LB = 1, omega_UB = 2,
#'    Prior_kappa_mean = 0.001, Prior_kappa_SD = 0.0001, kappa_LB = 0.0005,
#'                              kappa_UB = 0.002,
#'    Prior_theta1_mean = 3, Prior_theta1_SD = 0.2, theta1_LB = 2.5,
#'                              theta1_UB = 3.5,
#'    Prior_theta2_mean = 10, Prior_theta2_SD = 1, theta2_LB = 8,
#'                              theta2_UB = 12)
#'
#' Output = estinternsp(X, control,
#'                      x_left, x_right, y_bottom, y_top,
#'                      W_dil,
#'                      AreaW, AreaMRW, radius)
#'
#' @export

estinternsp <- function(X, control,
                        x_left, x_right, y_bottom, y_top,
                        W_dil,
                        AreaW, AreaMRW, radius){

  ### Initializing part

  NStep = control$NStep
  SamplingFreq = control$SamplingFreq
  BurnIn = control$BurnIn

  parameter <- matrix(0,NStep,5)
  CCdim <- rep(0,NStep)

  alpha <- control$Prior_alpha_mean
  alpha_sd <- control$Prior_alpha_SD
  alpha_lb <- control$alpha_LB # update restriction of alpha
  alpha_ub <- control$alpha_UB
  omega <- control$Prior_omega_mean
  omega_sd <- control$Prior_omega_SD
  omega_lb <- control$omega_LB # update restriction of omega
  omega_ub <- control$omega_UB
  kappa <- control$Prior_kappa_mean
  kappa_sd <- control$Prior_kappa_SD
  kappa_lb <- control$kappa_LB # update restriction of kappa
  kappa_ub <- control$kappa_UB
  theta1 <- control$Prior_theta1_mean
  theta1_sd <- control$Prior_theta1_SD
  theta1_lb <- control$theta1_LB # update restriction of theta1
  theta1_ub <- control$theta1_UB
  theta2 <- control$Prior_theta2_mean
  theta2_sd <- control$Prior_theta2_SD
  theta2_lb <- control$theta2_LB # update restriction of theta2
  theta2_ub <- control$theta2_UB

  set.seed(1)
  CCC <- rpoispp(kappa, win = W_dil)
  CC <- t(rbind(CCC$x,CCC$y));CCdim[1] <- dim(CC)[1]
  integral <- KumulaVsech(CC, omega, x_left, x_right, y_bottom, y_top)
  #logP <- logpXCbeta(X,CC,alpha,omega,area(W),integral)
  logP <- logpXCbeta(X, CC, alpha, omega, AreaW, integral)

  # interaction loglikelihood
  CCdist <- rdist(CC)
  r <- coeff(c(theta1,theta2));
  r1 <- r[1]; r2 <- r[2]; t1 <- theta1; t2 <- theta2; t3 <- 0.5; R <- 0;
  rho0sum <- rep(0,dim(CC)[1]);                # row sum of rho matrix
  res <- pCClik2(c(theta1,theta2), CC)
  rho0sum <- res$rhosum
  likelihoodprev <- res$likelihood

  parameter[1,] <- c(alpha,omega,kappa,theta1,theta2)



  ### Markov chain Monte Carlo

  for(step in 2:NStep){

    if (step%%100==0) {print(step)}

    ## Update offspring parameters
    Newalpha <- rnorm(1,alpha,alpha_sd)
    Newomega <- rnorm(1,omega,omega_sd)

    # update likelihoods
    if( Newalpha < alpha_lb | Newalpha > alpha_ub | Newomega < omega_lb | Newomega > omega_ub ){ logratio=-Inf }else{
      int2 <- KumulaVsech(CC, Newomega, x_left, x_right, y_bottom, y_top);
      logP2 <- logpXCbeta(X, CC, Newalpha, Newomega, AreaW, int2);

      # update parameters with following probabilities
      logratio <- logP2-logP
    }
    if(log(runif(1))<logratio){
      alpha <- Newalpha;
      omega <- Newomega;
      logP <- logP2; integral = int2;
    }


    ## Update parents parameters (DMH update part)
    Newkappa <- rnorm(1,kappa,kappa_sd)
    Newtheta1 <- rnorm(1,theta1,theta1_sd)
    Newtheta2 <- rnorm(1,theta2,theta2_sd)

    if( Newkappa < kappa_lb | Newkappa > kappa_ub | Newtheta1 < theta1_lb | Newtheta1 > theta1_ub | Newtheta2 < theta2_lb | Newtheta2 > theta2_ub ){ logratio=-Inf }else{

      # generate auxiliary variable for given proposed parameters
      aux <- AuxVarGen(Newkappa, c(Newtheta1,Newtheta2), likelihoodprev, rho0sum, CC, AreaMRW, W_dil, radius)
      logAux2 <- pCClik(c(Newtheta1,Newtheta2), aux[[1]])
      logAux <- pCClik(c(theta1,theta2), aux[[1]])
      logCC2res <- pCClik2(c(Newtheta1,Newtheta2), CC)
      logCC2 <- logCC2res$likelihood
      logCC <- pCClik(c(theta1,theta2), CC)

      # DMH ratio
      logratio <- (logCC2 + logAux) - (logCC + logAux2) + dim(CC)[1]*log(Newkappa/kappa) + dim(aux[[1]])[1]*log(kappa/Newkappa)
    }
    if(log(runif(1))<logratio){
      kappa <- Newkappa
      theta1 <- Newtheta1
      theta2 <- Newtheta2
      likelihoodprev <- logCC2
      rho0sum <- logCC2res$rhosum
    }

    # update parent process using BDMCMC
    r <- coeff(c(theta1,theta2))
    BDMCMCres <- StepMovePoint(kappa, alpha, omega, c(theta1,theta2), r, likelihoodprev, rho0sum, X, CC, logP, integral, AreaMRW, W_dil,
                               x_left, x_right, y_bottom, y_top,
                               AreaW)
    CC <- BDMCMCres[[1]];
    logP <- BDMCMCres[[2]]
    integral <- BDMCMCres[[3]]
    CCdim[step] <- dim(CC)[1]
    likelihoodprev <- BDMCMCres[[4]]
    rho0sum <- BDMCMCres[[5]]

    # save parameter results for each step
    parameter[step,] <- c(alpha,omega,kappa,theta1,theta2)
  }

  # Burn-in and thinning from MCMC results
  Postalpha = parameter[seq(BurnIn + 1, nrow(parameter), by = SamplingFreq), 1]
  Postomega = parameter[seq(BurnIn + 1, nrow(parameter), by = SamplingFreq), 2]
  Postkappa = parameter[seq(BurnIn + 1, nrow(parameter), by = SamplingFreq), 3]
  Posttheta1 = parameter[seq(BurnIn + 1, nrow(parameter), by = SamplingFreq), 4]
  Posttheta2 = parameter[seq(BurnIn + 1, nrow(parameter), by = SamplingFreq), 5]

  # Get parameter's hat value
  Alphahat = median(Postalpha)
  Omegahat = median(Postomega)
  Kappahat = median(Postkappa)
  Theta1hat = median(Posttheta1)
  Theta2hat = median(Posttheta2)

  # Get CI
  AlphaCI = quantile(Postalpha, probs = c(0.025, 0.975))
  OmegaCI = quantile(Postomega, probs = c(0.025, 0.975))
  KappaCI = quantile(Postkappa, probs = c(0.025, 0.975))
  Theta1CI = quantile(Posttheta1, probs = c(0.025, 0.975))
  Theta2CI = quantile(Posttheta2, probs = c(0.025, 0.975))

  # save results
  res = list(CC = CC,
             Alphahat = Alphahat, Omegahat = Omegahat, Kappahat = Kappahat, Theta1hat = Theta1hat, Theta2hat = Theta2hat,
             AlphaCI = AlphaCI, OmegaCI = OmegaCI, KappaCI = KappaCI, Theta1CI = Theta1CI, Theta2CI = Theta2CI,
             control = control, W_dil = W_dil,
             parameter = parameter)
  class(res) <- c("output_estinternsp", class(res))
  return(res)
}

##############################################################################################

#' Text output describing the posterior distributions
#' @export
#'
#' @description Output printing function of interaction neyman-scott process. It prints the estimated values and the posterior confidence intervals for each parameter.
#'
#' @param x Output of the main function \code{\link{estinternsp}}.
#' @param ... further arguments passed from generic print function.
#'
#' @return Text output summarizing the posterior distributions for each parameter.
#'
#' @seealso \code{\link{estinternsp}}
#' @examples
#' library(spatstat)
#' # library(spatstat.geom)
#' library(fields)
#' library(mvtnorm)
#' library(binspp)
#'
#' # Generate example window
#' W <- owin(xrange = c(0, 100), yrange = c(0, 50)) # example window (rectangle)
#' radius = 2
#' W_dil = dilation.owin(W, radius)
#' AreaW <- area(W)
#' AreaMRW <- area(W_dil)
#'
#' x_left = W$xrange[1]
#' x_right = W$xrange[2]
#' y_bottom = W$yrange[1]
#' y_top = W$yrange[2]
#'
#' # True parameters
#' alpha <- 15
#' omega <- 1.5
#' kappa <- 0.001
#' theta1 <- 3
#' theta2 <- 10
#'
#' # Generate parents process
#' CCC <- rpoispp(kappa, win = W_dil)
#' CC <- t(rbind(CCC$x,CCC$y))
#' CCdist <- rdist(CC)
#' r <- binspp::coeff(c(theta1,theta2))
#' r1 <- r[1]; r2 <- r[2]; t1 <- theta1; t2 <- theta2; t3 <- 0.5; R <- 0;
#' rho0sum <- rep(0,dim(CC)[1])                # row sum of rho matrix
#' res <- binspp::pCClik2(c(theta1,theta2), CC)
#' rho0sum <- res$rhosum
#' likelihoodprev <- res$likelihood
#'
#' # this is just for example, for a real estimate
#' # use value of at least niter = 10000
#' CC.true <- AuxVarGen(kappa, c(theta1,theta2), likelihoodprev, rho0sum, CC,
#'     AreaMRW, W_dil, 100)
#'
#' # Generate offsprings given parents
#' gaus <- function(n, omega) {
#'   matrix(rnorm(2 * n, mean=0, sd=omega), ncol=2)
#' }
#' parents <- ppp(CC.true[[1]][,1],CC.true[[1]][,2],window=W)
#' np <- npoints(parents)
#'
#' csize <- qpois(runif(np, min = dpois(0, alpha)), alpha)
#' noff <- sum(csize)
#' xparent <- parents$x
#' yparent <- parents$y
#' x0 <- rep.int(xparent, csize)
#' y0 <- rep.int(yparent, csize)
#' dd <- gaus(noff, omega)
#' xy <- xy.coords(dd)
#' dx <- xy$x; dy <- xy$y
#'
#' xoff <- x0 + dx; yoff <- y0 + dy
#' result <- ppp(xoff, yoff, window = W_dil, check = FALSE, marks = NULL)
#' X <- cbind(result$x,result$y) # Generate example data
#'
#' # this is just for example, for a real estimate
#' # use values of at least NStep = 20000 and BurnIn = 5000
#' control = list(NStep = 50, BurnIn = 25, SamplingFreq = 10,
#'    Prior_alpha_mean = 15, Prior_alpha_SD = 2, alpha_LB = 10, alpha_UB = 20,
#'    Prior_omega_mean = 1.5, Prior_omega_SD = 0.2, omega_LB = 1, omega_UB = 2,
#'    Prior_kappa_mean = 0.001, Prior_kappa_SD = 0.0001, kappa_LB = 0.0005,
#'                              kappa_UB = 0.002,
#'    Prior_theta1_mean = 3, Prior_theta1_SD = 0.2, theta1_LB = 2.5,
#'                              theta1_UB = 3.5,
#'    Prior_theta2_mean = 10, Prior_theta2_SD = 1, theta2_LB = 8,
#'                              theta2_UB = 12)
#'
#' Output = estinternsp(X, control,
#'                      x_left, x_right, y_bottom, y_top,
#'                      W_dil,
#'                      AreaW, AreaMRW, radius)
#' print(Output)
#'
#' @export
print.output_estinternsp <- function(x, ...) {
#print_outputs_internsp <- function(Output){
  Output = x

  cat("Estimate of alpha (controls the expected number of occurrences or events around specified focal points or regions): ")
  cat(Output$Alphahat)
  cat("\n")
  cat("Estimate of omega (controls the width of cluster events activity): ")
  cat(Output$Omegahat)
  cat("\n")
  cat("Estimate of kappa (controls the overall intensity of the parent process): ")
  cat(Output$Kappahat)
  cat("\n")
  cat("Estimate of theta1 (control the shape of the interaction function, provides the peak value): ")
  cat(Output$Theta1hat)
  cat("\n")
  cat("Estimate of theta2 (control the shape of the interaction function, provides the location of the peak value): ")
  cat(Output$Theta2hat)
  cat("\n")

  cat("Posterior confidence interval for alpha: \n")
  print(Output$AlphaCI)
  cat("\n")
  cat("Posterior confidence interval for omega: \n")
  print(Output$OmegaCI)
  cat("\n")
  cat("Posterior confidence interval for kappa: \n")
  print(Output$KappaCI)
  cat("\n")
  cat("Posterior confidence interval for theta1: \n")
  print(Output$Theta1CI)
  cat("\n")
  cat("Posterior confidence interval for theta2: \n")
  print(Output$Theta2CI)
  cat("\n")
}
