#' esthtpr
#'
#' @description Bayesian MCMC estimation of homogeneous Thomas process.
#'
#' This package function introduces the Bayesian MCMC estimate of parameters
#'  in homogeneous Thomas process. The process is set in observation window
#'   W which is a union of rectangles.
#'
#' Parameters of the model:
#'
#'   kappa – intensity of cluster centers
#'
#'   alpha – expected number of points in a cluster
#'
#'   omega – standard deviation of the bivariate normal which forms the spread of points in a cluster
#'
#' The estimation algorithm is the modified Bayesian MCMC algorithm described in Kopecký, Mrkvička (2016).
#'
#' @param A1 contains a list of left bottom corners
#' @param A2 contains a list of right bottom corners
#' @param B1 contains a list of left bottom corners
#' @param B2 contains a list of right bottom corners
#' @param NStep length of MCMC (recommended at least 50000 for homogeneous Thomas process)
#' @param DiscardStep discard phase of MCMC (Recommended at least 10000 for homogeneous Thomas process)
#' @param Jump record every Jump step in order to avoid autocorrelation between steps of MCMC
#' @param MR maximal range of cluster sizes we expect. Used to simulate extra centers outside of study region (square root of (content(W))/10)
#' @param PPalpha1 parameters of prior density for alpha - it is Gamma distribution (shape)
#' @param PPalpha2 the expected number of points per cluster (scale)
#' @param PPomega1 parameter of prior for omega - it is Gamma distribution  (shape)
#' @param PPomega2 the expected radii of the cluster (scale)
#' @param salpha standard deviation for update of alpha
#' @param somega standard deviation for update of omega
#' @param Silent If the program should print the results during the run
#'
#' @return
#' Kappahat – estimate of Kappa – median of aposterior distribution.
#' Alphahat – estimate of Alpha – median of aposterior distribution.
#' Omegahat – estimate of Omega – median of aposterior distribution.
#' KappaCI – 95\% confidence interval for kappa computed from aposterior distribution.
#' AlphaCI – 95\% confidence interval for alpha computed from aposterior distribution.
#' OmegaCI – 95\% confidence interval for omega computed from aposterior distribution.
#' Postkappa – the MCMC chain recorded every Jump step with discarded first Discard steps.
#' Postalpha – the MCMC chain recorded every Jump step with discarded first Discard steps.
#' Postomega – the MCMC chain recorded every Jump step with discarded first Discard steps.
#' PostLogP – the Log likelihood of MCMC chain recorded every Jump step with discarded first Discard steps.

#' @export
#' @importFrom graphics axis hist lines plot
#' @importFrom stats dgamma median pnorm quantile rnorm runif
#'
#' @examples
#'
#'# Specify the observation Window as a list of nonintersecting rectangles
#'# A1 contains a list of left bottom corners, A2 right bottom corners
#'# A1,A2,B1,B2 input - type - vectors
#'
#' A1 = c(0, 1); A2 = c(1, 2); B1 = c(0, 0); B2 = c(1, 0.5);
#'
#'#settings must be specified
#'NStep = 250      # length of MCMC (Recommended at least 50000 for homogeneous Thomas process)
#'DiscardStep = 50 # discard phase of MCMC (Recommended at least 10000 for homogeneous Thomas process)
#'Jump = 10        # record every Jump step in order to avoid autocorrelation between steps of MCMC
#'MR = 100         # Maximal range of cluster sizes we expect. Used to simulate extra centers outside of study region
#'PPalpha1 = 2     # parameters of prior density for alpha - it is Gamma distribution  - shape and scale
#'PPalpha2 = 15    # PPalha2 is the expected number of points per cluster
#'PPomega1 = 2     # parameter of prior for omega - it is Gamma distribution  - shape and scale
#'PPomega2 = 100   # PPomega2 is the expected radii of the cluster
#'salpha = PPalpha2 / 30  # standard deviation for update of alpha
#'somega = PPomega2 / 30  # standard deviation for update of omega
#'Silent = FALSE    # If the program should print the results during the run
#'
#'
#'
#'result = esthtpr(
#'  A1,               # Specify the observation Window as a list of nonintersecting rectangles
#'  A2,               # A1 contains a list of left bottom corners, A2 right bottom corners
#'  B1,
#'  B2,
#'  NStep = NStep,      # length of MCMC (Recommended at least 50000 for homogeneous Thomas process)
#'  DiscardStep = DiscardStep, # discard phase of MCMC (Recommended at least 10000 for homogeneous Thomas process)
#'  Jump = Jump,        # record every Jump step in order to avoid autocorrelation between steps of MCMC
#'  MR = MR,         # Maximal range of cluster sizes we expect. Used to simulate extra centers outside of study region
#'  PPalpha1 = PPalpha1,     # parameters of prior density for alpha - it is Gamma distribution  - shape and scale
#'  PPalpha2 = PPalpha2,    # PPalha2 is the expected number of points per cluster
#'  PPomega1 = PPomega1,     # parameter of prior for omega - it is Gamma distribution  - shape and scale
#'  PPomega2 = PPomega2,   # PPomega2 is the expected radii of the cluster
#'  salpha = salpha,  # standard deviation for update of alpha
#'  somega = somega,  # standard deviation for update of omega
#'  Silent = Silent    # If the program should print the results during the run
#')
#'
#'result
#'
#' @references
#' Kopecký J., Mrkvička T.: On the Bayesian estimation for the stationary Neyman-Scott point processes, Applications of Mathematics 61/4, 2016, 503-514.
#'
esthtpr <- function(
  A1,                  # Specify the observation Window as a list of nonintersecting rectangles
  A2,                  # A1 contains a list of left bottom corners, A2 right bottom corners
  B1,
  B2,
  NStep = 250000,      # length of MCMC (Recommended at least 50000 for homogeneous Thomas process)
  DiscardStep = 50000, # discard phase of MCMC (Recommended at least 10000 for homogeneous Thomas process)
  Jump = 10,           # record every Jump step in order to avoid autocorrelation between steps of MCMC
  MR = NULL,           # Maximal range of cluster sizes we expect. Used to simulate extra centers outside of study region (default: sqrt(area of observation window) / 10)
  PPalpha1 = 2,        # parameters of prior density for alpha - it is Gamma distribution (shape)
  PPalpha2 = 15,       # PPalpha2 is the expected number of points per cluster (scale)
  PPomega1 = 2,        # parameter of prior for omega - it is Gamma distribution  (shape)
  PPomega2 = NULL,     # PPomega2 is the expected radii of the cluster (scale) (default: sqrt(area of observation window) / 10)
  salpha = PPalpha2 / 30,  # standard deviation for update of alpha
  somega = NULL,       # standard deviation for update of omega (default: PPomega2 / 30)
  Silent = FALSE       # If the program should print the results during the run
) {

  # Specify the observation Window as a list of nonintersecting rectangles
  # A1 contains a list of left bottom corners, A2 right bottom corners
  # A1,A2,B1,B2 input - type - vectors

  # I compute the union of rectangles
  W = owin(c(A1[1], A2[1]), c(B1[1], B2[1]))
  if (length(A1) >= 2){
    for(i in 2:length(A1)){
      W2 = owin(c(A1[i], A2[i]), c(B1[i], B2[i]))
      W = union.owin(W, W2)
    }
  }

  #I precalculate some constants
  AreaW = area(W)
  if (is.null(MR)){
    MR = sqrt(AreaW) / 10
  }
  WMR = dilation.owin(W, MR)
  AreaMRW = area(WMR)

  if (is.null(PPomega2)){
    PPomega2 = sqrt(AreaW) / 10
  }

  if (is.null(somega)){
    somega = PPomega2 / 30
  }

  # import data here X must be the data as matrix of two collumns only
  X = rThomas(kappa = 0.00010, scale = 30, mu = 10, win = W)
#  plot(X)
  X = t(rbind(X$x, X$y))


# import data here X must be the data as matrix of two collumns only
  plot(X, main = "data X")

#  print("PPalpha1")
#  print(PPalpha1)
#  print("PPalpha2")
#  print(PPalpha2)

#Prior for alpha
Prioralpha = function(a){dgamma(a, shape = PPalpha1, scale = PPalpha2)};
#plot(Prioralpha(0:(8 * PPalpha2)), type = 'l', main = "Prior for alpha")
plot(Prioralpha(seq(0, 8 * PPalpha2, length.out = 100)), type = 'l', main = "Prior for alpha")

#Prior for omega
Prioromega = function(o){dgamma(o, shape = PPomega1, scale = PPomega2)};
plot(Prioromega(seq(0, 8 * PPomega2, length.out = 100)), type='l', main = "Prior for omega")


##############Bayesian estimation homogeneous Thomas process ##################
#omega is scale parameter here


#inner functions
lll = function(x, CC, omega){
  .prod = function(C){exp(-crossprod((x - C), (x - C)) / (2 * omega ^ 2))}
  sum(apply(CC, 1, .prod))
};
#x=X[1,]
#lll(x,CC,omega)

KumulaVsech = function(CC, omega, A1, A2, B1, B2){
  S=0;
  for(i in 1:length(A1)){
    .dist = function(C){pmvnorm(lower = c(A1[i], B1[i]), upper = c(A2[i], B2[i]), mean = C, sigma = diag(2) * omega ^ 2)}
    S = S + sum(apply(CC, 1, .dist))
  };
  S
};
#KumulaVsech(CC,omega,A1,A2,B1,B2)

logpXCBeta = function(X, CC, alpha, omega, AreaW, integral){
  .lll = function(x){lll(x,CC,omega)}
  AreaW - alpha*integral + dim(X)[1]*log(alpha/(2*pi*omega^2)) + sum(log(apply(X,1,.lll)))
};
#logpXCBeta(X,CC,1,omega,W, 10)

NewCenterPoint = function(WMR){P = runifpoint(1,WMR);c(P$x[1],P$y[1])}
#NewCenterPoint(WMR)

#Integralrhohat = function(kappa,W){kappa*AreaW}
#IntegralrhohatMR = function(kappa,WMR){kappa*(AreaMRW)}
#Integralrhohat(2,W)
#IntegralrhohatMR(2,WMR)


StepBeta = function(kappa, alpha, omega, X, CC, logP, integral, AreaMRW, PPalpha2, PPomega2, salpha, somega){

  while((Newalpha <- stats::rnorm(1, alpha, salpha)) < 0){};
  while((Newomega <- stats::rnorm(1, omega, somega)) < 0){};

  Newkappa = dim(X)[1]/AreaW/Newalpha;

  int2 = KumulaVsech(CC, Newomega, A1, A2, B1, B2);
  logP1 = logpXCBeta(X, CC, Newalpha, Newomega, AreaW, int2);
  #If[P1 == 0., P1 = 0];

  if(log(runif(1)) < (logP1 - logP + kappa*AreaMRW - Newkappa*AreaMRW +
                      dim(CC)[1]*log((alpha/Newalpha)) +
                      log(Prioralpha(Newalpha)*Prioromega(Newomega)) -
                      log(Prioralpha(alpha)*Prioromega(omega)) +
                      log((1 - pnorm(0,alpha, salpha))/(1-pnorm(0,Newalpha, salpha))*(1 -
                                                                                      pnorm(0,omega, somega))/(1-pnorm(0,Newomega, somega)))))
  {Vystup = c(Newkappa, Newalpha, Newomega, logP1, int2)}
  else {Vystup = c(kappa, alpha, omega, logP, integral)}

  Vystup
};
#integral = KumulaVsech(CC,omega,A1,A2,B1,B2)
#logP = logpXCBeta(X,CC,1,omega,AreaW, integral)
#StepBeta(kappa, alpha, omega, X, CC, logP, integral, AreaMRW)



StepMovePoint = function(kappa, alpha, omega, X, CC, logP, integral, AreaMRW, WMR){
  if(runif(1) < 1/3){
    Discard = ceiling(runif(1,min = 0, max = dim(CC)[1]));
    int2 = integral - KumulaVsech(t(as.matrix(CC[Discard,])), omega, A1, A2, B1, B2);
    CC1 = CC[-Discard,];
    NewCenter = NewCenterPoint(WMR);
    CC1 = rbind(CC1, NewCenter);
    int2 = int2 + KumulaVsech(t(as.matrix(NewCenter)), omega, A1, A2, B1, B2);
    logP2 = logpXCBeta(X, CC1, alpha, omega, AreaW, int2);
    if(log(runif(1)) < (logP2 - logP)) {Vystup = list(CC1, logP2, int2)} else {Vystup = list(CC, logP, integral)};
  } else {
    if(runif(1) < 1/2 || dim(CC)[1] < 3){
      NewCenter = NewCenterPoint(WMR);
      CC1 = rbind(CC, NewCenter);
      int2 = integral + KumulaVsech(t(as.matrix(NewCenter)), omega, A1, A2, B1, B2);
      logP2 = logpXCBeta(X, CC1, alpha, omega, AreaW, int2);
      if(log(runif(1)) < (logP2 - logP + log(kappa*(AreaMRW)) - log(dim(CC1)[1]))) {Vystup = list(CC1, logP2, int2)} else {Vystup = list(CC, logP, integral)};
    } else
    {  Discard = ceiling(runif(1,min = 0, max = dim(CC)[1]));
    int2 = integral - KumulaVsech(t(as.matrix(CC[Discard,])), omega, A1, A2, B1, B2);
    CC1 = CC[-Discard,];
    logP2 = logpXCBeta(X, CC1, alpha, omega, AreaW, int2);
    if(log(runif(1)) < (logP2 - logP - log(kappa*(AreaMRW)) + log(dim(CC1)[1]))) {Vystup = list(CC1, logP2, int2)} else {Vystup = list(CC, logP, integral)};
    }
  }
  Vystup
};
#StepMovePoint(kappa, alpha, omega, X, CC, logP, integral, AreaMRW, WMR)


pBetaX = function(X, NStep, Jump, A1, A2, B1, B2, AreaW, AreaMRW, WMR, PPalpha2, PPomega2, salpha, somega) {
  alpha = runif(1,sqrt(dim(X)[1])/10, sqrt(dim(X)[1]));
  omega = runif(1,sqrt(AreaMRW)/100, sqrt(AreaMRW)/10);
  lambda = dim(X)[1]/AreaW;
  kappa = lambda/alpha;
  CC = rpoispp(kappa, win = WMR);
  CC = t(rbind(CC$x,CC$y))
  integral = KumulaVsech(CC,omega,A1,A2,B1,B2)
  logP = logpXCBeta(X,CC,alpha,omega,AreaW, integral)

  pBX = c(kappa, alpha, omega, logP, integral);
  for(step in 1:NStep){
    S = StepBeta(kappa, alpha, omega, X, CC, logP, integral, AreaMRW, PPalpha2, PPomega2, salpha, somega);
    kappa = S[1];
    alpha = S[2];
    omega = S[3];
    logP = S[4];
    integral = S[5];

    T = StepMovePoint(kappa, alpha, omega, X, CC, logP, integral,AreaMRW, WMR);
    CC = T[[1]];
    logP = T[[2]];
    integral = T[[3]];

    if( (step %% Jump) == 0){ pBX = rbind(pBX, S); if(Silent == FALSE){print(step); print(S); print(dim(CC)[1]);}};
  };
  pBX
};

MCMCestThomas = function(X, A1, A2,B1, B2, MR = 0.1, NStep = 100000, DiscardStep = 10000, Jump=10, PPalpha2, PPomega2, salpha, somega)
{
  pBX = pBetaX(X, NStep, Jump, A1, A2,B1, B2, AreaW, AreaMRW, WMR, PPalpha2, PPomega2, salpha, somega)
  Postkappa = pBX[(round(DiscardStep/Jump)+1):(dim(pBX)[1]),1]
  Postalpha = pBX[(round(DiscardStep/Jump)+1):(dim(pBX)[1]),2]
  Postomega = pBX[(round(DiscardStep/Jump)+1):(dim(pBX)[1]),3]
  PostlogP = pBX[(round(DiscardStep/Jump)+1):(dim(pBX)[1]),5]

  Kappahat = median(Postkappa)
  Alphahat = median(Postalpha)
  Omegahat = median(Postomega)
  KappaCI = quantile(Postkappa, probs = c(0.025,0.975))
  AlphaCI = quantile(Postalpha, probs = c(0.025,0.975))
  OmegaCI = quantile(Postomega, probs = c(0.025,0.975))

  res <- list(Kappahat=Kappahat,
              Alphahat=Alphahat,
              Omegahat=Omegahat,
              KappaCI=KappaCI,
              AlphaCI=AlphaCI,
              OmegaCI=OmegaCI,
              pBX=pBX[(round(DiscardStep/Jump)+1):(dim(pBX)[1]),]
  )

  res
}








################################################################################################

Output = MCMCestThomas(X, A1, A2,B1, B2, MR = MR, NStep = NStep, DiscardStep = DiscardStep, Jump = Jump, PPalpha2 = PPalpha2, PPomega2 = PPomega2, salpha = salpha, somega = somega)

###################                 Results
###########################################################################
#Estimate for kappa (Intensity of parent process)
Output$Kappahat
#Estimate for alpha (Mean number of points in a cluster)
Output$Alphahat
#Estimate for omega (standard deviation of cluster size)
Output$Omegahat

#Aposterior confidence interval for kappa (Intensity of parent process)
Output$KappaCI
#Aposterior confidence interval for alpha (Mean number of points in a cluster)
Output$AlphaCI
#Aposterior confidence interval for omega (standard deviation of cluster size)
Output$OmegaCI


#Aposterior distributions
Output$Postkappa = Output$pBX[,1]
Output$Postalpha = Output$pBX[,2]
Output$Postomega = Output$pBX[,3]
Output$PostlogP  = Output$pBX[,5]

Output$Postkappa
Output$Postalpha
Output$Postomega
Output$PostlogP



#Aposterior histogram for kappa (Intensity of parent process)
hist(Output$Postkappa)
#Aposterior histogram for alpha (Mean number of points in a cluster)
hist(Output$Postalpha)
#Aposterior histogram for omega (standard deviation of cluster size)
hist(Output$Postomega)

#Aposterior MCMC trace for kappa
plot(Output$Postkappa,type = 'l',xaxt = "n", main = "Aposterior MCMC trace for kappa")
axis(1, at = (0:4) * (length(Output$Postkappa) / 4), labels = c(DiscardStep,DiscardStep+(NStep-DiscardStep)/4,DiscardStep+(NStep-DiscardStep)/4*2,DiscardStep+(NStep-DiscardStep)/4*3,DiscardStep+(NStep-DiscardStep)))
lines(rep(Output$Kappahat, length(Output$Postkappa)))
lines(rep(Output$KappaCI[1], length(Output$Postkappa)), lty=2)
lines(rep(Output$KappaCI[2], length(Output$Postkappa)), lty=2)

#Aposterior MCMC trace for alpha
plot(Output$Postalpha,type = 'l', xaxt = "n", main = "Aposterior MCMC trace for alpha")
axis(1, at = (0:4) * (length(Output$Postkappa) / 4), labels = c(DiscardStep, DiscardStep + (NStep - DiscardStep) / 4, DiscardStep + (NStep - DiscardStep) / 4 * 2, DiscardStep + (NStep - DiscardStep) / 4 * 3, DiscardStep + (NStep - DiscardStep)))
lines(rep(Output$Alphahat, length(Output$Postkappa)))
lines(rep(Output$AlphaCI[1], length(Output$Postkappa)), lty=2)
lines(rep(Output$AlphaCI[2], length(Output$Postkappa)), lty=2)

#Aposterior MCMC trace for omega
plot(Output$Postomega,type = 'l',xaxt = "n", main="Aposterior MCMC trace for omega")
axis(1, at=(0:4)*(length(Output$Postkappa)/4), labels=c(DiscardStep,DiscardStep+(NStep-DiscardStep)/4,DiscardStep+(NStep-DiscardStep)/4*2,DiscardStep+(NStep-DiscardStep)/4*3,DiscardStep+(NStep-DiscardStep)))
lines(rep(Output$Omegahat,length(Output$Postkappa)))
lines(rep(Output$OmegaCI[1], length(Output$Postkappa)), lty=2)
lines(rep(Output$OmegaCI[2], length(Output$Postkappa)), lty=2)

#Aposterior MCMC trace for log likelihood
plot(Output$PostlogP,type = 'l', xaxt = "n", main = "Aposterior MCMC trace for log likelihood")
axis(1, at = (0:4) * (length(Output$Postkappa) / 4), labels=c(DiscardStep,DiscardStep+(NStep-DiscardStep)/4,DiscardStep+(NStep-DiscardStep)/4*2,DiscardStep+(NStep-DiscardStep)/4*3,DiscardStep+(NStep-DiscardStep)))

Output
}
