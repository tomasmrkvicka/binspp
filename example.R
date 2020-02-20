require(binspp)

### True parameters
# lambda= -.5
# theta = 10
# kappa = 50
# omega = .01

### Hyperpriors
# a_kappa = 1.5; b_kappa = .02
# a_omega = 2; b_omega = 50
# a_theta = 2.5; b_theta = .25
# l_lambda <- -1; u_lambda <- .99


### Set working directory
# setwd("C:/path/to/directory")


X = rgtp(kappa= 10, omega = 0.1, lambda = 0.5, theta = 10, win = owin(c(0,1),c(0,1)))
plot(X$X)
plot(X$C)


a_kappa = 4
b_kappa = 1
x <- seq(0, 100, length=100)
hx <- dlnorm(x, a_kappa, b_kappa)
plot(x, hx, type="l", lty=1, xlab="x value",
     ylab="Density", main="Prior")

a_omega = -3
b_omega = 1
x <- seq(0, 1, length=100)

hx <- dlnorm(x,a_omega,b_omega)
plot(x, hx, type="l", lty=1, xlab="x value",
     ylab="Density", main="Prior")

l_lambda = -1
u_lambda = 0.99
x <- seq(-1, 1, length=100)
hx <- dunif(x,l_lambda,u_lambda)
plot(x, hx, type="l", lty=1, xlab="x value",
     ylab="Density", main="Prior")

a_theta = 4
b_theta = 1
x <- seq(0, 100, length=100)
hx <- dlnorm(x,a_theta,b_theta)
plot(x, hx, type="l", lty=1, xlab="x value",
     ylab="Density", main="Prior")


est = estgtp(X$X, skappa = exp(a_kappa+((b_kappa^2)/2))/100, somega = exp(a_omega+((b_omega^2)/2))/100,
              dlambda = 0.01,
              stheta = exp(a_theta+((b_theta^2)/2))/100, smove = 0.1,
              a_kappa = a_kappa, b_kappa = b_kappa, 
              a_omega = a_omega, b_omega = b_omega,
              l_lambda = l_lambda, u_lambda = u_lambda,
              a_theta = a_theta, b_theta = b_theta,
              iter = 1000, plot.step = 1000, save.step = 1e9,
              filename = "")


### Discard - number of iterations to be discarded as burn in for the estimation
discard = 100;

### Step - every step iteration is taken in the parameter estimation.  
step = 10; 

result = estgtpr(est, discard, step);

result
