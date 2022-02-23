# 1
# Prepare the dataset:
## please, set correct path to your dataset file and load it manually
##data("dataset_N4.RData")
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4
z_beta = list(refor=cov_refor, slope=cov_slope)
z_alpha = list(tmi=cov_tmi, tdensity=cov_tdensity)
z_omega = list(slope=cov_slope, reserv=cov_reserv)

# Determine the union of rectangles:
W = owin(c(x_left[1],x_right[1]),c(y_bottom[1],y_top[1]))
if(length(x_left)>=2){
  for(i in 2:length(x_left)){
    W2 = owin(c(x_left[i],x_right[i]),c(y_bottom[i],y_top[i]))
    W=union.owin(W,W2)
  }
}

# Dilated observation window:
W_dil = spatstat.geom::dilation.owin(W,100)


# User-specified hyperparameters for prior distributions:
# control = list(NStep=500, BurnIn=100, SamplingFreq=10, Prior_alpha_mean=3, Prior_alpha_SD=2,
#                Prior_omega_mean=5.5, Prior_omega_SD=5, Prior_alphavec_SD=c(4.25,0.012),
#                Prior_omegavec_SD=c(0.18,0.009))

# Default parameters for prior distributions:
control = list(NStep=500, BurnIn=100, SamplingFreq=10)


# MCMC estimation:
Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Which hyperparameter values for prior distribution were used?
Output$priorParameters


# Text output + series of figures:
print_outputs(Output)
plot_outputs(Output)


# Recompute the outputs when another value of burn-in is desired,
# without running the chain again:
Out2 <- re_estimate(Output, BurnIn=200)
print_outputs(Out2)
plot_outputs(Out2)


# 2
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov_refor, slope=cov_slope)
z_alpha = list(tmi=cov_tmi, tdensity=cov_tdensity)
z_omega = list(slope=cov_slope, reserv=cov_reserv)

# Determine the union of rectangles:
W = owin(c(x_left[1],x_right[1]),c(y_bottom[1],y_top[1]))
if(length(x_left)>=2){
  for(i in 2:length(x_left)){
    W2 = owin(c(x_left[i],x_right[i]),c(y_bottom[i],y_top[i]))
    W=union.owin(W,W2)
  }
}

# Dilated observation window:
W_dil = dilation.owin(W,100)


# Default parameters for prior distributions:
control = list(NStep=500, BurnIn=100, SamplingFreq=10)


# MCMC estimation:
Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Which hyperparameter values for prior distribution were used?
Output$priorParameters


# Text output + series of figures:
print_outputs(Output)
plot_outputs(Output)


# Recompute the outputs when another value of burn-in is desired,
# without running the chain again:
Out2 <- re_estimate(Output, BurnIn=200)
print_outputs(Out2)
plot_outputs(Out2)



#3
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov_refor, slope=cov_slope)

# Determine the union of rectangles:
W = owin(c(x_left[1],x_right[1]),c(y_bottom[1],y_top[1]))
if(length(x_left)>=2){
  for(i in 2:length(x_left)){
    W2 = owin(c(x_left[i],x_right[i]),c(y_bottom[i],y_top[i]))
    W=union.owin(W,W2)
  }
}

# Dilated observation window:
W_dil = dilation.owin(W,100)


# Estimating the intensity function of the parent process:
aux = first_step(X, z_beta, W_dil, plot=TRUE)



#4
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov_refor, slope=cov_slope)
z_alpha = list(tmi=cov_tmi, tdensity=cov_tdensity)
z_omega = list(slope=cov_slope, reserv=cov_reserv)

# Determine the union of rectangles:
W = owin(c(x_left[1],x_right[1]),c(y_bottom[1],y_top[1]))
if(length(x_left)>=2){
  for(i in 2:length(x_left)){
    W2 = owin(c(x_left[i],x_right[i]),c(y_bottom[i],y_top[i]))
    W=union.owin(W,W2)
  }
}

# Dilated observation window:
W_dil = dilation.owin(W,100)


# Default parameters for prior distributions:
control = list(NStep=500, BurnIn=100, SamplingFreq=10)


# MCMC estimation:
Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Text output + series of figures:
print_outputs(Output)
plot_outputs(Output)


#5
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov_refor, slope=cov_slope)
z_alpha = list(tmi=cov_tmi, tdensity=cov_tdensity)
z_omega = list(slope=cov_slope, reserv=cov_reserv)

# Determine the union of rectangles:
W = owin(c(x_left[1],x_right[1]),c(y_bottom[1],y_top[1]))
if(length(x_left)>=2){
  for(i in 2:length(x_left)){
    W2 = owin(c(x_left[i],x_right[i]),c(y_bottom[i],y_top[i]))
    W=union.owin(W,W2)
  }
}

# Dilated observation window:
W_dil = dilation.owin(W,100)


# Default parameters for prior distributions:
control = list(NStep=500, BurnIn=100, SamplingFreq=10)


# MCMC estimation:
Output = estintp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Text output + series of figures:
print_outputs(Output)
plot_outputs(Output)


#6
# Unit square observation window:
W <- owin()

# Dilation of the observation window:
W_dil <- dilation(W,0.1)
W_dil <- as.mask(W_dil)

# Define covariates:
f1 <- function(x,y){x}
f2 <- function(x,y){y}
f3 <- function(x,y){1-(y-0.5)^2}
cov1 <- as.im(f1, W=W_dil)
cov2 <- as.im(f2, W=W_dil)
cov3 <- as.im(f3, W=W_dil)


# Stationary Thomas process:
X <- rThomasInhom(kappa=50, alpha=log(10), omega=log(0.01), W=W, W_dil=W_dil)
plot(X)


# Thomas-type cluster process with inhomogeneity in all model components:
X <- rThomasInhom(kappa=10, betavec=c(1), z_beta=list(cov1), alpha=log(10), alphavec=c(1), z_alpha=list(cov2), omega=log(0.01), omegavec=c(1), z_omega=list(cov3), W=W, W_dil=W_dil)
plot(X)

# 7
# estgtpr.R

kappa = 10
omega = .1
lambda= .5
theta = 10

X = rgtp(kappa, omega, lambda, theta, win = owin(c(0, 1), c(0, 1)))
plot(X$X)
plot(X$C)

a_kappa = 4
b_kappa = 1
x <- seq(0, 100, length = 100)
hx <- dlnorm(x, a_kappa, b_kappa)
plot(x, hx, type = "l", lty = 1, xlab = "x value",
     ylab = "Density", main = "Prior")

a_omega = -3
b_omega = 1
x <- seq(0, 1, length = 100)
hx <- dlnorm(x, a_omega, b_omega)
plot(x, hx, type = "l", lty = 1, xlab = "x value",
     ylab = "Density", main = "Prior")

l_lambda = -1
u_lambda = 0.99
x <- seq(-1, 1, length = 100)

hx <- dunif(x, l_lambda, u_lambda)
plot(x, hx, type = "l", lty = 1, xlab = "x value",
     ylab = "Density", main = "Prior")

a_theta = 4
b_theta = 1
x <- seq(0, 100, length = 100)
hx <- dlnorm(x, a_theta, b_theta)
plot(x, hx, type = "l", lty = 1, xlab = "x value",
     ylab = "Density", main = "Prior")

est = estgtp(X$X,
          skappa = exp(a_kappa + ((b_kappa ^ 2) / 2)) / 100, somega = exp(a_omega + ((b_omega ^ 2) / 2)) / 100,
          dlambda = 0.01,
          stheta = exp(a_theta + ((b_theta ^ 2) / 2)) / 100, smove = 0.1,
          a_kappa = a_kappa, b_kappa = b_kappa,
          a_omega = a_omega, b_omega = b_omega,
          l_lambda = l_lambda, u_lambda = u_lambda,
          a_theta = a_theta, b_theta = b_theta,
          iter = 1000, plot.step = 1000, save.step = 1e9,
          filename = "")

discard = 100
step = 10

result = estgtpr(est, discard, step)
result
