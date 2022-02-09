# 1
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov.refor, slope=cov.slope)
z_alpha = list(tmi=cov.tmi, tdensity=cov.tdensity)
z_omega = list(slope=cov.slope, reserv=cov.reserv)

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


# User-specified hyperparameters for prior distributions:
# control = list(NStep=500, BurnIn=100, SamplingFreq=10, Prior_alpha_mean=3, Prior_alpha_SD=2,
#                Prior_omega_mean=5.5, Prior_omega_SD=5, Prior_alphavec_SD=c(4.25,0.012),
#                Prior_omegavec_SD=c(0.18,0.009))

# Default parameters for prior distributions:
control = list(NStep=500, BurnIn=100, SamplingFreq=10)


# MCMC estimation:
Output = binspp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Which hyperparameter values for prior distribution were used?
Output$priorParameters


# Text output + series of figures:
print.outputs(Output)
plot.outputs(Output)


# Recompute the outputs when another value of burn-in is desired,
# without running the chain again:
Out2 <- re_estimate(Output, BurnIn=200)
print.outputs(Out2)
plot.outputs(Out2)



# 2
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov.refor, slope=cov.slope)
z_alpha = list(tmi=cov.tmi, tdensity=cov.tdensity)
z_omega = list(slope=cov.slope, reserv=cov.reserv)

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
Output = binspp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Which hyperparameter values for prior distribution were used?
Output$priorParameters


# Text output + series of figures:
print.outputs(Output)
plot.outputs(Output)


# Recompute the outputs when another value of burn-in is desired,
# without running the chain again:
Out2 <- re_estimate(Output, BurnIn=200)
print.outputs(Out2)
plot.outputs(Out2)



#3
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov.refor, slope=cov.slope)

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
aux = first.step(X, z_beta, W_dil, plot=TRUE)



#4
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov.refor, slope=cov.slope)
z_alpha = list(tmi=cov.tmi, tdensity=cov.tdensity)
z_omega = list(slope=cov.slope, reserv=cov.reserv)

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
Output = binspp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Text output + series of figures:
print.outputs(Output)
plot.outputs(Output)



#5
# Prepare the dataset:
X = trees_N4
x_left = x_left_N4
x_right = x_right_N4
y_bottom = y_bottom_N4
y_top = y_top_N4

z_beta = list(refor=cov.refor, slope=cov.slope)
z_alpha = list(tmi=cov.tmi, tdensity=cov.tdensity)
z_omega = list(slope=cov.slope, reserv=cov.reserv)

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
Output = binspp(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Text output + series of figures:
print.outputs(Output)
plot.outputs(Output)


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

