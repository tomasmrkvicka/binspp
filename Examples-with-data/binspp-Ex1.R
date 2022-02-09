# 1
# Prepare the dataset:
## data("datasetN4")
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
W_dil = spatstat.geom::dilation.owin(W,100)


# User-specified hyperparameters for prior distributions:
# control = list(NStep=500, BurnIn=100, SamplingFreq=10, Prior_alpha_mean=3, Prior_alpha_SD=2,
#                Prior_omega_mean=5.5, Prior_omega_SD=5, Prior_alphavec_SD=c(4.25,0.012),
#                Prior_omegavec_SD=c(0.18,0.009))

# Default parameters for prior distributions:
control = list(NStep=500, BurnIn=100, SamplingFreq=10)


# MCMC estimation:
Output = binsppf(X, control, x_left, x_right, y_bottom, y_top, W_dil, z_beta, z_alpha, z_omega)


# Which hyperparameter values for prior distribution were used?
Output$priorParameters


# Text output + series of figures:
printOutputs(Output)
plotOutputs(Output)


# Recompute the outputs when another value of burn-in is desired,
# without running the chain again:
Out2 <- re_estimate(Output, BurnIn=200)
printOutputs(Out2)
plotOutputs(Out2)
