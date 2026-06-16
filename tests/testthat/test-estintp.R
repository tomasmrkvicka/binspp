# Prepare the dataset:
X <- trees_N4
x_left <- x_left_N4
x_right <- x_right_N4
y_bottom <- y_bottom_N4
y_top <- y_top_N4

z_beta <- list(refor = cov_refor, slope = cov_slope)
z_alpha <- list(tmi = cov_tmi, tdensity = cov_tdensity)
z_omega <- list(slope = cov_slope, reserv = cov_reserv)

# Determine the union of rectangles:
W <- owin(c(x_left[1], x_right[1]), c(y_bottom[1], y_top[1]))
if (length(x_left) >= 2) {
  for (i in 2:length(x_left)) {
    W2 <- owin(c(x_left[i], x_right[i]), c(y_bottom[i], y_top[i]))
    W <- union.owin(W, W2)
  }
}

# Dilated observation window:
W_dil <- dilation.owin(W, 100)


# User-specified hyperparameters for prior distributions:
control <- list(
  NStep = 100, BurnIn = 20, SamplingFreq = 5,
  Prior_alpha_mean = 3, Prior_alpha_SD = 2, Prior_omega_mean = 5.5,
  Prior_omega_SD = 5, Prior_alphavec_SD = c(4.25, 0.012),
  Prior_omegavec_SD = c(0.18, 0.009)
)


set.seed(1234)


Output <- estintp(
  X = X, control = control, x_left = x_left, x_right = x_right,
  y_bottom = y_bottom, y_top = y_top, W_dil = W_dil, z_beta = z_beta,
  z_alpha = z_alpha, z_omega = z_omega, verbose = FALSE
)


testthat::test_that("Output structure", {
  testthat::expect_s3_class(Output, "output_estintp")
  testthat::expect_type(Output, "list")
})

testthat::test_that("Output values", {
  testthat::expect_equal(Output$betahat.firststep[[1]], -0.001053888899525035)
  testthat::expect_equal(Output$betahat.firststep[[2]], -0.037289748386392770)

  testthat::expect_equal(Output$Pvalbetahat[[1]], 1.068480277588435e-07)
  testthat::expect_equal(Output$Pvalbetahat[[2]], 3.842628168775053e-02)

  testthat::expect_equal(Output$Alphahat, 3.262512207445075)
  testthat::expect_equal(Output$AlphaCI[[1]], 3.142701135709731)
  testthat::expect_equal(Output$AlphaCI[[2]], 3.461583874829456)

  testthat::expect_equal(Output$Alphabethat[1], -1.5229702132813773829)
  testthat::expect_equal(Output$Alphabethat[2], -0.0003667534484401336)

  testthat::expect_equal(Output$Omegahat, 7.520976616511347)

  testthat::expect_equal(Output$noOmega, 2)
})


raw <- rawMCMCoutput(Output)


testthat::test_that("Raw output structure", {
  testthat::expect_type(raw, "list")
})

testthat::test_that("Raw output values", {
  testthat::expect_equal(raw$raw.logP[21], 20364259.04628616)
  testthat::expect_equal(raw$raw.kappa[21], 5.448136154312056e-06)
  testthat::expect_equal(raw$raw.noParents[21], 75)
  testthat::expect_equal(raw$raw.acceptsAlpha[21], 0.050999999999999997)
})


set.seed(4321)

sim <- simulate(Output)


testthat::test_that("Simulation from the fitted model, output structure", {
  testthat::expect_s3_class(sim, "ppp")
  testthat::expect_type(sim, "list")
})

testthat::test_that("Simulation from the fitted model, output values", {
  testthat::expect_equal(sim$n, 370)
  testthat::expect_equal(sim$x[1], 6843.494772745912)
  testthat::expect_equal(sim$y[7], 3040.868843744716)
})
