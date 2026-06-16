kappa <- 10
omega <- .1
lambda <- .5
theta <- 10


set.seed(1234)

X <- rgtp(kappa, omega, lambda, theta, win = owin(c(0, 1), c(0, 1)))

testthat::test_that("Simulation output structure", {
  testthat::expect_type(X, "list")
})

testthat::test_that("Simulation output values", {
  testthat::expect_equal(X$X$n, 186)
  testthat::expect_equal(X$C$n, 25)
  testthat::expect_equal(X$X$x[10], 0.7692580174962445)
  testthat::expect_equal(X$C$y[21], 0.03907188922166815)
})


a_kappa <- 4
b_kappa <- 1

a_omega <- -3
b_omega <- 1

l_lambda <- -1
u_lambda <- 0.99

a_theta <- 4
b_theta <- 1


set.seed(4321)

pdf(file = NULL)

est <- estgtp(X$X,
  skappa = exp(a_kappa + ((b_kappa^2) / 2)) / 100,
  somega = exp(a_omega + ((b_omega^2) / 2)) / 100,
  dlambda = 0.01,
  stheta = exp(a_theta + ((b_theta^2) / 2)) / 100, smove = 0.1,
  a_kappa = a_kappa, b_kappa = b_kappa,
  a_omega = a_omega, b_omega = b_omega,
  l_lambda = l_lambda, u_lambda = u_lambda,
  a_theta = a_theta, b_theta = b_theta,
  iter = 50, plot.step = 50, save.step = 1e9,
  filename = ""
)

dev.off()


testthat::test_that("Estimation output structure", {
  testthat::expect_type(est, "list")
})

testthat::test_that("Simulation output values", {
  testthat::expect_equal(est$info[[2]], 480)
  testthat::expect_equal(est$n_C[50], 36)
  testthat::expect_equal(est$theta[10], 36.22727692626274)
  testthat::expect_equal(est$lambda[21], 0.006447300063446165)
})
