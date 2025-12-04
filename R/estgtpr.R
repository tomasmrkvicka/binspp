#' Results for Bayesian MCMC estimation of parameters of generalized Thomas process
#'
#' @description Calculates median values for \emph{kappa}, \emph{omega}, \emph{lambda}, \emph{theta}; calculates 2.5 and 97.5 quantile and draws trace plots.
#'
#' @param est Output from [estgtp()] function.
#' @param discard Number of iterations to be discarded as burn in for the estimation.
#' @param step Every \emph{step} iteration is taken in the parameter estimation.
#'
#' @return Median and quantile values and plots (\emph{kappa}, \emph{omega}, \emph{lambda}, \emph{theta}).
#'
#' @md
#' @examples
#'
#' library(spatstat)
#' kappa = 10
#' omega = .1
#' lambda= .5
#' theta = 10
#'
#' X = rgtp(kappa, omega, lambda, theta, win = owin(c(0, 1), c(0, 1)))
#' plot(X$X)
#' plot(X$C)
#'
#' a_kappa = 4
#' b_kappa = 1
#' x <- seq(0, 100, length = 100)
#' hx <- dlnorm(x, a_kappa, b_kappa)
#' plot(x, hx, type = "l", lty = 1, xlab = "x value",
#'      ylab = "Density", main = "Prior")
#'
#' a_omega = -3
#' b_omega = 1
#' x <- seq(0, 1, length = 100)
#' hx <- dlnorm(x, a_omega, b_omega)
#' plot(x, hx, type = "l", lty = 1, xlab = "x value",
#'      ylab = "Density", main = "Prior")
#'
#' l_lambda = -1
#' u_lambda = 0.99
#' x <- seq(-1, 1, length = 100)
#'
#' hx <- dunif(x, l_lambda, u_lambda)
#' plot(x, hx, type = "l", lty = 1, xlab = "x value",
#'      ylab = "Density", main = "Prior")
#'
#' a_theta = 4
#' b_theta = 1
#' x <- seq(0, 100, length = 100)
#' hx <- dlnorm(x, a_theta, b_theta)
#' plot(x, hx, type = "l", lty = 1, xlab = "x value",
#'      ylab = "Density", main = "Prior")
#'
#' est = estgtp(X$X,
#'           skappa = exp(a_kappa + ((b_kappa ^ 2) / 2)) / 100,
#'           somega = exp(a_omega + ((b_omega ^ 2) / 2)) / 100,
#'           dlambda = 0.01,
#'           stheta = exp(a_theta + ((b_theta ^ 2) / 2)) / 100, smove = 0.1,
#'           a_kappa = a_kappa, b_kappa = b_kappa,
#'           a_omega = a_omega, b_omega = b_omega,
#'           l_lambda = l_lambda, u_lambda = u_lambda,
#'           a_theta = a_theta, b_theta = b_theta,
#'           iter = 50, plot.step = 50, save.step = 1e9,
#'           filename = "")
#'
#' discard = 10
#' step = 10
#'
#' result = estgtpr(est, discard, step)
#'
#' @export
#'
estgtpr <- function(est, discard = 100, step = 10) {

        vyst=NULL
        v = seq(discard, length(est$omega), step)

        #median kappa
        vyst$median$kappa = median(est$kappa[v])
        #95% confidence interval for kappa
        vyst$quantiles$kappa = quantile(est$kappa[v], c(0.025, 0.975))

        plot(est$kappa[v], type = 'l', xaxt = "n")
        axis(1, at = (0:4) * (length(v) / 4),
             labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4,
                        v[1] + (v[length(v)] - v[1]) / 4 * 2,
                        v[1] + (v[length(v)] - v[1]) / 4 * 3,
                        v[length(v)]))

        lines(rep(median(est$kappa[v]), length(v)))
        lines(rep(quantile(est$kappa[v], c(0.025)), length(v)), lty=2)
        lines(rep(quantile(est$kappa[v], c(0.975)), length(v)), lty=2)

        #median omega
        vyst$median$omega=median(est$omega[v])
        #95% confidence interval for omega
        vyst$quantiles$omega=quantile(est$omega[v], c(0.025, 0.975))

        plot(est$omega[v], type = 'l', xaxt = "n")
        axis(1, at = (0:4) * (length(v) / 4),
             labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4,
                        v[1] + (v[length(v)] - v[1]) / 4 * 2,
                        v[1] + (v[length(v)] - v[1]) / 4 * 3,
                        v[length(v)]))

        lines(rep(median(est$omega[v]), length(v)))
        lines(rep(quantile(est$omega[v], c(0.025)), length(v)), lty=2)
        lines(rep(quantile(est$omega[v], c(0.975)), length(v)), lty=2)

        #median lambda
        vyst$median$lambda=median(est$lambda[v])
        #95% confidence interval for lambda
        vyst$quantiles$lambda=quantile(est$lambda[v], c(0.025, 0.975))

        plot(est$lambda[v], type = 'l', xaxt = "n")
        axis(1, at=(0:4) * (length(v) / 4),
             labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4,
                        v[1] + (v[length(v)] - v[1]) / 4 * 2,
                        v[1] + (v[length(v)] - v[1]) / 4 * 3,
                        v[length(v)]))
        lines(rep(median(est$lambda[v]), length(v)))
        lines(rep(quantile(est$lambda[v], c(0.025)), length(v)), lty=2)
        lines(rep(quantile(est$lambda[v], c(0.975)), length(v)), lty=2)


        #median theta
        vyst$median$theta=median(est$theta[v])
        #95% confidence interval for theta
        vyst$quantiles$theta=quantile(est$theta[v],c(0.025,0.975))


        plot(est$theta[v],type = 'l',xaxt = "n")
        axis(1, at=(0:4) * (length(v) / 4),
             labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4,
                        v[1] + (v[length(v)] - v[1]) / 4 * 2,
                        v[1] + (v[length(v)] - v[1]) / 4 * 3,
                        v[length(v)]))
        lines(rep(median(est$theta[v]), length(v)))
        lines(rep(quantile(est$theta[v], c(0.025)), length(v)), lty=2)
        lines(rep(quantile(est$theta[v], c(0.975)), length(v)), lty=2)

        return (vyst)
}
