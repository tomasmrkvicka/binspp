#' estgtpr
#' 
#' @description Calculates median values for kappa, omega, lambda, theta; calculates 2.5 and 97.5 quantile and draws trace plots
#'
#' @param est: Output from estgtp function. 
#' @param discard: Number of iterations to be discarded as burn in for the estimation.
#' @param step: Every step iteration is taken in the parameter estimation. 
#'
#' @return median and quantile values and plots (kappa, omega, lambda, theta) 
#' @export
#'
#' @examples
estgtpr <- function(est, discard = 100, step = 10) {

        v = seq(discard, length(est$omega), step)

        #median kappa
        median(est$kappa[v])
        
        #95% confidence interval for kappa 
        quantile(est$kappa[v], c(0.025, 0.975))
        

        plot(est$kappa[v], type = 'l', xaxt = "n")
        axis(1, at = (0:4) * (length(v) / 4), 
                labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4, 
                           v[1] + (v[length(v)] - v[1]) / 4 * 2, 
                           v[1] + (v[length(v)] - v[1]) / 4 * 3, 
                           v[length(v)]))
        lines(rep(median(est$kappa[v]),  length(v)))
        lines(rep(quantile(est$kappa[v], c(0.025)), length(v)), lty=2)     
        lines(rep(quantile(est$kappa[v], c(0.975)), length(v)), lty=2)  
        
        #median omega
        median(est$omega[v])
        
        #95% confidence interval for omega 
        quantile(est$omega[v], c(0.025, 0.975))
        
        plot(est$omega[v], type = 'l',xaxt = "n")
        axis(1, at = (0:4) * (length(v) / 4), 
                labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4,
                           v[1] + (v[length(v)] - v[1]) / 4 * 2,
                           v[1] + (v[length(v)] - v[1]) / 4 * 3,
                           v[length(v)]))
        lines(rep(median(est$omega[v]), length(v)))
        lines(rep(quantile(est$omega[v], c(0.025)), length(v)), lty=2)     
        lines(rep(quantile(est$omega[v], c(0.975)), length(v)), lty=2)  
        
        
        #median lambda
        median(est$lambda[v])
        
        #95% confidence interval for lambda 
        quantile(est$lambda[v], c(0.025, 0.975))
        
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
        median(est$theta[v])
        
        #95% confidence interval for theta 
        quantile(est$theta[v], c(0.025, 0.975))
        
        
        plot(est$theta[v],type = 'l',xaxt = "n")
        axis(1, at=(0:4) * (length(v) / 4), 
             labels = c(v[1], v[1] + (v[length(v)] - v[1]) / 4,
                        v[1] + (v[length(v)] - v[1]) / 4 * 2,
                        v[1] + (v[length(v)] - v[1]) / 4 * 3,
                        v[length(v)]))
        lines(rep(median(est$theta[v]), length(v)))
        lines(rep(quantile(est$theta[v], c(0.025)), length(v)), lty=2)     
        lines(rep(quantile(est$theta[v], c(0.975)), length(v)), lty=2)  

        return (v)
}
