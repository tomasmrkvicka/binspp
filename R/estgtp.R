plot.conn <- function( X, C ){

  plot( C, pch = 19, col = 2, asp = 1, xlim = C$xrange, ylim = C$yrange)
  points( X, pch = 19, cex = .5)
  for( i in 1:length(X$x)){

    lines( c(X$x[i],C$x[X$P[i]]), c(X$y[i],C$y[X$P[i]]))

  }

}

#' estgtp
#'
#' @description Bayesian MCMC estimation of parameters of generalized Thomas process.
#' @param X: A point pattern dataset (object of class "ppp") to which the model should be fitted.
#' @param kappa0: Initial value for kappa, by default it will be set as expectation of prior for kappa.
#' @param omega0 Initial value for omega, by default it will be set as expectation of prior for omega.
#' @param lambda0 Initial value for lambda, by default it will be set as expectation of prior for lambda.
#' @param theta0 Initial value for theta, by default it will be set as expectation of prior for theta.
#' @param skappa variability of proposal for kappa: second parameter of log-normal distribution
#' @param somega variability of proposal for omega: second parameter of log-normal distribution
#' @param dlambda variability of proposal for lambda: half of range of uniform distribution
#' @param stheta variability of proposal for theta: second parameter of log-normal distribution
#' @param smove variability of proposal for moving center point: SD of normal distribution 
#' @param a_kappa: First parameter of prior distribution for kappa, which is log-normal distribution. 
#' @param b_kappa: Second parameter of prior distribution for kappa, which is log-normal distribution.
#' @param a_omega: First parameter of prior distribution for omega, which is log-normal distribution.
#' @param b_omega: Second parameter of prior distribution for omega, which is log-normal distribution.
#' @param l_lambda: First parameter of prior distribution for lambda, which is uniform distribution.
#' @param u_lambda: Second parameter of prior distribution for lambda, which is uniform distribution.
#' @param a_theta: First parameter of prior distribution for theta, which is log-normal distribution.
#' @param b_theta: Second parameter of prior distribution for theta, which is log-normal distribution.
#' @param iter: Number of iterations of MCMC.
#' @param plot.step: Step for the graph plotting
#' @param save.step: Step for the parameters saving. The file must be specified or has to be set to larger than iter.
#' @param filename: The name of the output RDS file
#'
#' @return The output is an estimated MCMC chain of parameters, centers and connections.
#' @export
#'
#' @examples
estgtp <- function( X,
                           kappa0 = exp(a_kappa + ((b_kappa ^ 2) / 2)),
                           omega0 = exp(a_omega + ((b_omega ^ 2) / 2)),
                           lambda0 = (l_lambda + u_lambda) / 2,
                           theta0 = exp(a_theta + ((b_theta ^ 2) / 2)),
                           skappa,
                           somega,
                           dlambda,
                           stheta,
                           smove,
                           a_kappa,
                           b_kappa,
                           a_omega,
                           b_omega,
                           l_lambda,
                           u_lambda,
                           a_theta,
                           b_theta,
                           iter = 500000,
                           plot.step = 1000,
                           save.step = 1000,
                           filename
){

#  require( looptimer)
  require( spatstat)
  require( FNN )

  inwin <- function(x,y){

    sum( (x>=0)*(x<=1)*(y>=0)*(y<=1))


  }


  # Get initial guess for centers
  C.step = 100; n_upd = 100; n_bdm = 100;
  C0 = NULL
  initial.guess = 'clever'
  if( is.null( C0) ){
    if( initial.guess != 'naive'){
      prep <- C_prep( X, lambda0, theta0, expand = 5*omega0 )
    }
    if( initial.guess == 'naive'){
      prep <- C_prep_naive( X,kappa0, lambda0, theta0, expand = 5*omega0 )
    }

    # Get X and C from the prepatation (connectiosn are created etc)
    X <- prep$X
    C <- prep$C
  }
  else{ C = C0 }
  # Count the number of daughters for each parent and compute
  # the integral of the gaussian kernel over the window
  # for each parent
  C$n <- sapply( 1:length(C$x), function(k) sum( X$P == k))
  X$P <- X$P-1
  C$mus <- intfun_cols_cpp(cbind(C$x, C$y), omega0, X$xrange, X$yrange )



  if( (lambda0<0)&(max( C$n )) >= floor(theta0/-lambda0) ){

    warning( "lambda0 value incompatible with daughter counts, setting lambda0 to 0" )
    lambda0 = 0
  }

  # Initialize the parameter vectors
  kappas <- thetas <- lambdas <- omegas <- vector( length = iter )
  otmp <- omegas[1] <- omega0
  ttmp <- thetas[1] <- theta0
  ltmp <- lambdas[1] <- lambda0
  ktmp <- kappas[1] <- kappa0
  int.hat <- length(X$x)/(diff(X$xrange)*diff(X$yrange))

  # Compute the generalized Poisson density (to use in dgenpoisbin)
  dgp <- dgp_comp_cpp( lambda0, theta0 )

  # Variables that count the number of deaths, births and update counts

  move <- switch <- birth_count <- death_count <- 0

  # Lists to store the simulated parent patterns and connections
  Plist <- Clist <- list()
  Clist[[1]] <- C
  Plist[[1]] <- X$P

  Cwar <- diff( C$xrange )*diff( C$yrange )
  # Vector that contains the number of daughters in each step
  n_Ciw <- n_Cow <- n_C <- vector( length = iter )
  n_C[1] <- length( C$x )
  n_Ciw[1] <- inwin( C$x, C$y )
  n_Cow[1] <- (n_C[1]-n_Ciw[1])/(Cwar-1)

  Cwin_ar = diff(C$xrange)*diff(C$yrange)

  # Time the loop
  #t0 <- looptimer()
  if (plot.step > iter) {
    warning("The value of 'plot.step' is greater than 'iter' parameter, no plots will be visible.", noBreaks. = TRUE)
  }

  for( i in 2:iter ){

    # Set the current parameter values to the most recent accepted
    ltmp <- lambdas[i-1]
    ttmp <- thetas[i-1]
    otmp <- omegas[i-1]
    ktmp <- kappas[i-1]

    # Update omega
    omegas[i] <- otmp
    omegas[i] <- update_omega(otmp, somega, X, C, ltmp, ttmp, dgp, a_omega, b_omega )

    # Update lambda and theta jointly
    upd_lt <- update_lt(ltmp, ttmp, stheta, dlambda, X, C, otmp, dgp, l_lambda, u_lambda, a_theta, b_theta )
    lambdas[i] <- upd_lt[1]
    thetas[i] <- upd_lt[2]
    dgp <- dgp_comp_cpp( lambdas[i], thetas[i] )
    # lambdas[i] <- ltmp
    # thetas[i] <- ttmp


    # Update kappa
    kappas[i]<- update_kappa( ktmp, skappa, length( C$x),Cwin_ar, a_kappa, b_kappa )
    # kappas[i] <- ktmp

    if( max( X$P )>(length(C$x)-1) ) browser()

    # Update connections n_upd times
    if( n_upd >0 ){
      switch = update_P_cpp_n(X, C, otmp, ltmp, ttmp, C$mus, dgp, n_upd, switch)
    }

    # Do birth/death/move n_bdm times
    # bdm( X, C, ktmp, otmp, ltmp, ttmp, smove, sbirth, n_bdm )
    sbirth = smove;
    #sbirth= 4*omega;
    if( n_bdm > 0 ){
      R = runif( n_bdm)
      for( j in 1:n_bdm ){
        if( R[j]<= .5) move = move_C_cpp(X, C, otmp, ltmp, ttmp, smove, dgp, move)
        if( (R[j]>.5)&(R[j]<=.75))  birth_count <- c_birth_cpp(X, C, ktmp, otmp, ltmp, ttmp, sbirth, dgp, birth_count )
        if( R[j]>.75 ) death_count = c_death_cpp(X, C, ktmp, otmp, ltmp, ttmp, sbirth, dgp, death_count )

      }
    }

    # Store parent count
    n_C[i] <- length( C$x )
    n_Ciw[i] <- inwin( C$x, C$y )
    n_Cow[i] <- (n_C[i]-n_Ciw[i])/(Cwar-1)

    # Plot, if the time is right
    if( i%%plot.step == 0 )  {
      Xtmp = X
      Xtmp$P = X$P+1
      par( mfrow = c(3,3))
      plot.conn(Xtmp,C)
      plot( kappas[1:(i-1)], type = 'l', main = expression(kappa))
      abline( h = kappa0, col = 2, lty = 2)
      plot( omegas[1:(i-1)], type = 'l', main = expression(omega))
      abline( h = omega0, lty = 2, col = 2)
      plot( thetas[1:(i-1)], type = 'l', main = expression(theta))
      abline( h = theta0, lty = 2, col = 2)
      plot( lambdas[1:(i-1)], type = 'l', main = expression(lambda))
      abline( h = lambda0, lty = 2, col = 2)
      plot( thetas[1:(i-1)]/(1-lambdas[1:(i-1)])*kappas[1:(i-1)], type = 'l', main = 'Intensity' )
      abline( h = int.hat, col = 2, lty = 2)
      plot( thetas[1:(i-1)]/(1-lambdas[1:(i-1)]), type = 'l', main = expression(mu) )
      abline( h = theta0/(1-lambda0), col = 2, lty = 2)
      # plot( n_C[1:(i-1)], type = 'l')
      plot( n_Ciw[1:(i-1)], type = 'l', ylim = range( c(n_Ciw[1:i], n_Cow[1:i]) ) )
      lines( n_Cow[1:(i-1)], col = 4)
      abline( h = n_Ciw[1], col = 2)
      abline( h = kappas[1], col = 2, lty = 2)

      unique( n1 <- sapply( 1:length( C$x ), function( k ) sum( X$P == k )) )
      mus.end <- intfun_cols_cpp(cbind(C$x, C$y ), omega, c(0,1), c(0,1))
      in_win <- which( mus.end >=.95 )
      barplot( (table( c(n1[in_win], 0:8) )-1)/sum( table( c(n1[in_win], 0:8) )-1), ylim = c(0,.5))
      barplot( dgenpois( 0:8, lambda, theta), col = rgb( 0,1,0,.1), add = T)

    }

    # Store parent pattern and connections
    if( i%%C.step == 0 ){Ctmp = C; Ctmp$tmp = 0;
    Clist[[i/C.step+1]]<-Ctmp;
    Plist[[i/C.step + 1]] <- X$P}

    # Save progress
    if( i%%save.step==0) saveRDS(object = list( C = C, X= X, C0 = C0, pars = list( kappas, omegas, lambdas, thetas )), filename)

    # Print the timing
 #   print( t0 <- looptimer(t0,iter,i,T) )
 #   cat("i= "); cat (i); cat("  iter = "); print(iter);

  }


  list( kappa = kappas, omega = omegas, lambda = lambdas, theta = thetas, Clist = Clist, Plist = Plist, n_C = n_C,
        n_Ciw = n_Ciw, n_Cow = n_Cow,
        info = c(birth = birth_count, death = death_count, move = move, updates = switch))

}
