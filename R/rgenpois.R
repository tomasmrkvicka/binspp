rgencdf <- function( lambda, theta, tol = 1e-6 ){

#  require( VGAM )
  m <- Inf
  if( lambda < 0 ){
    m <- floor( theta/ -lambda )
  }


  if( m<2 ){

    stop( 'm is too small')

  }


  if( m>100 ){
    m <- ceiling( theta/(1-lambda) + 5*theta/(1-lambda)^3 )
#    while( dgenpois( m, lambda, theta)>tol){
    while (dgenpois0( m, lambda = lambda, theta = theta) > tol){
      m <- ceiling( m + theta / (1 - lambda)^3 )
    }
  }

#  cdf_tmp <- cumsum( dgenpois( x = 0:m,lambda, theta ))
  cdf_tmp <- cumsum( dgenpois0( x = 0:m,lambda = lambda, theta = theta ))
  cdf_tmp <- cdf_tmp/cdf_tmp[m+1]

}

rgenpois <- function( n, lambda, theta, tol = 1e-6, cdf_tmp = NULL ){
#  require( VGAM )
  m <- Inf
  if( lambda < 0 ){
    m <- floor( theta/ -lambda )
  }


  if( m<2 ){

    stop( 'm is too small')

  }

  if( is.null( cdf_tmp) ){
    if( m>100 ){
      m <- ceiling( theta/(1-lambda) + 5*theta/(1-lambda)^3 )
#      while( dgenpois( m, lambda, theta)>tol){
      while( dgenpois0( m, lambda = lambda, theta = theta) > tol){
        m <- ceiling( m+ theta/(1-lambda)^3 )
      }
    }

#    cdf_tmp <- cumsum( dgenpois( x = 0:m,lambda, theta ))
    cdf_tmp <- cumsum( dgenpois0( x = 0:m,lambda = lambda, theta = theta ))
    cdf_tmp <- cdf_tmp/cdf_tmp[m+1]
  }
  x <- runif( n )

  sapply( x, function( k ) min( which(cdf_tmp>k))-1 )


}
