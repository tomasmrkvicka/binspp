C_prep <- function( X, lambda, theta, expand = .05 ){
#  require( cluster )
#  require( FNN )
  m <- Inf
  if( lambda < 0) m <- theta/-lambda 
  
  Xxy <- cbind( X$x, X$y)
  
  mu <- theta/(1-lambda) 
  n_C <- ceiling(length(X$x)/mu)
  C_sugg <- pam( Xxy, k = n_C)
  
  P <- C_sugg$clustering
  D <- table( P )
  
  cents <- t( sapply( 1:n_C, function(k) colMeans( matrix( Xxy[P==k,], ncol = 2 ))) )
  
  
  C <- list( x = cents[,1], y = cents[,2], D = D, xrange = X$xrange + c(-expand,expand), 
             yrange = X$yrange + c(-expand,expand))
  X <- list( x = X$x, y = X$y, P = P, xrange = X$xrange, yrange= X$yrange)
  list( X = X, C = C)
}