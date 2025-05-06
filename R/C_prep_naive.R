C_prep_naive <- function( X, kappa, lambda, theta, expand = .05 ){
  
  xr <- X$xrange
  yr <- X$yrange
  
  xrc <- X$xrange + c(-expand,expand)
  yrc <- X$yrange + c(-expand,expand)
  
  
  Xxy <- cbind( X$x, X$y)
  
  m <- 1e9; nmax <- Inf
  if( lambda<0) m <- theta/-lambda
  
  while( nmax>m){
    C <- rpoispp(kappa, win = owin( xrange = xrc, yrange = yrc ))
    n_C <- C$n
    nmax <- ceiling( length(X$x)/n_C)
  }
  
  D <- rep( nmax-1, n_C)
  D[1:(length(X$x)-sum(D))] <- D[1:(length(X$x)-sum(D))] + 1
  P <- rep(1:n_C, D)
  
  cents <- cbind( C$x, C$y )
  
  C <- list( x = cents[,1], y = cents[,2], D = D, xrange = X$xrange + c(-expand,expand), 
             yrange = X$yrange + c(-expand,expand), n = D)
  X <- list( x = X$x, y = X$y, P = P, xrange = X$xrange, yrange= X$yrange)
  list( X = X, C = C)
}