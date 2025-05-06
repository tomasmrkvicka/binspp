#' Simulation of generalized Thomas process
#'
#' @description Simulation of generalized Thomas process.
#'
#' @param kappa intensity of cluster centers.
#' @param omega standard deviation of normal distribution specifying the clusters spread.
#' @param lambda parameter of generalised Poisson distribution controlling over or under dispersion.
#' @param theta parameter of generalised Poisson distribution controlling the mean number of points in a cluster.
#' @param win window in which to simulate the pattern. An object in the \cr [spatstat.geom::owin()] format of the \pkg{spatstat} package.
#' @param nsim number of simulations.
#' @param expand the size of expansion of window to simulate the centers of clusters.
#'
#' @return A list(X, C), where \emph{X} is Generalized Thomas process, and \emph{C} is Process of cluster centers for Generalized Thomas process.
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
#' @export
#'
rgtp <- function( kappa, omega, lambda, theta, win = owin(c(0, 1), c(0, 1)), nsim = 1, expand = 4 * omega ){

  marks = FALSE; cdf_tmp = NULL

  winar <- (diff( win$xrange) + 2 * expand) * (diff( win$yrange ) + 2 * expand)
  # Cwin <- boundingbox(expand.owin(win, distance = expand ))
  Nc <- rpois( nsim, kappa * winar )
  # C <- rpoispp( kappa, win = Cwin, nsim = nsim, drop = F )
  C <- lapply( 1:nsim, function(x) cbind( runif( Nc[x], win$xrange[1]-expand, win$xrange[2] + expand),
                runif( Nc[x], win$yrange[1]-expand, win$yrange[2] + expand) ) )

  C <- C[Nc>0]
  Nc <- Nc[Nc>0]
  Cwin <- list( xrange = win$xrange + c(-1,1)*expand, yrange = win$yrange + c(-1,1)*expand)
  N_daught <- lapply( Nc, function(x) rgenpois( x, lambda, theta, cdf_tmp = cdf_tmp) )
  if( any( sapply( N_daught, length) == 0 ) ){browser()}
  daught <- lapply( 1:length( C ), function( k ) do.call( rbind, lapply( 1:length( N_daught[[k]] ),
                                                                         function(z) cbind( rnorm(N_daught[[k]][z],C[[k]][z,1],omega),
                                                                                            rnorm(N_daught[[k]][z],C[[k]][z,2],omega ) ) ) ))
  P <- list()
  for( i in 1:length( C )){
    P[[i]] <- rep( 1:Nc[i], N_daught[[i]])
    P[[i]] <- P[[i]][which( (daught[[i]][,1]>=win$xrange[1])&(daught[[i]][,1]<=win$xrange[2])&
                              (daught[[i]][,2]>=win$yrange[1])&(daught[[i]][,2]<=win$yrange[2]) )]
  }

  # daught <- lapply( daught, as.ppp, W = win)
  # daught <- lapply( daught, function(x) {attr( x,'rejects') <- NULL; x} )
  daught <- lapply( daught, inwin, win = win )

  if( any( sapply( daught, function(x) length(dim(x)))<2)){browser()}

  daught <- lapply( 1:length( daught), function(k) list( x = daught[[k]][,1], y = daught[[k]][,2], P = P[[k]], xrange = win$xrange, yrange = win$yrange))
  parent <- lapply( 1:length( C), function(k) list( x = C[[k]][,1], y = C[[k]][,2], D = N_daught[[k]], xrange = Cwin$xrange, yrange = Cwin$yrange))

  if( nsim == 1){
    daught <- daught[[1]]
    parent <- parent[[1]]
  }
  res=list( X = daught, C = parent)
  return(res)
}

inwin <- function( X, win ){
 X <- matrix( X[ (X[,1]>=win$xrange[1])&((X[,1]<=win$xrange[2]))&
           (X[,2]>=win$yrange[1])&((X[,2]<=win$yrange[2])), ], ncol = 2)
 return(X)
}

