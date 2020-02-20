# rgtp <- function( kappa = NULL, omega, lambda, theta, win, nsim = 1, expand = 4*omega, C = NULL ){
#
#   if( !is.null( C )){
#     nsim <- length( C )
#   }
#
#   if( is.null( C )){
#     Cwin <- boundingbox(expand.owin(win, distance = expand ))
#     C <- rpoispp( kappa, win = Cwin, nsim = nsim, drop = F )
#   }
#
#   C <- C[sapply( C, getElement, 'n')>0]
#   Cwin <- C[[1]]$win
#   N_daught <- lapply( C, function(x) rgenpois( x$n, lambda, theta, cdf_tmp = cdf_tmp) )
#   daught <- lapply( 1:length( C ), function( k ) do.call( rbind, lapply( 1:length( N_daught[[k]] ), function(z) cbind( rnorm(N_daught[[k]][z],C[[k]]$x[z],omega),rnorm(N_daught[[k]][z],C[[k]]$y[z],omega ) ) ) ))
#   P <- list()
#   for( i in 1:length( C )){
#     P[[i]] <- rep( 1:C[[i]]$n, N_daught[[i]])
#     P[[i]] <- P[[i]][which( (daught[[i]][,1]>=win$xrange[1])&(daught[[i]][,1]<=win$xrange[2])&
#                               (daught[[i]][,2]>=win$yrange[1])&(daught[[i]][,2]<=win$yrange[2]) )]
#   }
#
#   daught <- lapply( daught, as.ppp, W = win)
#   daught <- lapply( daught, function(x) {attr( x,'rejects') <- NULL; x} )
#
#   daught <- lapply( 1:length( daught), function(k) list( x = daught[[k]]$x, y = daught[[k]]$y, P = P[[k]], xrange = win$xrange, yrange = win$yrange))
#   parent <- lapply( 1:length( C), function(k) list( x = C[[k]]$x, y = C[[k]]$y, D = N_daught[[k]], xrange = Cwin$xrange, yrange = Cwin$yrange))
#
#   if( nsim == 1){
#     daught <- daught[[1]]
#     parent <- parent[[1]]
#   }
#   list( X = daught, C = parent)
#
# }


#' RGTP
#' 
#' @description Simulation of generalized Thomas process.
#' @param kappa: Intensity of the Poisson process of cluster centres. A single positive number.
#' @param omega: Standard deviation of random displacement (along each coordinate axis) of a point from its cluster centre.
#' @param lambda: Shape parameter of Generalized Poisson distribution for number of points in a cluster.
#' @param theta: Location parameter of Generalized Poisson distribution for number of points in a cluster.
#' @param win: Window in which to simulate the pattern. An object of class "owin" or something acceptable to as.owin.
#' @param nsim: Number of simulations.
#' @param expand: The size of expansion of window to simulate the centers of clusters. 
#' @param C: Process of center points.
#'
#' @return A list(X, C), where X is Generalized Poisson process, and C is Process of cluster centers for Generalized Poisson process.
#' @export
#'
#' @examples
rgtp <- function( kappa, omega, lambda, theta, win = owin(c(0,1),c(0,1)), nsim = 1, expand = 4*omega, C = NULL ){

  marks = FALSE; cdf_tmp = NULL
  if( !is.null( C )){
    nsim <- length( C )
  }

  if( is.null( C )){
    winar <- (diff( win$xrange) + 2*expand)*(diff( win$yrange ) +2*expand)
    # Cwin <- boundingbox(expand.owin(win, distance = expand ))
    Nc <- rpois( nsim, kappa*winar )
    # C <- rpoispp( kappa, win = Cwin, nsim = nsim, drop = F )
    C <- lapply( 1:nsim, function(x) cbind( runif( Nc[x], win$xrange[1]-expand, win$xrange[2] + expand),
                runif( Nc[x], win$yrange[1]-expand, win$yrange[2] + expand) ) )
  }

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
  list( X = daught, C = parent)

}

inwin <- function( X, win ){

 X <- matrix( X[ (X[,1]>=win$xrange[1])&((X[,1]<=win$xrange[2]))&
           (X[,2]>=win$yrange[1])&((X[,2]<=win$yrange[2])), ], ncol = 2)

 X

}

