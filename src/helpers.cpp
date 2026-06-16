#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "math.h"

using namespace Rcpp;
using namespace Eigen;
using Eigen::MatrixXd;

NumericVector arma2nvec( arma::vec x){
  return( NumericVector(x.begin(), x.end()));
}

IntegerVector arma2ivec( arma::vec x){
  return( IntegerVector(x.begin(), x.end()));
}

NumericVector arma2nvec( arma::uvec x){
  return( NumericVector(x.begin(), x.end()));
}

IntegerVector arma2ivec( arma::uvec x){
  return( IntegerVector(x.begin(), x.end()));
}



// Compute the generalized poisson likelihood for each of observations x with parameters lambda and theta.
// Pass log_ret = True to get log likelihoods
// [[Rcpp::export]]
NumericVector dgenpois_cpp( NumericVector x, double lambda, double theta, bool log_ret ){

  NumericVector ans(x.size(), NA_REAL) ;

  if( (lambda<0)&(theta/-lambda < 2)){
    //Rcout<< "Invalid genpois parameters" << std::endl;
    return ans;
  }

  ans = -x*lambda - theta + (x-1)*log( theta + x*lambda ) +
    log( theta ) - lgamma( x + 1 );


  if( !log_ret ){
    ans = exp( ans );
  }

  return ans;

}

// Same as above but for single observation x
double dgenpois_cpp( double x, double lambda, double theta, bool log_ret ){

  double ans =  NA_REAL;

  if( (lambda<0)&(theta/-lambda < 2)){
    //Rcout<< "Invalid genpois parameters" << std::endl;
    return ans;
  }

  ans = -x*lambda - theta + (x-1)*log( theta + x*lambda ) +
    log( theta ) - lgamma( x + 1 );


  if( !log_ret ){
    ans = exp( ans );
  }

  return ans;

}


// Compute the pmf for the genarlized Poisson (with some truncation)
// [[Rcpp::export]]
NumericVector dgp_comp_cpp( double lambda, double theta ){


  int m = 1e9;
  double k = 10.0;
  double M = 1/(1-lambda);
  if( lambda < 0 ){
    m = floor( theta/-lambda );
  }

  if( (m > 1e3)||(m<0)){
    m = ceil(theta*M + k*sqrt( theta*pow(M,3)) );
    double eps = 1e-6;
    double curr_val = dgenpois_cpp( m, lambda, theta, false );
    while( curr_val > eps ){
      k++;
      m = ceil(theta*M + k*sqrt( theta*pow(M,3)) );
      curr_val = dgenpois_cpp( m, lambda, theta, false );
    }
  }

  IntegerVector to_comp_int= seq_len( m+1 )-1;
  NumericVector to_comp = as<NumericVector>(to_comp_int);
  NumericVector out = dgenpois_cpp(to_comp, lambda, theta, false);
  if( lambda < 0 ){ out = out/sum(out); }
  return out;

}

/* Compute the likelihood for observed daughter counts x. mu is the integrals over the observation window
 * for each parent. dgp is the generalized Poisson pmf.
 *
 * For a parent to have x observed daughters, it must have at least x daughters in total. But, it could also have
 * x observed and, say, z unobserved. Given that it has x+z in total, the likelihood of observing x is
 * dbinom(x, x+z, mu), i.e. a binomial with n=x+z, p=mu. To get the unconditional (not conditioned on the total,
 * but only on the location, that is) probability of observing x we must take
 * sum( p(x+k in total)*p(only x observed) ). That is what this function does.
 */

// [[Rcpp::export]]
NumericVector dgenpoisbin_cpp( IntegerVector x, double lambda, double theta, NumericVector mu, NumericVector dgp ){

  NumericVector out( x.size ());
  //Rcout<< max(x) << ' ' << dgp.size() << std::endl;
  for( int i = 0; i<mu.size();i++){
    //Rcout<< i << std::endl;
    for( int j = x(i); j<dgp.size(); j++){
      //Rcout<< j << std::endl;
      out(i) += dgp(j)*R::dbinom( x(i),j,mu(i), 0);
    }
  }

  return( out );
}

// Same as above, but for single observation.
double dgenpoisbin_cpp( int x, double lambda, double theta, double mu, NumericVector dgp ){

  double out = 0;


  for( int j = x; j<dgp.size(); j++){
    out += dgp(j)*R::dbinom( x,j,mu, 0);
  }

  return( out );
}


// Compute the kernel integral of the observation window for parents c
// [[Rcpp::export]]
NumericVector intfun_cols_cpp(NumericMatrix c, double omega, NumericVector xr, NumericVector yr){

  NumericVector xc = c(_,0);
  NumericVector yc = c(_,1);
  arma::vec xl = 1-pnorm(xc,xr(0),omega);
  arma::vec xu = 1-pnorm(xc,xr(1),omega, FALSE);
  arma::vec yl = 1-pnorm(yc,yr(0),omega);
  arma::vec yu = 1-pnorm(yc,yr(1),omega, FALSE);

  return( arma2nvec( 1-xl-xu-yl-yu+xl%yl+xl%yu+xu%yl+xu%yu ) );
  // return( wrap( 1-xl + xu%yu ) );


}


// Same as above, but for single parent.
// [[Rcpp::export]]
double intfun_cpp(NumericVector c, double omega, NumericVector xr, NumericVector yr){

  double xc = c(0);
  double yc = c(1);
  double xl = 1-R::pnorm(xc,xr(0),omega, 1, 0);
  double xu = 1-R::pnorm(xc,xr(1),omega, 0,0);
  double yl = 1-R::pnorm(yc,yr(0),omega, 1, 0);
  double yu = 1-R::pnorm(yc,yr(1),omega, 0,0);

  return( 1-xl-xu-yl-yu+xl*yl+xl*yu+xu*yl+xu*yu  );
  // return( wrap( 1-xl + xu%yu ) );


}
