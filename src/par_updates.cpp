// #define NDEBUG 1
// #include <RcppEigen.h>
#include "math.h"
#include "header.h"
// #include <RcppArmadillo.h>
// #define NDEBUG 1

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace Eigen;
using Eigen::MatrixXd;

// Compute the log Hastings ratio for a proposal of omega
//[[Rcpp::export]]
double r_omega_cpp( double ostar, double omega, List X, List C, double lambda, double theta, NumericVector dgp,
                    double a_omega, double b_omega){

  double rat = 0;

  NumericVector P = X["P"];
  ostar = exp( ostar );
  // omega = exp( omega );

  NumericVector cx = C["x"];
  NumericVector cy = C["y"];

  NumericMatrix Cxy(cx.size(),2);
  Cxy(_,0) = cx;
  Cxy(_,1) = cy;

  NumericVector xx = X["x"];
  NumericVector xy = X["y"];

  NumericMatrix Xxy(xx.size(),2);
  Xxy(_,0) = xx;
  Xxy(_,1) = xy;


  NumericVector xr = X["xrange"];
  NumericVector yr = X["yrange"];

  NumericVector mus = C["mus"];
  NumericVector mus_star = intfun_cols_cpp( Cxy, ostar, xr, yr );

  IntegerVector O = C["n"];

  double n = xx.size();

  // double m = cx.size();

  NumericMatrix diff_mat(Xxy.nrow(),2);

  for( int i = 0; i<Xxy.nrow();i++){
    diff_mat(i,_) = Cxy(P(i),_)-Xxy(i,_);
  }

  ArrayXXd diff_mat_Xd = as<ArrayXXd>(diff_mat);
  ArrayXXd toSum = diff_mat_Xd*diff_mat_Xd;

  rat = (2*n-1)*log( omega/ostar) -
    .5*toSum.sum()*(1/(ostar*ostar)-1/(omega*omega)) +
  sum( log( dgenpoisbin_cpp(O, lambda, theta, mus_star, dgp)))-
  sum( log( dgenpoisbin_cpp(O, lambda, theta, mus, dgp)))+
  (a_omega - 1)*log( ostar/omega)+b_omega*(omega - ostar);

  return rat;

}

//Compute the joint log Hastings ratio for proposed lambda and theta
//[[Rcpp::export]]
double r_lt_cpp( double lstar, double tstar, double lambda, double theta,
                 double dl, List X, List C, double omega, NumericVector dgp,
                 NumericVector dgpstar, double l_lambda, double u_lambda,
                 double a_theta, double b_theta){

  IntegerVector n = C["n"];
  if( lstar <= -tstar/max(n)){
    return -1e25;
  }

  NumericVector xr = X["xrange"];
  NumericVector yr = X["yrange"];

  NumericVector mus = C["mus"];

  NumericVector ul = NumericVector::create( u_lambda, lambda+dl);
  NumericVector ulstar = NumericVector::create( u_lambda, lstar+dl);

  NumericVector ll = NumericVector::create( l_lambda, lambda-dl, -tstar/max(n));
  NumericVector llstar = NumericVector::create( l_lambda, lstar-dl, -theta/max(n));

  double IL = min( ul )- max( ll );
  double ILstar = min( ulstar ) - max( llstar );

  double out = sum( log( dgenpoisbin_cpp(n, lstar, tstar, mus, dgpstar)))-
    sum( log( dgenpoisbin_cpp(n, lambda, theta, mus, dgp)))+
    log( IL/ILstar) + log( tstar/theta) +
    (a_theta-1)*log( tstar/theta) + (theta-tstar)*b_theta;

  return out;
}

// Compute the log Hastings ratio for kappa
//[[Rcpp::export]]
double r_kappa_cpp( double kstar, double kappa, double m, double win_size, double a_kappa, double b_kappa){

  return (m+1)*log( kstar/kappa ) + win_size*(kappa-kstar) + (a_kappa-1)*log( kstar/kappa) + (kappa-kstar)*b_kappa;

}


// Update omega
//[[Rcpp::export]]
double update_omega( double omega, double somega, List X,
                     List& C, double lambda, double theta, NumericVector dgp,
                     double a_omega, double b_omega){

  NumericVector cx = C["x"];
  NumericVector cy = C["y"];
  NumericMatrix Cxy = cbind( cx, cy );

  NumericVector xr= X["xrange"];
  NumericVector yr = X["yrange"];

  double ostar = R::rnorm(log(omega), somega );
  double orat = r_omega_cpp(ostar, omega, X, C, lambda, theta, dgp, a_omega, b_omega );
  double U = R::runif(0,1);
  if( log(U) <= orat ){
    omega = exp( ostar );
    C["mus"]= intfun_cols_cpp(Cxy, omega, xr, yr);
  }
  return (omega);

}

// Jointly update lambda and theta
//[[Rcpp::export]]
NumericVector update_lt( double lambda, double theta, double stheta,
                         double dl, List X, List& C, double omega, NumericVector dgp,
                         double l_lambda, double u_lambda,
                         double a_theta, double b_theta){

  IntegerVector n = C["n"];

  double tstar = exp( R::rnorm(  log( theta ), stheta ));
  if( -tstar/max(n)<(lambda+dl)){
    double ll = max( NumericVector::create(l_lambda,lambda-dl,-tstar/max(n)));
    double ul = min( NumericVector::create(u_lambda,lambda+dl));
    double lstar = R::runif( ll, ul);
    // //Rcout<< "dgp" << std::endl;
    NumericVector dgpstar = dgp_comp_cpp(lstar, tstar );
    // //Rcout<< "dgp done" << std::endl;
    double ltrat = r_lt_cpp(lstar, tstar, lambda, theta, dl, X, C, omega, dgp, dgpstar,l_lambda, u_lambda, a_theta, b_theta);
    // //Rcout<< "dgp done" << std::endl;
    double U = log( R::runif(0,1));
    if( U<=ltrat ){
      lambda = lstar;
      theta = tstar;
      dgp = dgpstar;
    }
  }

  return (NumericVector::create(lambda, theta));
}

// Update kappa
//[[Rcpp::export]]
double update_kappa( double kappa, double skappa, double m, double win_size, double a_kappa, double b_kappa){

  double kstar = exp( R::rnorm( log(kappa), skappa));
  double krat = r_kappa_cpp(kstar, kappa, m, win_size, a_kappa, b_kappa);
  double U = log( R::runif(0,1));
  if( U<= krat){kappa = kstar;}
  return (kappa);

}

// Jointly
//[[Rcpp::export]]
NumericVector mcmc_lt( double lambda, double theta, double stheta, NumericVector mus,
                       double dl, IntegerVector n, double l_lambda, double u_lambda,
                       double a_theta, double b_theta){


  double tstar = exp( R::rnorm(  log( theta ), stheta ));
  if( -tstar/max(n)<(lambda+dl)){
    double ll = max( NumericVector::create(l_lambda,lambda-dl,-tstar/max(n)));
    double ul = min( NumericVector::create(u_lambda,lambda+dl));
    double lstar = R::runif( ll, ul);
    // //Rcout<< "dgp" << std::endl;
    NumericVector dgp = dgp_comp_cpp( lambda, theta );
    NumericVector dgpstar = dgp_comp_cpp(lstar, tstar );
    // //Rcout<< "dgp done" << std::endl;
    if( lstar <= -tstar/max(n)){
      return NumericVector::create(lambda, theta);
    }


    NumericVector ulstar = NumericVector::create( u_lambda, lstar+dl);
    NumericVector llstar = NumericVector::create( l_lambda, lstar-dl, -theta/max(n));

    double IL =ul-ll;
    double ILstar = min( ulstar ) - max( llstar );

    double out = sum( log( dgenpoisbin_cpp(n, lstar, tstar, mus, dgpstar)))-
      sum( log( dgenpoisbin_cpp(n, lambda, theta, mus, dgp)))+
      log( IL/ILstar) + log( tstar/theta) +
      (a_theta-1)*log( tstar/theta) + (theta-tstar)*b_theta;

    // //Rcout<< "dgp done" << std::endl;
    double U = log( R::runif(0,1));
    if( U<=out ){
      lambda = lstar;
      theta = tstar;

    }
  }



  return NumericVector::create(lambda, theta);
}
