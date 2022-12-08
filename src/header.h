#include <RcppArmadillo.h>
#define NDEBUG // added in v0.1.23 for MacOS builds
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace Eigen;
using Eigen::MatrixXd;


// HELPERS
NumericVector arma2nvec( arma::vec x);
IntegerVector arma2ivec( arma::vec x);
NumericVector arma2nvec( arma::uvec x);
IntegerVector arma2ivec( arma::uvec x);

NumericVector dgenpois_cpp( NumericVector x, double lambda, double theta, bool log_ret );

double dgenpois_cpp( double x, double lambda, double theta, bool log_ret );

NumericVector dgp_comp_cpp( double lambda, double theta );

NumericVector dgenpoisbin_cpp( IntegerVector x, double lambda, double theta, NumericVector mu, NumericVector dgp );

double dgenpoisbin_cpp( int x, double lambda, double theta, double mu, NumericVector dgp );

NumericVector intfun_cols_cpp(NumericMatrix c, double omega, NumericVector xr, NumericVector yr);

double intfun_cpp(NumericVector c, double omega, NumericVector xr, NumericVector yr);

double r_omega_cpp( double ostar, double omega, List X, List C, double lambda, double theta, NumericVector dgp);

double r_lt_cpp( double lstar, double tstar, double lambda, double theta,
                 double dl, List X, List C, double omega, NumericVector dgp,
                 NumericVector dgpstar);

double r_kappa_cpp( double kstar, double kappa, double m, double win_size);

void update_P_cpp( List& X, List& C, double omega, double lambda, double theta,
                   NumericVector mus,
                   NumericVector dgp);

void update_P_cpp_n( List X, List C, double omega, double lambda, double theta,
                     NumericVector mus,
                     NumericVector dgp, int n_upd);

void move_C_cpp( List X, List& C, double omega, double lambda, double theta,
                 double sd_c, NumericVector dgp );

void c_birth_cpp( List& X, List& C, double kappa, double omega, double lambda, double theta, double s_birth,
                  NumericVector dgp);

void c_death_cpp( List& X, List& C, double kappa, double omega,
                  double lambda, double theta, double s_birth,
                  NumericVector dgp, int kill_c, arma::vec R);

double update_omega( double omega, double somega, List X,
                     List& C, double lambda, double theta, NumericVector dgp);
NumericVector update_lt( double lambda, double theta, double stheta,
                         double dl, List X, List& C, double omega, NumericVector dgp);
double update_kappa( double kappa, double skappa, double m, double win_size);
