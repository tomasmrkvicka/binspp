// #include "../inst/include/binspp.h"
// #include <Rcpp.h>

// // #define ARMA_DONT_USE_WRAPPER

#include <RcppArmadillo.h>

// [[Rcpp::depends("RcppArmadillo")]]

#include <limits>

#ifdef _OPENMP
#include <omp.h>
#else
// for machines with compilers void of openmp support
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#define omp_get_wtime()        0
#define EIGEN_DONT_PARALLELIZE
#endif

#define min(x,y) (((x) < (y)) ? (x) : (y))

// I prefer not to use namespaces for RcppArmadillo
// using namespace std;
// using namespace Rcpp;


// calculate likelihood part
// [[Rcpp::export]]
double logpXCbeta(arma::mat X, arma::mat CC, double alpha, double omega, double AreaW, double integral){

  int nXrows = X.n_rows;
  int nCCrows = CC.n_rows;
  arma::rowvec d;
  arma::vec res = arma::zeros(nXrows);
  arma::vec dummy;

  for(int i = 0; i< nXrows; i++){
      for(int j = 0; j< nCCrows; j++){
  	    d =  X.row(i) - CC.row(j);
  		dummy =  (d*trans(d));
          res[i] = res[i] + exp(-dummy[0]/(2*omega*omega));
      }
}
double lik = AreaW - alpha*integral + nXrows*log(alpha/(2*(M_PI)*omega*omega)) + sum(log(res));

  return(lik);
}


// calculate r1 and r2
//' @title Calculate parameters for Birth and Death Interaction likelihood functions.
//' @description Calculates r1 and r2, which serve as distance thresholds for Birth and Death Interaction likelihood functions.
//' @param theta A vector of theta1 and theta2.
//' @return A vector of r1 and r2.
//' @export
// [[Rcpp::export]]
arma::vec coeff(arma::vec theta){

  arma::vec r(2);
  double t1 = theta[0], t2 = theta[1], t3 = 0.5, R = 0;
  double a = (t2-R)*(t2-R)/(t1*t3*t3);
  double b = t3*t3*(t1 - 1);
  double d = 27*a*b*b;
  double c = pow(d + sqrt((d+2)*(d+2)-4) + 2,1.0/3.0);
  double deltar = 1/(3*b) * ( c/pow(2,1.0/3.0) + pow(2,1.0/3.0)/c + 1 );
  r[0] = a / pow(deltar,1.5)+ t2;
  r[1] = r[0] - sqrt(deltar);

  return(r);
}


// Birth step attraction-repulsion likelihood evaluation
// [[Rcpp::export]]
Rcpp::List BirthInteractionLik(arma::vec xprop, arma::mat CC, arma::rowvec rho0sum, arma::vec theta, arma::vec r){

  int count = CC.n_rows;
  double r1 = r[0], r2 = r[1], t1 = theta[0], t2 = theta[1], t3 = 0.5, R = 0;
  double distance;
  arma::rowvec logrb = arma::zeros<arma::rowvec>(count);
  for(int i = 0; i< count; i++){
  	distance = sqrt( (xprop[0] - CC(i,0))*(xprop[0] - CC(i,0)) + (xprop[1] - CC(i,1))*(xprop[1] - CC(i,1)) );
  	if(distance>3000){
      logrb[i] = log(1);
    }else if ( (distance>r1) && (distance<=3000) ){
      logrb[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
    }else if ( (distance>R) && (distance<=r1) ){
      logrb[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
  	}
	}

	arma::rowvec rhosum = rho0sum + logrb;
    rhosum.insert_cols(count,1);
    rhosum[count] = sum(trans(logrb));

    double likelihood = 0;
	for(int i = 0; i<count+1; i++){
    	likelihood = min(rhosum[i],2) + likelihood;
    }

  return Rcpp::List::create(Rcpp::Named("likelihood") = likelihood, Rcpp::Named("rhosum") = rhosum);
}



// Death step attraction-repulsion likelihood evaluation
// [[Rcpp::export]]
Rcpp::List DeathInteractionLik(arma::vec xprop, int Death, arma::mat CC, arma::rowvec rho0sum, arma::vec theta, arma::vec r){

  Death = Death - 1; // for C code matching index
  int count = CC.n_rows;
  double r1 = r[0], r2 = r[1], t1 = theta[0], t2 = theta[1], t3 = 0.5, R = 0;
  double distance;
  arma::rowvec logrd = arma::zeros<arma::rowvec>(count);
  for(int i = 0; i< count; i++){
  	distance = sqrt( (xprop[0] - CC(i,0))*(xprop[0] - CC(i,0)) + (xprop[1] - CC(i,1))*(xprop[1] - CC(i,1)) );
    if(distance>3000){
      logrd[i] = log(1);
    }else if ( (distance>r1) && (distance<=3000) ){
      logrd[i] = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
    }else if ( (distance>R) && (distance<=r1) ){
      logrd[i] = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
    }
	}

	arma::rowvec rhosum = rho0sum - logrd;
	rhosum.shed_col(Death);

	double likelihood = 0;
	for(int i = 0; i<count-1; i++){
    	likelihood = min(rhosum[i],2) + likelihood;
  }

  return Rcpp::List::create(Rcpp::Named("likelihood") = likelihood, Rcpp::Named("rhosum") = rhosum);
}



// evaluate unnormalized likelihood for auxiliary variable
// [[Rcpp::export]]
double pCClik(arma::vec thetaprop, arma::mat CC){

    int count = CC.n_rows;
  	arma::vec r = coeff(thetaprop);
    double r1 = r[0], r2 = r[1], t1 = thetaprop[0], t2 = thetaprop[1], t3 = 0.5, R = 0;


    // lhY and lhYp
    double distance;
    arma::mat resultCCp = arma::zeros(count,count);
    for(int i = 0; i < count; i++){
      for(int j = 0; j <= i; j++){
        distance = sqrt( (CC(i,0) - CC(j,0))*(CC(i,0) - CC(j,0)) + (CC(i,1) - CC(j,1))*(CC(i,1) - CC(j,1)) );

        if(distance>3000){
            resultCCp(i,j) = resultCCp(j,i) = log(1);
        }else if ( (distance>r1) && (distance<=3000) ){
            resultCCp(i,j) = resultCCp(j,i) = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
        }else if ( (distance>R) && (distance<=r1) ){
            resultCCp(i,j) = resultCCp(j,i) = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
        }else {
            resultCCp(i,j) = resultCCp(j,i) = 0;
        }

      }
    }

    arma::rowvec rhoCCpsum = sum(resultCCp);
    double lhCCp = 0;
    for(int i = 0; i<count; i++){
    	lhCCp = min(rhoCCpsum[i],2) + lhCCp;
    }


return(lhCCp);
}



// evaluate unnormalized likelihood for auxiliary variable
//' @title Evaluate unnormalized likelihood for auxiliary variable
//' @description Calculates the unnormalized likelihood for an auxiliary variable by evaluating pairwise interaction between points. The interaction thresholds are derived from the input theta vector.
//' @param thetaprop A vector of theta1 and theta2.
//' @param CC A coordinate matrix of points.
//' @return A list of the computed likelihood and a row vector of summed interaction between points.
//' @export
// [[Rcpp::export]]
Rcpp::List pCClik2(arma::vec thetaprop, arma::mat CC){

    int count = CC.n_rows;
	  arma::vec r = coeff(thetaprop);
    double r1 = r[0], r2 = r[1], t1 = thetaprop[0], t2 = thetaprop[1], t3 = 0.5, R = 0;


    // lhY and lhYp
    double distance;
    arma::mat resultCCp = arma::zeros(count,count);
    for(int i = 0; i < count; i++){
      for(int j = 0; j <= i; j++){
        distance = sqrt( (CC(i,0) - CC(j,0))*(CC(i,0) - CC(j,0)) + (CC(i,1) - CC(j,1))*(CC(i,1) - CC(j,1)) );

        if(distance>3000){
          resultCCp(i,j) = resultCCp(j,i) = log(1);
        }else if ( (distance>r1) && (distance<=3000) ){
          resultCCp(i,j) = resultCCp(j,i) = log( 1 + 1/(t3*t3*(distance-r2)*(distance-r2)) );
        }else if ( (distance>R) && (distance<=r1) ){
          resultCCp(i,j) = resultCCp(j,i) = log( t1 - (sqrt(t1)*(distance-t2)/(t2-R))*(sqrt(t1)*(distance-t2)/(t2-R)) );
        }else {
          resultCCp(i,j) = resultCCp(j,i) = 0;
        }

      }
    }

    arma::rowvec rhoCCpsum = sum(resultCCp);
    double lhCCp = 0;
    for(int i = 0; i<count; i++){
    	lhCCp = min(rhoCCpsum[i],2) + lhCCp;
    }

  return Rcpp::List::create(Rcpp::Named("likelihood") = lhCCp, Rcpp::Named("rhosum") = rhoCCpsum);
}



// calculate temporal likelihood part
// [[Rcpp::export]]
double logpCCtemp(arma::mat X, arma::mat CC, double alpha, double omega, double numratio, double AreaW, double integral){

  int nXrows = X.n_rows;
  int nCCrows = CC.n_rows;
  arma::rowvec d;
  arma::vec res = arma::zeros(nXrows);
  arma::vec dummy;

  for(int i = 0; i< nXrows; i++){
      for(int j = 0; j< nCCrows; j++){
  	    d =  X.row(i) - CC.row(j);
    		dummy =  (d*trans(d));
        res[i] = res[i] + numratio*exp(-dummy[0]/(2*omega*omega));
      }
  }
  double lik = AreaW - alpha*integral + nXrows*log(alpha/(2*(M_PI)*omega*omega)) + sum(log(res));

  return(lik);
}
