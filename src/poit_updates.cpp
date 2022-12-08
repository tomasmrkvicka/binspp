// #include <RcppEigen.h>
#include "math.h"
#include "header.h"
// #include <RcppArmadillo.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Update parent--daughter connections.
//[[Rcpp::export]]
double update_P_cpp( List& X, List& C, double omega, double lambda, double theta,
                   NumericVector mus,
                   NumericVector dgp, double& p_count){

  NumericVector P = X["P"];
  NumericVector n = C["n"];

  NumericVector xx = X["x"];
  NumericVector xy = X["y"];
  NumericMatrix Xxy = cbind( xx, xy);

  NumericVector cx = C["x"];
  NumericVector cy = C["y"];
  NumericMatrix Cxy = cbind( cx, cy);

  int switch_idx = floor( R::runif(0, 1)*xx.size());
  int new_P = floor( R::runif(0, 1)*cx.size());
  int old_P = P(switch_idx);

  NumericVector mus_switch = NumericVector::create(mus(old_P), mus(new_P));
  NumericVector n_switch= NumericVector::create( n(old_P), n(new_P));

  NumericVector P_star(clone( P ));
  P_star( switch_idx ) = new_P;

  IntegerVector n_old = IntegerVector::create(n(old_P),n(new_P));

  IntegerVector n_star(clone( n_old ));
  n_star = n_star + IntegerVector::create( -1, 1);

  if( (lambda<0)&(n_star(1)>floor( theta/-lambda))){
    return p_count;
  }


  double diffNew =  ( Xxy(switch_idx,1)-Cxy(new_P,1))*(Xxy(switch_idx,1)-Cxy(new_P,1))+
    ( Xxy(switch_idx,0)-Cxy(new_P,0))*(Xxy(switch_idx,0)-Cxy(new_P,0));


  double diffOld =  ( Xxy(switch_idx,1)-Cxy(old_P,1))*(Xxy(switch_idx,1)-Cxy(old_P,1))+
    ( Xxy(switch_idx,0)-Cxy(old_P,0))*(Xxy(switch_idx,0)-Cxy(old_P,0));



  double rat = -.5/(omega*omega)*( diffNew - diffOld) +
    sum( log( dgenpoisbin_cpp( n_star, lambda, theta, mus_switch, dgp)))+
    sum( lgamma( n_star + 1)) -
    sum( log( dgenpoisbin_cpp( n_old, lambda, theta, mus_switch, dgp)))-
    sum( lgamma( n_old+1));



  double U = R::runif(0,1);
  if( (log(U)<=rat)&(old_P!=new_P)){

    n(old_P) -=1;
    n(new_P) +=1;
    C["n"]= n;
    X["P"] = P_star;
    p_count +=1;
  }

  return p_count;

}

// Same as above, but for single parent
//[[Rcpp::export]]
double update_P_cpp_n( List& X, List& C, double omega, double lambda, double theta,
                     NumericVector mus,
                     NumericVector dgp, int n_upd, double& p_count){

  for( int i = 0; i < n_upd; i++){
    p_count =  update_P_cpp(X, C, omega, lambda, theta, mus, dgp, p_count );
  }
  return p_count;

}

// MH parent moves.
//[[Rcpp::export]]
double move_C_cpp( List X, List& C, double omega, double lambda, double theta,
                 double sd_c, NumericVector dgp, double& m_count ){


  NumericVector xr = X["xrange"];
  NumericVector yr = X["yrange"];

  NumericVector cxr = C["xrange"];
  NumericVector cyr = C["yrange"];

  NumericVector cx = C["x"];
  NumericVector cy = C["y"];

  NumericVector xx = X["x"];
  NumericVector xy = X["y"];

  int move_id = floor( R::runif(0, 1)*cx.size());

  NumericVector Cxy(2);
  Cxy(0) = cx(move_id);
  Cxy(1) = cy(move_id);

  NumericVector mus = C["mus"];
  NumericVector Os = C["n"];
  double O = Os(move_id);
  double mu = mus(move_id);
  double muc = intfun_cpp( Cxy, sd_c, cxr, cyr);

  NumericVector Cxy_star( 2 );
  Cxy_star(0) = 2*cxr(1);
  Cxy_star(1) = 2*cyr(1);

  while( Cxy_star(0)>cxr(1)||Cxy_star(0)<cxr(0)||
         Cxy_star(1)>cyr(1)||Cxy_star(1)<cyr(0)){

    Cxy_star = rnorm( 2, 0,sd_c)+Cxy;

  }

  double mu_star = intfun_cpp( Cxy_star, omega, xr, yr);
  double muc_star = intfun_cpp( Cxy_star, sd_c, cxr, cyr);

  arma::vec P = X["P"];
  IntegerVector ids = arma2ivec(arma::find( P == move_id));


  double rat = log( muc/muc_star) +log( dgenpoisbin_cpp(0, lambda, theta, mu_star, dgp ) /
     dgenpoisbin_cpp(0, lambda, theta, mu, dgp ) ) ;


  if( ids.size()>0){
    double new_dist = 0, old_dist = 0;

    for( int i = 0; i < ids.size(); i++){

      new_dist += (xx(ids(i))-Cxy_star(0))*(xx(ids(i))-Cxy_star(0)) +
        (xy(ids(i))-Cxy_star(1))*(xy(ids(i))-Cxy_star(1));

      old_dist += (xx(ids(i))-Cxy(0))*(xx(ids(i))-Cxy(0)) +
        (xy(ids(i))-Cxy(1))*(xy(ids(i))-Cxy(1));

    }

    rat = -.5/(omega*omega)*( new_dist - old_dist) +
      log( dgenpoisbin_cpp(O, lambda, theta, mu_star, dgp )/
        dgenpoisbin_cpp(O, lambda, theta, mu, dgp ))+
          log( muc/muc_star);
  }

  double U = Rf_runif(0,1);
  if( log(U)<=rat){
    cx(move_id) = Cxy_star(0);
    C["x"] = cx;
    cy(move_id) = Cxy_star(1);
    C["y"] = cy;
    mus(move_id) = mu_star;
    C["mus"] = mus;
    m_count += 1;
    // //Rcout<< "move" << std::endl;
  }
  return m_count;
}



// MH add parent
//[[Rcpp::export]]

double c_birth_cpp( List& X, List& C, double kappa, double omega, double lambda, double theta, double s_birth,
                  NumericVector dgp, double& b_count){


  IntegerVector n = C["n"];
  NumericVector xrc =  C["xrange"];
  NumericVector yrc = C["yrange"];
  NumericVector xr = X["xrange"];
  NumericVector yr = X["yrange"];
  NumericVector mus = C["mus"];
  IntegerVector P = X["P"];

  NumericVector xx = X["x"];
  NumericVector xy = X["y"];
  // NumericMatrix Xxy = cbind( xx, xy);

  NumericVector cx = C["x"];
  NumericVector cy = C["y"];
  NumericMatrix Cxy = cbind( cx, cy);

  NumericVector new_c(2);
  new_c(0) = R::runif(xrc(0),xrc(1));
  new_c(1) = R::runif(yrc(0),yrc(1));

  NumericVector new_dist( xx.size());
  NumericVector old_dist( xx.size());
  arma::vec dist_rat(xx.size());


  for( int i = 0; i < xx.size(); i++ ){

    new_dist(i) = exp( -.5*(  (xx(i)-new_c(0))*(xx(i)-new_c(0))+
      (xy(i)-new_c(1))*(xy(i)-new_c(1)))/(s_birth*s_birth) );

    old_dist(i) = exp( -.5*(  (xx(i)-cx(P(i)))*(xx(i)-cx(P(i)))+
      (xy(i)-cy(P(i)))*(xy(i)-cy(P(i))))/(s_birth*s_birth) );


    dist_rat(i) = new_dist(i)/(new_dist(i)+old_dist(i));

  }

  arma::vec U = Rcpp::runif( xx.size() );
  IntegerVector steal = arma2ivec(find( U <= dist_rat ));
  IntegerVector keep = arma2ivec(find( U > dist_rat ));
  // arma::uvec steal= find( U <= dist_rat );
  // arma::uvec keep = find( U > dist_rat  );

  if( (lambda<0)&(steal.size()>theta/-lambda) ){return b_count;}
  if( steal.size()>(dgp.size()-1)){return b_count;}

  IntegerVector P_star(clone(P));
  // P_star[steal] = Cxy.nrow();

  NumericVector cx_star = cx;
  cx_star.push_back(new_c(0));
  NumericVector cy_star = cy;
  cy_star.push_back(new_c(1));
  NumericMatrix Cxy_star = cbind( cx_star, cy_star);

  NumericVector mus_star = clone(mus);
  mus_star.push_back(intfun_cpp(new_c, omega, xr, yr));

  IntegerVector n_star = clone(n);
  for( int i = 0; i<steal.size(); i++ ){

    n_star(P(steal(i))) -= 1;
    P_star(steal(i)) = Cxy.nrow();
  }

  n_star.push_back(steal.size());
  NumericVector dist_rat0 = arma2nvec( dist_rat);

  double steal_prob = -sum( log( 1-pmin(dist_rat0, 1.0)));
  double mutmp = mus_star(mus_star.size()-1);
  //
  double rat = log( dgenpoisbin_cpp(0, lambda, theta, mutmp, dgp ) )+
    log( kappa ) - log( Cxy_star.nrow()) +
    log( (xrc(1)-xrc(0))*(yrc(1)-yrc(0)) ) + steal_prob;

  if( steal.size()>0){
    NumericMatrix dist_pmf( Cxy.nrow(), steal.size());

    for( int j = 0; j < steal.size(); j++){
      for( int i = 0; i < Cxy.nrow(); i++ ){

        dist_pmf(i,j) = exp(-.5*((cx(i)-xx(steal(j)))*(cx(i)-xx(steal(j)))+
          (cy(i)-xy(steal(j)))*(cy(i)-xy(steal(j))))/(s_birth*s_birth));

      }
      dist_pmf(_,j) = dist_pmf(_,j)/sum(dist_pmf(_,j));
    }

    // NumericVector steal_dist = as<NumericVector>(wrap(dist_rat.elem(steal)));
    // NumericVector keep_dist = as<NumericVector>(wrap(dist_rat.elem(keep)));

    NumericVector steal_dist = dist_rat0[steal];
    NumericVector keep_dist = dist_rat0[keep];

    steal_prob = -sum( log( pmin( 1, steal_dist))) - sum( log( pmin( 1, 1-keep_dist)));
    IntegerVector steal_pars =  P[steal];
    IntegerVector un_steal_pars = unique( steal_pars );

    double dist_diff = 0;

    for( int i = 0; i < xx.size(); i++){

      dist_diff += -.5*( (xx(i)-cx_star(P_star(i)))*(xx(i)-cx_star(P_star(i)))+
        (xy(i)-cy_star(P_star(i)))*(xy(i)-cy_star(P_star(i)))-
        (xx(i)-cx(P(i)))*(xx(i)-cx(P(i)))-
        (xy(i)-cy(P(i)))*(xy(i)-cy(P(i))) )/(omega*omega);

    }

    IntegerVector new_n =  n_star[un_steal_pars];
    new_n.push_back( n_star(cx.size()));
    IntegerVector old_n =  n[un_steal_pars];

    NumericVector new_mus = mus_star[un_steal_pars];
    new_mus.push_back( mus_star(cx.size()));
    NumericVector old_mus = mus[un_steal_pars];

    NumericVector pmfs( steal.size());
    for(int i = 0; i<steal.size(); i++){
      pmfs(i) = dist_pmf( P(steal(i)),i);
    }

    // //Rcout<< old_n << std::endl;

    rat =
      dist_diff +
      sum( log( dgenpoisbin_cpp(new_n, lambda, theta, new_mus, dgp )))+
      sum( lgamma( new_n+1))-
      sum( log( dgenpoisbin_cpp(old_n, lambda, theta, old_mus, dgp )))-
      sum( lgamma( old_n+1))+
      log( kappa )-
      log( cx_star.size()) +
      sum(log(pmfs))+
      log( (xrc(1)-xrc(0))*(yrc(1)-yrc(0)) ) + steal_prob;
  }
  double R = log( R::runif(0,1));
  if( R <= rat ){
    C["x"] = cx_star;
    C["y"] = cy_star;
    C["n"] = n_star;
    C["mus"] = mus_star;
    X["P"] = P_star;
    b_count += 1;

  }
  return b_count;
}

// MH kill parent
//[[Rcpp::export]]
double c_death_cpp( List& X, List& C, double kappa, double omega,
                  double lambda, double theta, double s_birth,
                  NumericVector dgp, double& d_count){
  NumericVector xrc =  C["xrange"];
  NumericVector yrc = C["yrange"];
  NumericVector xr = X["xrange"];
  NumericVector yr = X["yrange"];
  std::vector<double> mus = C["mus"];
  NumericVector xx = X["x"];
  NumericVector xy = X["y"];
  NumericMatrix Xxy = cbind( xx, xy);
  std::vector<double> cx = C["x"];
  std::vector<double> cy = C["y"];
  // NumericMatrix Cxy = cbind( cx, cy);

  int kill_c = floor( R::runif(0, 1)*cx.size());

  arma::vec P = X["P"];
  arma::vec P_star =  P;
  arma::uvec P_shift = find( P_star > kill_c);
  // IntegerVector P_shift = wrap(find( P_star > kill_c ));
  for( unsigned int i = 0; i < P_shift.size(); i++){
    P_star(P_shift(i)) -= 1;
  }

  std::vector<double> cx_star =cx;
  cx_star.erase(cx_star.begin() + kill_c);
  std::vector<double> cy_star = cy;
  cy_star.erase(cy_star.begin() + kill_c);
  std::vector<double> mus_star = mus;
  mus_star.erase(mus_star.begin()+kill_c);
  std::vector<int> n = C["n"];
  std::vector<int> n_star = n;
  n_star.erase( n_star.begin() + kill_c);

  // Causes tracemem error for some reason
  // IntegerVector move_x = as<IntegerVector>(wrap( find( P == kill_c ) ));
  // IntegerVector stay_x = as<IntegerVector>(wrap( find( P != kill_c ) ));

  arma::uvec mx = find( P == kill_c );
  arma::uvec sx = find( P != kill_c );
  IntegerVector move_x = IntegerVector(mx.begin(), mx.end());
  IntegerVector stay_x = IntegerVector(sx.begin(), sx.end());

  IntegerVector new_parents( move_x.size());
  IntegerVector old_parents( move_x.size());
  // //Rcout<< move_x.size() << std::endl;
  NumericMatrix dist_pmf( cx_star.size(), move_x.size());
  NumericMatrix dist_cdf( cx_star.size(), move_x.size());

  if( move_x.size()>0){

    for( int j = 0; j < move_x.size(); j++){
      for( unsigned int i = 0; i < cx_star.size(); i++ ){

        dist_pmf(i,j) = exp(-.5*((cx_star[i]-xx(move_x(j)))*(cx_star[i]-xx(move_x(j)))+
          (cy_star[i]-xy(move_x(j)))*(cy_star[i]-xy(move_x(j))))/(s_birth*s_birth));

      }
      dist_pmf(_,j) = dist_pmf(_,j)/sum(dist_pmf(_,j));
    }


    dist_cdf(0,_) = dist_pmf( 0,_ );

    for( int i = 1; i<dist_cdf.nrow(); i++){

      dist_cdf( i,_ ) = dist_cdf( i-1,_) + dist_pmf(i,_);

    }
    arma::vec R = runif( move_x.size());
    int j = 0;
    for( int i = 0; i<move_x.size(); i++ ){
      j = 0;
      while( dist_cdf(j,i)<= R(i) ){
        j++;
      }
      new_parents(i) = j;
      old_parents(i) = j;
      if(j>=kill_c){old_parents(i) = j+1;}
      n_star[j] += 1;
      P_star(move_x(i)) = j;
    }

  }

  std::vector<double> new_dist( xx.size());
  std::vector<double> old_dist( xx.size());
  NumericVector dist_rat( xx.size());

  for( int i = 0; i<xx.size(); i++ ){
    new_dist[i] = exp( -.5/(s_birth*s_birth)*(
      (xx(i)-cx[kill_c])*(xx(i)-cx[kill_c]) +
      (xy(i)-cy[kill_c])*(xy(i)-cy[kill_c])
    ));

    old_dist[i] = exp( -.5/(s_birth*s_birth)*(
      (xx(i)-cx_star[P_star[i]])*(xx(i)-cx_star[P_star[i]]) +
      (xy(i)-cy_star[P_star[i]])*(xy(i)-cy_star[P_star[i]])
    ));

    dist_rat(i) = new_dist[i]/(new_dist[i]+old_dist[i]);

  }
  double steal_prob = sum( log( 1-pmin(dist_rat, 1)));
  // //Rcout<< move_x << '\n' << steal_prob <<  std::endl;

  if( move_x.size()>0){

    NumericVector move_rat = dist_rat[move_x];
    NumericVector stay_rat = dist_rat[stay_x];
    steal_prob = sum( log( pmin( 1, move_rat)) ) + sum( log(1-pmin( 1, stay_rat)));
  }

  if( (lambda<0)&(*max_element(n_star.begin(), n_star.end())>floor( theta/-lambda))){ return d_count;}

  double rat = 0;
  if( move_x.size() == 0){
    rat = -log( dgenpoisbin_cpp(0, lambda, theta, mus[kill_c], dgp))-
      log( kappa ) +
      log( cx_star.size()+1)-
      log( (xrc(1)-xrc(0))*(yrc(1)-yrc(0))) +
      steal_prob;
    // //Rcout<<  steal_prob <<  std::endl;

  }
  if( move_x.size()>0){
    double dist_new = 0, dist_old = 0;

    for( int i = 0; i < move_x.size(); i++){

      dist_new += (xx[move_x(i)]-cx_star[new_parents(i)])*(xx[move_x(i)]-cx_star[new_parents(i)])+
        (xy[move_x(i)]-cy_star[new_parents(i)])*(xy[move_x(i)]-cy_star[new_parents(i)]);
      dist_old += (xx[move_x(i)]-cx[kill_c])*(xx[move_x(i)]-cx[kill_c])+
        (xy[move_x(i)]-cy[kill_c])*(xy[move_x(i)]-cy[kill_c]);
    }

    IntegerVector new_pars = unique( new_parents );
    IntegerVector old_pars = unique( old_parents );
    old_pars.push_back(kill_c);

    IntegerVector n_new( new_pars.size());
    NumericVector mus_new( new_pars.size() );
    IntegerVector n_old( old_pars.size());
    NumericVector mus_old( old_pars.size() );

    for( int i = 0; i < new_pars.size(); i++ ){
      n_new(i) = n_star[new_pars(i)];
      mus_new(i) = mus_star[new_pars(i)];
    }
    for( int i = 0; i < old_pars.size(); i++ ){
      n_old(i) = n[old_pars(i)];
      mus_old(i) = mus[old_pars(i)];
    }

    double pmf_tmp = 0;
    for( int i =0; i<new_parents.size(); i++ ){ pmf_tmp += log( dist_pmf(new_parents(i),i) );}

    rat = -.5/(omega*omega)*( dist_new - dist_old) +
      sum( log( dgenpoisbin_cpp(n_new, lambda, theta, mus_new, dgp)))+
      sum( lgamma( n_new +1 ))-
      sum( log( dgenpoisbin_cpp(n_old, lambda, theta, mus_old, dgp)))-
      sum( lgamma( n_old +1 ))-
      log( kappa )+
      log( cx.size())-
      pmf_tmp-
      log( (xrc(1)-xrc(0))*(yrc(1)-yrc(0))) +
      steal_prob;
  }
  double U = R::runif(0,1);
  if( log(U)<= rat ){
    C["x"] = cx_star;
    C["y"] = cy_star;
    C["mus"] = mus_star;
    X["P"] = P_star;
    C["n"] = n_star;
    d_count += 1;
  }
  return d_count;
}


// MH birth-death-move
//[[Rcpp::export]]
void bdm( List& X,
          List& C,
          double ktmp,
          double otmp,
          double ltmp,
          double ttmp,
          double smove,
          double sbirth,
          int n_bdm,
          double& m_count,
          double& b_count,
          double& d_count){
  double U = 0;
  NumericVector dgp = dgp_comp_cpp(ltmp, ttmp );
  NumericVector U0 = Rcpp::runif( n_bdm );
  for( int j = 0; j < n_bdm; j++){

    U = U0(j);
    // double U = R::runif(0,1);
    if( U<=0.5){
      // //Rcout<< "move" << std::endl;
      m_count = move_C_cpp(X, C, otmp, ltmp, ttmp, smove, dgp, m_count);
    }
    if( (U>0.5)&(U<=0.75) ){
      // //Rcout<< "birth" << std::endl;
      b_count = c_birth_cpp(X, C, ktmp, otmp, ltmp, ttmp, sbirth, dgp, b_count);
    }
    if( U>0.75){
      // //Rcout<< "death" << std::endl;
      d_count = c_death_cpp(X, C, ktmp, otmp, ltmp, ttmp, sbirth, dgp, d_count);
    }
    // //Rcout<< j << std::endl;
  }

}
