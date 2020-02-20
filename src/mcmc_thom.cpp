#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <ctime>
using namespace Rcpp;


// Point pattern class to handle data and moms

class pp{
  std::vector<std::vector<double> > xyz;
  int dim;
  int n;
  std::vector<std::vector<double> > bbox;
  double Area;
  double temp;
public:
  pp() {temp = -1;};
  ~pp(){};
  void setxy(NumericMatrix m){
    dim = m.ncol();
    xyz.resize(dim);
    n = m.nrow();
    for(int i = 0; i < dim; i++) {
      xyz.at(i).resize(n);
      for(int j =0; j < n; j++) xyz.at(i).at(j) = m(j,i);
    }
  }
  void setbb(NumericMatrix bb){
    bbox.resize(dim);
    Area = 1;
    for(int i=0; i < dim; i++) {
      bbox.at(i).push_back(bb(0,i));
      bbox.at(i).push_back(bb(1,i));
      Area *= bb(1,i)-bb(0,i);
    }
  }
  double runif(int j) {
    return Rcpp::runif(1, bbox.at(j).at(0), bbox.at(j).at(1))[0];
  }
  double getx(int i, int d) {return xyz.at(d).at(i);}
  void remove(int *j) {
    for(int i=0; i < dim; i++) xyz.at(i).erase(xyz.at(i).begin() + *j);
    n--;
  }
  void addxy(double x, double y){
    xyz.at(0).push_back(x);
    xyz.at(1).push_back(y);
    n++;
  }
  void movexy(int *i, double x, double y){
    xyz.at(0).at(*i) = x;
    xyz.at(1).at(*i) = y;
  }

  int size() {return n;}
  double area() {return Area;}

  double intfun(double x, double y, double s2) {
    // s2 = sqrt(s2);
    double xl = R::pnorm(bbox.at(0).at(0), x, s2, 1, 0);
    double xu = R::pnorm(bbox.at(0).at(1), x, s2, 0, 0);
    double yl = R::pnorm(bbox.at(1).at(0), y, s2, 1, 0);
    double yu = R::pnorm(bbox.at(1).at(1), y, s2, 0, 0);
    return 1.0 - xl - xu - yl - yu + xl*yl + xl*yu + xu*yl + xu*yu;
  }

  // convert to NumericMatrix
  NumericMatrix tomatrix() {
    NumericMatrix out(n, dim);
    for(int i = 0; i < dim; i++)
      for(int j = 0; j < n; j++)
        out(j,i) = xyz.at(i).at(j);
    return out;
  }


};



// general functions for BDM part
// update mom-pattern for one data pattern, in-place.
void bdm(pp pat, pp *mom, double kappa, double alpha, double omega) {
  // par's shouldn't be in log scale
  double p, d2, xksum, xn, yn, intold, intnew, temp;
  double Rm0 = runif(1)[0];
  double Rm1 = runif(1)[0];
  double Rm2 = runif(1)[0];
  int i,j,k;
  // the denominator
  std::vector<double > den(pat.size());
  // omega = exp(omega);
  // kappa = exp(kappa);
  // alpha = exp(alpha);
  //
  for(i=0; i < pat.size(); i++){
    for(j = 0; j < mom->size(); j++){
      d2 = pow( pat.getx(i,0)-mom->getx(j,0) , 2) + pow( pat.getx(i,1)-mom->getx(j,1) , 2);
      den.at(i) += exp(-d2/(2.0*omega*omega));
    }
  }
  // then we do the action:
  if(Rm0 < 0.5) { // birth-death
    if(Rm1 < 0.5) {// birth
      xn = mom->runif(0);
      yn = mom->runif(1);
      //Rprintf("(%f,%f)",xn,yn);
      intnew = pat.intfun(xn, yn, omega);
      //Rprintf("(%f)", intnew);
      xksum = 0;
      for(i = 0; i < pat.size(); i++){
        d2 = pow(pat.getx(i,0) - xn,2) + pow(pat.getx(i,1)-yn, 2); // not nice way
        temp = exp(-d2/(2.0*omega*omega));
        xksum += log(1 + temp/den.at(i));
      }
      p = log(kappa) - alpha * intnew + xksum - log(mom->size() + 1) + log(mom->area());
      if(log(Rm2) < p) {
        mom->addxy(xn, yn);
      }
      //Rprintf("b");
      //Rprintf("b(%f)[%f,%f]", p, xn, yn);
    }else{ // death
      j = sample(mom->size(), 1)[0] - 1;
      intold = pat.intfun(mom->getx(j,0), mom->getx(j,1), omega);
      xksum = 0;
      for(i = 0; i < pat.size(); i++){
        d2 = 0;
        for(k = 0; k < 2; k++) d2 += pow(pat.getx(i,k)-mom->getx(j,k), 2);
        temp = exp(-d2/(2.0*omega*omega));
        xksum += log(1 - temp/den.at(i));
      }
      p = -log(kappa) + alpha * intold + xksum - log(mom->area()) + log(mom->size());
      if(log(Rm2) < p) { // death
        if(mom->size()>2)
          mom->remove(&j);
      }
      //Rprintf("d");
      //Rprintf("d(%f)[%i]", p, j+1);
    }
  }else{ // move
    j = sample(mom->size(), 1)[0] - 1;
    double Rm = runif(1)[0];
    // new location
    xn = mom->runif(0);
    yn = mom->runif(1);
    intold = pat.intfun(mom->getx(j,0), mom->getx(j,1), omega);
    intnew = pat.intfun(xn, yn, omega);
    xksum = 0;
    // old to new change
    for(i=0; i < pat.size(); i++) {
      d2 = pow(pat.getx(i,0) - xn,2) + pow(pat.getx(i,1)-yn, 2); // not nice way
      temp = exp(-d2/(2.0*omega*omega));
      d2 = pow(pat.getx(i,0)-mom->getx(j,0), 2) + pow(pat.getx(i,1)-mom->getx(j,1), 2);
      temp -= exp(-d2/(2.0*omega*omega));
      xksum += log(1 + temp/den.at(i));
    }
    p = -alpha * (intnew - intold) + xksum;
    if(log(Rm) < p) {
      mom->movexy(&j, xn, yn);
    }
    //Rprintf("m");
    //Rprintf("m(%f)",p);
  }
}



/*  Main function */


// [[Rcpp::export]]
List mcmc_thom (List patlist, // data patterns
                List momslist, // initial guess for parents
                List pat_bboxes, // bounding boxes for data patterns
                List mom_bboxes, // bounding boxes for moms
                int nsim, // Number of iterations

                double kappa0,
                double a_kappa,
                double b_kappa,
                double proposal_s_kappa,

                double alpha0,
                double a_alpha,
                double b_alpha,
                double proposal_s_alpha,

                double omega0,
                double a_omega,
                double b_omega,
                double proposal_s_omega,

                char name

) {
  RNGScope scope;
  int i, j, l, iter;

  double a, b, c, p, star, ros, temp, intint, intint_star, low, upp, oo;

  int nsamples = patlist.size();

  NumericVector kappa(nsim);
  NumericVector alpha(nsim);
  NumericVector omega(nsim);

  // initialise
  kappa(0) = kappa0;
  alpha(0) = alpha0;
  omega(0) = omega0;

  // precompute. somewhat cumbersome, need to change from List to vector<vector>
  std::vector<pp *> moms(nsamples);
  std::vector<pp> pats(nsamples);

  SEXP lts;
  NumericMatrix y;

  // prepare mom pattern objects
  for(j = 0; j < nsamples; j++) {
    // moms
    lts = momslist(j);
    y = lts;
    moms.at(j) = new pp();
    moms.at(j)->setxy(y);
    lts = mom_bboxes(j);
    y = lts;
    moms.at(j)->setbb(y);
    // data
    lts = patlist(j);
    y = lts;
    pats.at(j).setxy(y);

    lts = pat_bboxes(j);
    y = lts;
    pats.at(j).setbb(y);
  }
  // group and subject indicators

  // For timing iterations
  time_t t = time(0);
  time_t ttmp;

  // main loop
  for(iter = 0; iter < nsim-1; iter++) {

    checkUserInterrupt();

    // update mother point pattern
    for(i = 0; i < nsamples; i++) {
      //Rprintf("\n");
      bdm(pats.at(i),
          moms.at(i),
          kappa(iter),
          alpha(iter),
          omega(iter));
    }

    // update subject level kappa

    star = exp( rnorm(1, log( kappa(iter) ), proposal_s_kappa)[0] );
    temp = 0;
    for(i = 0; i < nsamples; i++) {
      temp += moms.at(i)->area() * ( kappa(iter)-star ) + (moms.at(i)->size()+1) * log(star / kappa(iter));
    }
    p = temp + (a_kappa-1)*log( star/kappa(iter)) + (kappa(iter)-star)*b_kappa;
    kappa(iter + 1) = log(runif(1)[0]) < p ? star : kappa(iter);

    // subject level alphas

    star = exp( rnorm(1, log( alpha(iter) ), proposal_s_alpha )[0] );
    temp = 0;
    oo = omega(iter);

    for(j = 0; j < nsamples; j++){

      intint = 0;
      for(i = 0; i < moms.at(j)->size(); i++) {
        intint += pats.at(j).intfun(moms.at(j)->getx(i,0), moms.at(j)->getx(i,1), oo);
      }
      temp += (alpha(iter) - star) * intint + pats.at(j).size() * log( (star / alpha(iter)) );
    }
    p = temp + (a_alpha +1)*log( star/alpha(iter)) + (alpha(iter)-star)*b_alpha;
    alpha(iter + 1) = log(runif(1)[0]) < p ? star : alpha(iter);
    //Rprintf("%f[%f]\n", p, temp);

    //if(iter >= 130) Rprintf("\n");
    // subject level omegas


    star = exp( rnorm(1, log( omega(iter) ), proposal_s_omega )[0] );
    oo = omega(iter);
    a = oo/star;
    ros = 0;
    for(j = 0; j < nsamples; j++){
      intint = 0;
      intint_star = 0;
      for(i = 0; i < moms.at(j)->size(); i++) {
        intint      += pats.at(j).intfun(moms.at(j)->getx(i,0), moms.at(j)->getx(i,1), oo); //NOW integrating over data window
        intint_star += pats.at(j).intfun(moms.at(j)->getx(i,0), moms.at(j)->getx(i,1), star);//NOW integrating over data window
      }
      temp = 0;
      for(l = 0; l < pats.at(j).size(); l++) {
        low = 0;
        upp = 0;
        for(i = 0; i < moms.at(j)->size(); i++) {
          b = pow( moms.at(j)->getx(i,0) - pats.at(j).getx(l,0), 2) + pow( moms.at(j)->getx(i,1) - pats.at(j).getx(l,1), 2);
          c = exp(-b/(2.0 * oo*oo));
          low += c;
          upp += pow(c, a*a);
        }
        temp += log(upp/low) + 2*log(a);
      }
      ros += -alpha(iter+1) * (intint_star - intint) + temp;
    }
    p = ros + (a_omega +1)*log( star/omega(iter)) + (omega(iter)-star)*b_omega;
    omega(iter + 1) = log(runif(1)[0]) < p ? star : omega(iter);



    if (iter % (nsim/10) == 0)
    {
      ttmp = time(0)  - t;
      t = time(0);
      Rcout << name << "says:\n";
      Rcout <<  ttmp << " seconds for last " << nsim/10 << " iterations.\n";
      Rcout <<  nsim-(iter+2) << "iterations and estimated " << (10*ttmp*(nsim-(iter+2)))/nsim <<" seconds " << " remaining.\n";
      // std::cout << "\r"<< percent << "%" << " done";
      // std::cout.flush();

    }
  }// end of main loop

  Rprintf("\n");

  List momlist(nsamples);
  for(i=0; i < nsamples; i++) {
    momlist(i) = moms.at(i)->tomatrix();
  }
  return List::create(Named("moms") = momlist, Named("kappa")=kappa, Named("alpha")=alpha, Named("omega")=omega);

}
