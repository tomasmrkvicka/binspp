#include <Rcpp.h>

using namespace Rcpp;




// [[Rcpp::export]]
NumericVector addMat(const NumericMatrix & m1, const NumericMatrix & m2) {
  NumericVector addMatrix = m1 + m2;
  addMatrix.attr("dim") = Dimension(m1.nrow(), m1.ncol());
  return addMatrix;
}


// [[Rcpp::export]]
NumericVector multMat(const NumericMatrix & m1, const NumericMatrix & m2) {
  NumericVector multMatrix = m1 * m2;
  multMatrix.attr("dim") = Dimension(m1.nrow(), m1.ncol());
  return multMatrix;
}


// [[Rcpp::export]]
double intalphaC(const List & z_beta, const double & alpha, const Nullable<NumericVector> & alphabet, const List & Wpix)
{
  NumericMatrix auxWpix = as<NumericMatrix>(Wpix["v"]);
  double xstep = as<double>(Wpix["xstep"]);
  double ystep = as<double>(Wpix["ystep"]);
  
  NumericMatrix aux = alpha * auxWpix;
  
  if (alphabet.isNotNull()) {
    NumericVector aalphabet(alphabet);
    
    List z0;
    NumericMatrix v;
    double abet;
    
    for (int i = 0; i < aalphabet.length(); i++) {
      z0 = z_beta[i];
      v = as<NumericMatrix>(z0["v"]);
      abet = aalphabet[i];
      aux = as<NumericMatrix>(addMat(aux, abet * v));
    } // for aalphabet.length
  } // if alphabet.isNotNull
  
  NumericVector aux2 = exp(as<NumericVector>(aux)) * as<NumericVector>(auxWpix);
  double result = sum(aux2)*xstep*ystep;
  
  return(result);
}


// [[Rcpp::export]]
double aozC(const List & z, const double & alpha, const Nullable<NumericVector> & alphabet, const NumericVector & u)
{
  double s = alpha;
  if (alphabet.isNotNull()) {
    NumericVector aalphabet(alphabet);
    
    List z0;
    double xstep, ystep; // int length, vlength;
    NumericVector xrange, yrange, xcol, yrow;
    NumericMatrix v; // int vcols, vrows;
    double p1, p2; int p1int, p2int;
    double w, abet;
    for (int i = 0; i < aalphabet.length(); i++) {
      z0 = z[i];
      xstep       = as<double>(z0["xstep"]);
      ystep       = as<double>(z0["ystep"]);
      // length      = z0.length();
      xrange = as<NumericVector>(z0["xrange"]);
      yrange = as<NumericVector>(z0["yrange"]);
      xcol = as<NumericVector>(z0["xcol"]);
      yrow = as<NumericVector>(z0["yrow"]);
      v = as<NumericMatrix>(z0["v"]);
      // vlength = v.length();
      // vcols = xcol.length();
      // vrows = yrow.length();
      p1 = (u[1] - yrange[0] + ystep / 2) / ystep;
      p2 = (u[0] - xrange[0] + xstep / 2) / xstep;
      p1int = round(p1);
      p2int = round(p2);
      if (fabs(p1 - p1int) == 0.5) {
        p1int = 2.0 * round(0.5 * p1);
      }
      if (fabs(p2 - p2int) == 0.5) {
        p2int = 2.0 * round(0.5 * p2);
      }
      --p1int; --p2int;
      w = v(p1int, p2int);
      abet = aalphabet[i];
      s = s + abet * w;
    } // for aalphabet.length
  } // if alphabet.isNotNull
  double result = (double)(exp(s));
  return(result);
}


// [[Rcpp::export]]
double logpXCbetC(const NumericMatrix & Y, const NumericMatrix & CC, const List & z_alpha, const List & z_omega, const double & alpha, const Nullable<NumericVector> & alphabet, const double & omega, const Nullable<NumericVector> & omegabet, const double & AreaW, const double & integral)
{
  NumericVector auxAlpha(CC.nrow());
  NumericVector auxOmega(CC.nrow());
  
  for (int i = 0; i < CC.nrow(); i++) {
    auxAlpha(i) = aozC(z_alpha, alpha, alphabet, CC(i,_));
    auxOmega(i) = aozC(z_omega, omega, omegabet, CC(i,_));
  } // for CC.nrow()
  
  NumericVector auxLLL(Y.nrow());
  NumericVector auxDist2(Y.nrow());
  
  for (int i = 0; i < Y.nrow(); i++) {
    auxDist2 = (CC(_,0) - Y(i,0)) * (CC(_,0) - Y(i,0)) + (CC(_,1) - Y(i,1)) * (CC(_,1) - Y(i,1));
    auxLLL(i) = sum(auxAlpha / (2 * M_PI * auxOmega * auxOmega) * exp( -auxDist2 / (2 * auxOmega*auxOmega)));
  } // for Y.nrow()
  
  double result = AreaW - integral + sum(log(auxLLL));
  return(result);
}



// [[Rcpp::export]]
double KumulaVsechC(const NumericMatrix & CC, const List & z_alpha, const List & z_omega, const double & alpha, const Nullable<NumericVector> & alphabet, const double & omega, const Nullable<NumericVector> & omegabet, const NumericVector & x_left, const NumericVector & x_right, const NumericVector & y_bottom, const NumericVector & y_top)
{
  NumericVector auxAlpha(CC.nrow());
  NumericVector auxOmega(CC.nrow());
  
  for (int i = 0; i < CC.nrow(); i++) {
    auxAlpha(i) = aozC(z_alpha, alpha, alphabet, CC(i,_));
    auxOmega(i) = aozC(z_omega, omega, omegabet, CC(i,_));
  } // for CC.nrow()
  
  double S = 0;
  NumericVector auxC(CC.nrow());
  
  for (int i = 0; i < x_left.length(); i++) {
    
    for (int j = 0; j < CC.nrow(); j++) {
      auxC(j) = auxAlpha(j) * ( R::pnorm(x_right(i), CC(j,0), auxOmega(j), true, false) - R::pnorm(x_left(i), CC(j,0), auxOmega(j), true, false) ) * ( R::pnorm(y_top(i), CC(j,1), auxOmega(j), true, false) - R::pnorm(y_bottom(i), CC(j,1), auxOmega(j), true, false) );
    } // for CC.nrow()
    
    S = S + sum(auxC);
  } // for x_left.length()
  
  return(S);
}



// [[Rcpp::export]]
double PrioralphaC(const double & a, const double & Prior_alpha_mean, const double & Prior_alpha_SD)
{
  double res = R::dnorm(a, Prior_alpha_mean, Prior_alpha_SD, 0);
  return(res);
}

// [[Rcpp::export]]
double PrioromegaC(const double & o, const double & Prior_omega_mean, const double & Prior_omega_SD)
{
  double res = R::dnorm(o, Prior_omega_mean, Prior_omega_SD, 0);
  return(res);
}

// [[Rcpp::export]]
double PrioralphabetC(const Nullable<NumericVector> & a, const Nullable<NumericVector> & Prior_alphavec_SD)
{
  double res = 1;
  if (a.isNotNull() & Prior_alphavec_SD.isNotNull()) {
    NumericVector aa(a);
    NumericVector bb(Prior_alphavec_SD);
    for (int i = 0; i < aa.length(); i++) {
      res = res * R::dnorm(aa(i), 0, bb(i), 0);
    } // for aa.length
  } // if a.isNotNull
  return(res);
}

// [[Rcpp::export]]
double PrioromegabetC(const Nullable<NumericVector> & o, const Nullable<NumericVector> & Prior_omegavec_SD)
{
  double res = 1;
  if (o.isNotNull() & Prior_omegavec_SD.isNotNull()) {
    NumericVector oo(o);
    NumericVector bb(Prior_omegavec_SD);
    for (int i = 0; i < oo.length(); i++) {
      res = res * R::dnorm(oo(i), 0, bb(i), 0);
    } // for oo.length
  } // if o.isNotNull
  return(res);
}



// [[Rcpp::export]]
List StepbetC(const double & kappa, const List & z_alpha, const List & z_omega, 
              const double & alpha, const double & salpha, const Nullable<NumericVector> & alphabet, const Nullable<NumericVector> & salphabet, 
              const double & omega, const double & somega, const Nullable<NumericVector> & omegabet, const Nullable<NumericVector> & somegabet,
              const NumericMatrix & Y, const NumericMatrix & CC, const double & logP, const double & integral, const double & integralrho,
              const NumericVector & x_left, const NumericVector & x_right, const NumericVector & y_bottom, const NumericVector & y_top, const List & Wpix, const double & AreaW, const double & FScoef,
              const double & Prior_alpha_mean, const double & Prior_alpha_SD, const double & Prior_omega_mean, const double & Prior_omega_SD,
              const Nullable<NumericVector> & Prior_alphavec_SD, const Nullable<NumericVector> & Prior_omegavec_SD)
{
  // update alpha
  double Newalpha = R::rnorm(alpha, salpha);
  NumericVector Newalphabet;
  
  if (alphabet.isNotNull() & salphabet.isNotNull()) {
    NumericVector aalphabet(alphabet);
    NumericVector ssalphabet(salphabet);
    NumericVector auxNewalphabet(aalphabet.length());
    
    for (int i = 0; i < aalphabet.length(); i++) {
      auxNewalphabet(i) = R::rnorm(aalphabet(i), ssalphabet(i));
      // Rcout << auxNewalphabet(i) << "\n";
    } // for aalphabet.length
    
    Newalphabet = auxNewalphabet;
    
  } // if alphabet.isNotNull & salphabet.isNotNull
  
  double int2, integralalpha, logP1, Newkappa, logProbAcceptA, logProbAcceptO;
  
  if (alphabet.isNotNull()) { // alphabet was updated
    int2 = KumulaVsechC(CC, z_alpha, z_omega, Newalpha, Newalphabet, omega, omegabet, x_left, x_right, y_bottom, y_top);
    integralalpha = intalphaC(z_alpha, Newalpha, Newalphabet, Wpix);
    logP1 = logpXCbetC(Y, CC, z_alpha, z_omega, Newalpha, Newalphabet, omega, omegabet , AreaW, int2);
    Newkappa = exp(FScoef) / integralalpha * AreaW;
    
    logProbAcceptA = logP1 - logP + integralrho*kappa - integralrho*Newkappa + CC.nrow()*(log(Newkappa) - log(kappa))
      + log(PrioralphabetC(Newalphabet, Prior_alphavec_SD)) - log(PrioralphabetC(alphabet, Prior_alphavec_SD))
      // + log(PrioromegabetC(omegabet, Prior_omegavec_SD)) - log(PrioromegabetC(omegabet, Prior_omegavec_SD))
      // + log(PrioromegaC(omega, Prior_omega_mean, Prior_omega_SD)) - log(PrioromegaC(omega, Prior_omega_mean, Prior_omega_SD))
      + log(PrioralphaC(Newalpha, Prior_alpha_mean, Prior_alpha_SD)) - log(PrioralphaC(alpha, Prior_alpha_mean, Prior_alpha_SD));
      
  } else { // alphabet is NULL and was not updated
    int2 = KumulaVsechC(CC, z_alpha, z_omega, Newalpha, alphabet, omega, omegabet, x_left, x_right, y_bottom, y_top);
    integralalpha = intalphaC(z_alpha, Newalpha, alphabet, Wpix);
    logP1 = logpXCbetC(Y, CC, z_alpha, z_omega, Newalpha, alphabet, omega, omegabet , AreaW, int2);
    Newkappa = exp(FScoef) / integralalpha * AreaW;
    
    logProbAcceptA = logP1 - logP + integralrho*kappa - integralrho*Newkappa + CC.nrow()*(log(Newkappa) - log(kappa))
      // + log(PrioralphabetC(Newalphabet, Prior_alphavec_SD)) - log(PrioralphabetC(alphabet, Prior_alphavec_SD))
      // + log(PrioromegabetC(omegabet, Prior_omegavec_SD)) - log(PrioromegabetC(omegabet, Prior_omegavec_SD))
      // + log(PrioromegaC(omega, Prior_omega_mean, Prior_omega_SD)) - log(PrioromegaC(omega, Prior_omega_mean, Prior_omega_SD))
      + log(PrioralphaC(Newalpha, Prior_alpha_mean, Prior_alpha_SD)) - log(PrioralphaC(alpha, Prior_alpha_mean, Prior_alpha_SD));
  }
  
  double accept = R::runif(0,1);
  List output;
  
  if (logP1==R_NegInf) { // all cluster centers far away from an observed point, reject proposal
    output["kappa"] = kappa;
    output["alpha"] = alpha;
    output["alphabet"] = alphabet;
    output["omega"] = omega;
    output["omegabet"] = omegabet;
    output["logP"] = logP;
    output["integral"] = integral;
    output["alphaAccept"] = 0;
  } else { // acceptance probability is positive
    if (log(accept) < logProbAcceptA){ // accept proposal
      // Rcout << "accept update alpha\n";
      output["kappa"] = Newkappa;
      output["alpha"] = Newalpha;
      if (alphabet.isNotNull()) { // alphabet was updated
        output["alphabet"] = Newalphabet;
      } else { // alphabet is NULL and was not updated
        output["alphabet"] = alphabet;
      }
      output["omega"] = omega;
      output["omegabet"] = omegabet;
      output["logP"] = logP1;
      output["integral"] = int2;
      output["alphaAccept"] = 1;
    } else { // reject proposal
      // Rcout << "reject update alpha\n";
      output["kappa"] = kappa;
      output["alpha"] = alpha;
      output["alphabet"] = alphabet;
      output["omega"] = omega;
      output["omegabet"] = omegabet;
      output["logP"] = logP;
      output["integral"] = integral;
      output["alphaAccept"] = 0;
    }
  }
  
  
  // update omega
  
  double Newomega = R::rnorm(omega, somega);
  // Rcout << "Newomega: " << Newomega << "\n";
  NumericVector Newomegabet;

  if (omegabet.isNotNull() & somegabet.isNotNull()) {
    NumericVector oomegabet(omegabet);
    NumericVector ssomegabet(somegabet);
    NumericVector auxNewomegabet(oomegabet.length());

    for (int i = 0; i < oomegabet.length(); i++) {
      auxNewomegabet(i) = R::rnorm(oomegabet(i), ssomegabet(i));
    } // for oomegabet.length

    Newomegabet = auxNewomegabet;

  } // if omegabet.isNotNull & somegabet.isNotNull


  integralalpha = intalphaC(z_alpha, output["alpha"], output["alphabet"], Wpix);
  Newkappa = exp(FScoef) / integralalpha * AreaW;

  if (omegabet.isNotNull()) { // omegabet is NOT null
    int2 = KumulaVsechC(CC, z_alpha, z_omega, output["alpha"], output["alphabet"], Newomega, Newomegabet, x_left, x_right, y_bottom, y_top);
    logP1 = logpXCbetC(Y, CC, z_alpha, z_omega, output["alpha"], output["alphabet"], Newomega, Newomegabet , AreaW, int2);

    logProbAcceptO = logP1 - as<double>(output["logP"]) + integralrho*as<double>(output["kappa"]) - integralrho*Newkappa + CC.nrow()*(log(Newkappa) - log(as<double>(output["kappa"])))
      // + log(PrioralphabetC(Newalphabet, Prior_alphavec_SD)) - log(PrioralphabetC(alphabet, Prior_alphavec_SD))
      // + log(PrioralphaC(Newalpha, Prior_alpha_mean, Prior_alpha_SD)) - log(PrioralphaC(alpha, Prior_alpha_mean, Prior_alpha_SD))
      + log(PrioromegabetC(Newomegabet, Prior_omegavec_SD)) - log(PrioromegabetC(omegabet, Prior_omegavec_SD))
      + log(PrioromegaC(Newomega, Prior_omega_mean, Prior_omega_SD)) - log(PrioromegaC(omega, Prior_omega_mean, Prior_omega_SD));

  } else { // omegabet is NULL
    int2 = KumulaVsechC(CC, z_alpha, z_omega, output["alpha"], output["alphabet"], Newomega, omegabet, x_left, x_right, y_bottom, y_top);
    logP1 = logpXCbetC(Y, CC, z_alpha, z_omega, output["alpha"], output["alphabet"], Newomega, omegabet , AreaW, int2);

    logProbAcceptO = logP1 - as<double>(output["logP"]) + integralrho*as<double>(output["kappa"]) - integralrho*Newkappa + CC.nrow()*(log(Newkappa) - log(as<double>(output["kappa"])))
      // + log(PrioralphabetC(Newalphabet, Prior_alphavec_SD)) - log(PrioralphabetC(alphabet, Prior_alphavec_SD))
      // + log(PrioralphaC(Newalpha, Prior_alpha_mean, Prior_alpha_SD)) - log(PrioralphaC(alpha, Prior_alpha_mean, Prior_alpha_SD))
      // + log(PrioromegabetC(Newomegabet, Prior_omegavec_SD)) - log(PrioromegabetC(omegabet, Prior_omegavec_SD))
      + log(PrioromegaC(Newomega, Prior_omega_mean, Prior_omega_SD)) - log(PrioromegaC(omega, Prior_omega_mean, Prior_omega_SD));
  }

  accept = R::runif(0,1);
  output["logProbAcceptAlpha"] = logProbAcceptA;
  output["logProbAcceptOmega"] = logProbAcceptO;

  if (logP1==R_NegInf) { // all cluster centers far away from an observed point, reject proposal
    // do not do anything, list "output" already contains required information
    output["omegaAccept"] = 0;
  } else { // acceptance probability is positive
    if (log(accept) < logProbAcceptO){ // accept proposal
      // Rcout << "accept update omega\n";
      output["kappa"] = Newkappa;
      output["omega"] = Newomega;
      if (omegabet.isNotNull()) { // omgebet was updated
        output["omegabet"] = Newomegabet;
      } else { // omegabet is NULL and was not updated
        output["omegabet"] = omegabet;
      }
      output["logP"] = logP1;
      output["integral"] = int2;
      output["omegaAccept"] = 1;
    } else { // reject proposal
      // Rcout << "reject update omega\n";
      // do not do anything, list "output" already contains required information
      output["omegaAccept"] = 0;
    }
  }

  return(output);
}



// [[Rcpp::export]]
NumericMatrix row_add (const NumericMatrix& x, const NumericVector& extraRow) {
  
  NumericMatrix x2(Dimension(x.nrow()+1, x.ncol()));
  
  for (int i = 0; i < x.nrow(); i++) {x2.row(i) = x.row(i);}
  x2.row(x.nrow()) = extraRow;
  
  return x2;
}


// [[Rcpp::export]]
int rand_int (const int & min, const int & max) {
  
  int res = floor(R::runif(min,max+1));
  
  return res;
}





// [[Rcpp::export]]
List StepMovePointC(const double & kappa, const List & z_alpha, const List & z_omega, 
                    const double & alpha, const Nullable<NumericVector> & alphabet, 
                    const double & omega, const Nullable<NumericVector> & omegabet,
                    const NumericMatrix & Y, const NumericMatrix & CC, const double & logP, const double & integral, const double & integralrho,
                    const NumericVector & x_left, const NumericVector & x_right, const NumericVector & y_bottom, const NumericVector & y_top, const double & AreaW,
                    const NumericVector & NewCenter)
{
  double choose = R::runif(0,1);
  List output;
  
  // Rcout << choose << "\n";
  
  if (choose < 1.0/3){ // move a parent point
    // Rcout << "move a parent point\n";
    int Discard = rand_int(0,(CC.nrow()-1));
    NumericMatrix auxCC = clone(CC);
    double int2 = integral - KumulaVsechC(auxCC(Range(Discard,Discard),_), z_alpha, z_omega, alpha, alphabet, omega, omegabet, x_left, x_right, y_bottom, y_top);
    NumericMatrix CC1 = clone(CC);
    CC1(Discard,_) = NewCenter;
    int2 = int2 + KumulaVsechC(CC1(Range(Discard,Discard),_), z_alpha, z_omega, alpha, alphabet, omega, omegabet, x_left, x_right, y_bottom, y_top);
    double logP2 = logpXCbetC(Y, CC1, z_alpha, z_omega, alpha, alphabet, omega, omegabet , AreaW, int2);
    
    double accept = R::runif(0,1);
    if (log(accept) < (logP2 - logP)){ // accept move proposal
      // Rcout << "move proposal accepted\n";
      output["CC"] = CC1;
      output["logP"] = logP2;
      output["integral"] = int2;
      output["parentAccept"] = 1;
    } else { // reject move proposal
      // Rcout << "move proposal rejected\n";
      output["CC"] = CC;
      output["logP"] = logP;
      output["integral"] = integral;
      output["parentAccept"] = 0;
    }
  } else { // do not move a parent point
    choose = R::runif(0,1);
    // Rcout << choose << "\n";
    if ((choose < 0.5) | (CC.nrow() < 3)){ // birth of a parent point
      // Rcout << "birth a parent point\n";
      NumericMatrix CC1 = clone(CC);
      CC1 = row_add(CC1, NewCenter);
      double int2 = integral + KumulaVsechC(CC1(Range(CC1.nrow()-1,CC1.nrow()-1),_), z_alpha, z_omega, alpha, alphabet, omega, omegabet, x_left, x_right, y_bottom, y_top);
      double logP2 = logpXCbetC(Y, CC1, z_alpha, z_omega, alpha, alphabet, omega, omegabet , AreaW, int2);
      
      double accept = R::runif(0,1);
      if (log(accept) < (logP2 - logP + log(kappa*integralrho) - log(CC1.nrow()))){ // accept birth proposal
        // Rcout << "birth proposal accepted\n";
        output["CC"] = CC1;
        output["logP"] = logP2;
        output["integral"] = int2;
        output["parentAccept"] = 1;
      } else { // reject birth proposal
        // Rcout << "birth proposal rejected\n";
        output["CC"] = CC;
        output["logP"] = logP;
        output["integral"] = integral;
        output["parentAccept"] = 0;
      }
    } else { // death of a parent point
      // Rcout << "death a parent point\n";
      int Discard = rand_int(0,(CC.nrow()-1));
      NumericMatrix auxCC2 = clone(CC);
      double int2 = integral - KumulaVsechC(auxCC2(Range(Discard,Discard),_), z_alpha, z_omega, alpha, alphabet, omega, omegabet, x_left, x_right, y_bottom, y_top);
      
      // remove the appropriate row of the matrix
      NumericMatrix CC1(Dimension(CC.nrow() - 1, CC.ncol()));
      int iter = 0; 
      for (int i = 0; i < CC.nrow(); i++) {
        if (i != Discard) {
          CC1.row(iter) = CC.row(i);
          iter++;
        }
      }
      
      double logP2 = logpXCbetC(Y, CC1, z_alpha, z_omega, alpha, alphabet, omega, omegabet , AreaW, int2);
      
      if (logP2==R_NegInf) { // all cluster centers far away from an observed point, reject death proposal
        // Rcout << "death proposal rejected\n";
        output["CC"] = CC;
        output["logP"] = logP;
        output["integral"] = integral;
        output["parentAccept"] = 0;
      } else {
        double accept = R::runif(0,1);
        if (log(accept) < (logP2 - logP - log(kappa*integralrho) + log(CC1.nrow()))){ // accept birth proposal
          // Rcout << "death proposal accepted\n";
          output["CC"] = CC1;
          output["logP"] = logP2;
          output["integral"] = int2;
          output["parentAccept"] = 1;
        } else { // reject death proposal
          // Rcout << "death proposal rejected\n";
          output["CC"] = CC;
          output["logP"] = logP;
          output["integral"] = integral;
          output["parentAccept"] = 0;
        }
      }
    }
  } // end move / do not move a parent point
  
  return output;
  
}

/*** R

*/

