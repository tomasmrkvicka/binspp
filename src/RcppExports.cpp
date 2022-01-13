// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// addMat
NumericVector addMat(const NumericMatrix& m1, const NumericMatrix& m2);
RcppExport SEXP _binspp_addMat(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(addMat(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// multMat
NumericVector multMat(const NumericMatrix& m1, const NumericMatrix& m2);
RcppExport SEXP _binspp_multMat(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(multMat(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// intalphaC
double intalphaC(const List& z_beta, const double& alpha, const Nullable<NumericVector>& alphabet, const List& Wpix);
RcppExport SEXP _binspp_intalphaC(SEXP z_betaSEXP, SEXP alphaSEXP, SEXP alphabetSEXP, SEXP WpixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type z_beta(z_betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< const List& >::type Wpix(WpixSEXP);
    rcpp_result_gen = Rcpp::wrap(intalphaC(z_beta, alpha, alphabet, Wpix));
    return rcpp_result_gen;
END_RCPP
}
// aozC
double aozC(const List& z, const double& alpha, const Nullable<NumericVector>& alphabet, const NumericVector& u);
RcppExport SEXP _binspp_aozC(SEXP zSEXP, SEXP alphaSEXP, SEXP alphabetSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(aozC(z, alpha, alphabet, u));
    return rcpp_result_gen;
END_RCPP
}
// logpXCbetC
double logpXCbetC(const NumericMatrix& Y, const NumericMatrix& CC, const List& z_alpha, const List& z_omega, const double& alpha, const Nullable<NumericVector>& alphabet, const double& omega, const Nullable<NumericVector>& omegabet, const double& AreaW, const double& integral);
RcppExport SEXP _binspp_logpXCbetC(SEXP YSEXP, SEXP CCSEXP, SEXP z_alphaSEXP, SEXP z_omegaSEXP, SEXP alphaSEXP, SEXP alphabetSEXP, SEXP omegaSEXP, SEXP omegabetSEXP, SEXP AreaWSEXP, SEXP integralSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_alpha(z_alphaSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_omega(z_omegaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type omegabet(omegabetSEXP);
    Rcpp::traits::input_parameter< const double& >::type AreaW(AreaWSEXP);
    Rcpp::traits::input_parameter< const double& >::type integral(integralSEXP);
    rcpp_result_gen = Rcpp::wrap(logpXCbetC(Y, CC, z_alpha, z_omega, alpha, alphabet, omega, omegabet, AreaW, integral));
    return rcpp_result_gen;
END_RCPP
}
// KumulaVsechC
double KumulaVsechC(const NumericMatrix& CC, const List& z_alpha, const List& z_omega, const double& alpha, const Nullable<NumericVector>& alphabet, const double& omega, const Nullable<NumericVector>& omegabet, const NumericVector& x_left, const NumericVector& x_right, const NumericVector& y_bottom, const NumericVector& y_top);
RcppExport SEXP _binspp_KumulaVsechC(SEXP CCSEXP, SEXP z_alphaSEXP, SEXP z_omegaSEXP, SEXP alphaSEXP, SEXP alphabetSEXP, SEXP omegaSEXP, SEXP omegabetSEXP, SEXP x_leftSEXP, SEXP x_rightSEXP, SEXP y_bottomSEXP, SEXP y_topSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_alpha(z_alphaSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_omega(z_omegaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type omegabet(omegabetSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_left(x_leftSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_right(x_rightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_bottom(y_bottomSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_top(y_topSEXP);
    rcpp_result_gen = Rcpp::wrap(KumulaVsechC(CC, z_alpha, z_omega, alpha, alphabet, omega, omegabet, x_left, x_right, y_bottom, y_top));
    return rcpp_result_gen;
END_RCPP
}
// PrioralphaC
double PrioralphaC(const double& a, const double& Prior_alpha_mean, const double& Prior_alpha_SD);
RcppExport SEXP _binspp_PrioralphaC(SEXP aSEXP, SEXP Prior_alpha_meanSEXP, SEXP Prior_alpha_SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_alpha_mean(Prior_alpha_meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_alpha_SD(Prior_alpha_SDSEXP);
    rcpp_result_gen = Rcpp::wrap(PrioralphaC(a, Prior_alpha_mean, Prior_alpha_SD));
    return rcpp_result_gen;
END_RCPP
}
// PrioromegaC
double PrioromegaC(const double& o, const double& Prior_omega_mean, const double& Prior_omega_SD);
RcppExport SEXP _binspp_PrioromegaC(SEXP oSEXP, SEXP Prior_omega_meanSEXP, SEXP Prior_omega_SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type o(oSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_omega_mean(Prior_omega_meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_omega_SD(Prior_omega_SDSEXP);
    rcpp_result_gen = Rcpp::wrap(PrioromegaC(o, Prior_omega_mean, Prior_omega_SD));
    return rcpp_result_gen;
END_RCPP
}
// PrioralphabetC
double PrioralphabetC(const Nullable<NumericVector>& a, const Nullable<NumericVector>& Prior_alphavec_SD);
RcppExport SEXP _binspp_PrioralphabetC(SEXP aSEXP, SEXP Prior_alphavec_SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type Prior_alphavec_SD(Prior_alphavec_SDSEXP);
    rcpp_result_gen = Rcpp::wrap(PrioralphabetC(a, Prior_alphavec_SD));
    return rcpp_result_gen;
END_RCPP
}
// PrioromegabetC
double PrioromegabetC(const Nullable<NumericVector>& o, const Nullable<NumericVector>& Prior_omegavec_SD);
RcppExport SEXP _binspp_PrioromegabetC(SEXP oSEXP, SEXP Prior_omegavec_SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type o(oSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type Prior_omegavec_SD(Prior_omegavec_SDSEXP);
    rcpp_result_gen = Rcpp::wrap(PrioromegabetC(o, Prior_omegavec_SD));
    return rcpp_result_gen;
END_RCPP
}
// StepbetC
List StepbetC(const double& kappa, const List& z_alpha, const List& z_omega, const double& alpha, const double& salpha, const Nullable<NumericVector>& alphabet, const Nullable<NumericVector>& salphabet, const double& omega, const double& somega, const Nullable<NumericVector>& omegabet, const Nullable<NumericVector>& somegabet, const NumericMatrix& Y, const NumericMatrix& CC, const double& logP, const double& integral, const double& integralrho, const NumericVector& x_left, const NumericVector& x_right, const NumericVector& y_bottom, const NumericVector& y_top, const List& Wpix, const double& AreaW, const double& FScoef, const double& Prior_alpha_mean, const double& Prior_alpha_SD, const double& Prior_omega_mean, const double& Prior_omega_SD, const Nullable<NumericVector>& Prior_alphavec_SD, const Nullable<NumericVector>& Prior_omegavec_SD);
RcppExport SEXP _binspp_StepbetC(SEXP kappaSEXP, SEXP z_alphaSEXP, SEXP z_omegaSEXP, SEXP alphaSEXP, SEXP salphaSEXP, SEXP alphabetSEXP, SEXP salphabetSEXP, SEXP omegaSEXP, SEXP somegaSEXP, SEXP omegabetSEXP, SEXP somegabetSEXP, SEXP YSEXP, SEXP CCSEXP, SEXP logPSEXP, SEXP integralSEXP, SEXP integralrhoSEXP, SEXP x_leftSEXP, SEXP x_rightSEXP, SEXP y_bottomSEXP, SEXP y_topSEXP, SEXP WpixSEXP, SEXP AreaWSEXP, SEXP FScoefSEXP, SEXP Prior_alpha_meanSEXP, SEXP Prior_alpha_SDSEXP, SEXP Prior_omega_meanSEXP, SEXP Prior_omega_SDSEXP, SEXP Prior_alphavec_SDSEXP, SEXP Prior_omegavec_SDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_alpha(z_alphaSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_omega(z_omegaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type salpha(salphaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type salphabet(salphabetSEXP);
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const double& >::type somega(somegaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type omegabet(omegabetSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type somegabet(somegabetSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< const double& >::type logP(logPSEXP);
    Rcpp::traits::input_parameter< const double& >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const double& >::type integralrho(integralrhoSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_left(x_leftSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_right(x_rightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_bottom(y_bottomSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_top(y_topSEXP);
    Rcpp::traits::input_parameter< const List& >::type Wpix(WpixSEXP);
    Rcpp::traits::input_parameter< const double& >::type AreaW(AreaWSEXP);
    Rcpp::traits::input_parameter< const double& >::type FScoef(FScoefSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_alpha_mean(Prior_alpha_meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_alpha_SD(Prior_alpha_SDSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_omega_mean(Prior_omega_meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type Prior_omega_SD(Prior_omega_SDSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type Prior_alphavec_SD(Prior_alphavec_SDSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type Prior_omegavec_SD(Prior_omegavec_SDSEXP);
    rcpp_result_gen = Rcpp::wrap(StepbetC(kappa, z_alpha, z_omega, alpha, salpha, alphabet, salphabet, omega, somega, omegabet, somegabet, Y, CC, logP, integral, integralrho, x_left, x_right, y_bottom, y_top, Wpix, AreaW, FScoef, Prior_alpha_mean, Prior_alpha_SD, Prior_omega_mean, Prior_omega_SD, Prior_alphavec_SD, Prior_omegavec_SD));
    return rcpp_result_gen;
END_RCPP
}
// row_add
NumericMatrix row_add(const NumericMatrix& x, const NumericVector& extraRow);
RcppExport SEXP _binspp_row_add(SEXP xSEXP, SEXP extraRowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type extraRow(extraRowSEXP);
    rcpp_result_gen = Rcpp::wrap(row_add(x, extraRow));
    return rcpp_result_gen;
END_RCPP
}
// rand_int
int rand_int(const int& min, const int& max);
RcppExport SEXP _binspp_rand_int(SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const int& >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(rand_int(min, max));
    return rcpp_result_gen;
END_RCPP
}
// StepMovePointC
List StepMovePointC(const double& kappa, const List& z_alpha, const List& z_omega, const double& alpha, const Nullable<NumericVector>& alphabet, const double& omega, const Nullable<NumericVector>& omegabet, const NumericMatrix& Y, const NumericMatrix& CC, const double& logP, const double& integral, const double& integralrho, const NumericVector& x_left, const NumericVector& x_right, const NumericVector& y_bottom, const NumericVector& y_top, const double& AreaW, const NumericVector& NewCenter);
RcppExport SEXP _binspp_StepMovePointC(SEXP kappaSEXP, SEXP z_alphaSEXP, SEXP z_omegaSEXP, SEXP alphaSEXP, SEXP alphabetSEXP, SEXP omegaSEXP, SEXP omegabetSEXP, SEXP YSEXP, SEXP CCSEXP, SEXP logPSEXP, SEXP integralSEXP, SEXP integralrhoSEXP, SEXP x_leftSEXP, SEXP x_rightSEXP, SEXP y_bottomSEXP, SEXP y_topSEXP, SEXP AreaWSEXP, SEXP NewCenterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_alpha(z_alphaSEXP);
    Rcpp::traits::input_parameter< const List& >::type z_omega(z_omegaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type omegabet(omegabetSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< const double& >::type logP(logPSEXP);
    Rcpp::traits::input_parameter< const double& >::type integral(integralSEXP);
    Rcpp::traits::input_parameter< const double& >::type integralrho(integralrhoSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_left(x_leftSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x_right(x_rightSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_bottom(y_bottomSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y_top(y_topSEXP);
    Rcpp::traits::input_parameter< const double& >::type AreaW(AreaWSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type NewCenter(NewCenterSEXP);
    rcpp_result_gen = Rcpp::wrap(StepMovePointC(kappa, z_alpha, z_omega, alpha, alphabet, omega, omegabet, Y, CC, logP, integral, integralrho, x_left, x_right, y_bottom, y_top, AreaW, NewCenter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_binspp_addMat", (DL_FUNC) &_binspp_addMat, 2},
    {"_binspp_multMat", (DL_FUNC) &_binspp_multMat, 2},
    {"_binspp_intalphaC", (DL_FUNC) &_binspp_intalphaC, 4},
    {"_binspp_aozC", (DL_FUNC) &_binspp_aozC, 4},
    {"_binspp_logpXCbetC", (DL_FUNC) &_binspp_logpXCbetC, 10},
    {"_binspp_KumulaVsechC", (DL_FUNC) &_binspp_KumulaVsechC, 11},
    {"_binspp_PrioralphaC", (DL_FUNC) &_binspp_PrioralphaC, 3},
    {"_binspp_PrioromegaC", (DL_FUNC) &_binspp_PrioromegaC, 3},
    {"_binspp_PrioralphabetC", (DL_FUNC) &_binspp_PrioralphabetC, 2},
    {"_binspp_PrioromegabetC", (DL_FUNC) &_binspp_PrioromegabetC, 2},
    {"_binspp_StepbetC", (DL_FUNC) &_binspp_StepbetC, 29},
    {"_binspp_row_add", (DL_FUNC) &_binspp_row_add, 2},
    {"_binspp_rand_int", (DL_FUNC) &_binspp_rand_int, 2},
    {"_binspp_StepMovePointC", (DL_FUNC) &_binspp_StepMovePointC, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_binspp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}