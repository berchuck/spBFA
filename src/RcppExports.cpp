// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GetRooti
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye);
RcppExport SEXP _spBFA_GetRooti(SEXP CovSEXP, SEXP EyeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Cov(CovSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Eye(EyeSEXP);
    rcpp_result_gen = Rcpp::wrap(GetRooti(Cov, Eye));
    return rcpp_result_gen;
END_RCPP
}
// H
arma::mat H(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu);
RcppExport SEXP _spBFA_H(SEXP PsiSEXP, SEXP TempCorIndSEXP, SEXP TimeDistSEXP, SEXP NuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< int >::type TempCorInd(TempCorIndSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type TimeDist(TimeDistSEXP);
    Rcpp::traits::input_parameter< int >::type Nu(NuSEXP);
    rcpp_result_gen = Rcpp::wrap(H(Psi, TempCorInd, TimeDist, Nu));
    return rcpp_result_gen;
END_RCPP
}
// SpEXP
arma::mat SpEXP(double Rho, arma::mat const& SpDist, int M);
RcppExport SEXP _spBFA_SpEXP(SEXP RhoSEXP, SEXP SpDistSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type SpDist(SpDistSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(SpEXP(Rho, SpDist, M));
    return rcpp_result_gen;
END_RCPP
}
// GetLogLik
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spBFA_GetLogLik(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogLik(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// GetLogLikMean
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List);
RcppExport SEXP _spBFA_GetLogLikMean(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogLikMean(DatObj_List, Para_List));
    return rcpp_result_gen;
END_RCPP
}
// SamplePPD
arma::mat SamplePPD(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spBFA_SamplePPD(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(SamplePPD(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// bfa_sp_Rcpp
Rcpp::List bfa_sp_Rcpp(Rcpp::List DatObj_List, Rcpp::List HyPara_List, Rcpp::List MetrObj_List, Rcpp::List Para_List, Rcpp::List DatAug_List, Rcpp::List McmcObj_List, arma::mat RawSamples, bool Interactive);
RcppExport SEXP _spBFA_bfa_sp_Rcpp(SEXP DatObj_ListSEXP, SEXP HyPara_ListSEXP, SEXP MetrObj_ListSEXP, SEXP Para_ListSEXP, SEXP DatAug_ListSEXP, SEXP McmcObj_ListSEXP, SEXP RawSamplesSEXP, SEXP InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type HyPara_List(HyPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MetrObj_List(MetrObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type DatAug_List(DatAug_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type McmcObj_List(McmcObj_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RawSamples(RawSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type Interactive(InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(bfa_sp_Rcpp(DatObj_List, HyPara_List, MetrObj_List, Para_List, DatAug_List, McmcObj_List, RawSamples, Interactive));
    return rcpp_result_gen;
END_RCPP
}
// EtaKrigging
arma::mat EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spBFA_EtaKrigging(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(EtaKrigging(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// YKrigging
arma::cube YKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose);
RcppExport SEXP _spBFA_YKrigging(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP EtaKrigSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type EtaKrig(EtaKrigSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(YKrigging(DatObj_List, Para_List, EtaKrig, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// Play
void Play();
RcppExport SEXP _spBFA_Play() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Play();
    return R_NilValue;
END_RCPP
}
// GetLStarJ
arma::colvec GetLStarJ(arma::mat const& U, arma::cube const& Weights, int K, int M, int O);
RcppExport SEXP _spBFA_GetLStarJ(SEXP USEXP, SEXP WeightsSEXP, SEXP KSEXP, SEXP MSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::cube const& >::type Weights(WeightsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLStarJ(U, Weights, K, M, O));
    return rcpp_result_gen;
END_RCPP
}
// GetLambda
arma::mat GetLambda(arma::mat const& Theta, arma::umat const& Xi, int K, int M, int O);
RcppExport SEXP _spBFA_GetLambda(SEXP ThetaSEXP, SEXP XiSEXP, SEXP KSEXP, SEXP MSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< arma::umat const& >::type Xi(XiSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLambda(Theta, Xi, K, M, O));
    return rcpp_result_gen;
END_RCPP
}
// GetWeights
arma::cube GetWeights(arma::cube const& Alpha, int K, int M, int L, int O);
RcppExport SEXP _spBFA_GetWeights(SEXP AlphaSEXP, SEXP KSEXP, SEXP MSEXP, SEXP LSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(GetWeights(Alpha, K, M, L, O));
    return rcpp_result_gen;
END_RCPP
}
// GetlogWeights
arma::cube GetlogWeights(arma::cube const& Alpha, int K, int M, int L, int O);
RcppExport SEXP _spBFA_GetlogWeights(SEXP AlphaSEXP, SEXP KSEXP, SEXP MSEXP, SEXP LSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(GetlogWeights(Alpha, K, M, L, O));
    return rcpp_result_gen;
END_RCPP
}
// CholInv
arma::mat CholInv(arma::mat const& Cov);
RcppExport SEXP _spBFA_CholInv(SEXP CovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Cov(CovSEXP);
    rcpp_result_gen = Rcpp::wrap(CholInv(Cov));
    return rcpp_result_gen;
END_RCPP
}
// Inv3
arma::mat Inv3(arma::mat const& A);
RcppExport SEXP _spBFA_Inv3(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(Inv3(A));
    return rcpp_result_gen;
END_RCPP
}
// makeSymm
arma::mat makeSymm(arma::mat const& A);
RcppExport SEXP _spBFA_makeSymm(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(makeSymm(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spBFA_GetRooti", (DL_FUNC) &_spBFA_GetRooti, 2},
    {"_spBFA_H", (DL_FUNC) &_spBFA_H, 4},
    {"_spBFA_SpEXP", (DL_FUNC) &_spBFA_SpEXP, 3},
    {"_spBFA_GetLogLik", (DL_FUNC) &_spBFA_GetLogLik, 4},
    {"_spBFA_GetLogLikMean", (DL_FUNC) &_spBFA_GetLogLikMean, 2},
    {"_spBFA_SamplePPD", (DL_FUNC) &_spBFA_SamplePPD, 4},
    {"_spBFA_bfa_sp_Rcpp", (DL_FUNC) &_spBFA_bfa_sp_Rcpp, 8},
    {"_spBFA_EtaKrigging", (DL_FUNC) &_spBFA_EtaKrigging, 4},
    {"_spBFA_YKrigging", (DL_FUNC) &_spBFA_YKrigging, 5},
    {"_spBFA_Play", (DL_FUNC) &_spBFA_Play, 0},
    {"_spBFA_GetLStarJ", (DL_FUNC) &_spBFA_GetLStarJ, 5},
    {"_spBFA_GetLambda", (DL_FUNC) &_spBFA_GetLambda, 5},
    {"_spBFA_GetWeights", (DL_FUNC) &_spBFA_GetWeights, 5},
    {"_spBFA_GetlogWeights", (DL_FUNC) &_spBFA_GetlogWeights, 5},
    {"_spBFA_CholInv", (DL_FUNC) &_spBFA_CholInv, 1},
    {"_spBFA_Inv3", (DL_FUNC) &_spBFA_Inv3, 1},
    {"_spBFA_makeSymm", (DL_FUNC) &_spBFA_makeSymm, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_spBFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
