#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive) {

  //Set MCMC object
  int BarLength = McmcObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  0%..  ";
  }

}



//Function to pilot adapt tuning parameter--------------------------------------------------------------
double PilotAdaptFunc(double TuningParameter, double AcceptancePct) {

  //Adjust tuning parameter using scaling based on size of acceptance rate
  if (AcceptancePct >= 0.90) TuningParameter *= 1.3;
  if ( (AcceptancePct >= 0.75 ) & (AcceptancePct < 0.90 ) ) TuningParameter *= 1.2;
  if ( (AcceptancePct >= 0.45 ) & (AcceptancePct < 0.75 ) ) TuningParameter *= 1.1;
  if ( (AcceptancePct <= 0.25 ) & (AcceptancePct > 0.15 ) ) TuningParameter *= 0.9;
  if ( (AcceptancePct <= 0.15 ) & (AcceptancePct > 0.10 ) ) TuningParameter *= 0.8;
  if (AcceptancePct <= 0.10) TuningParameter *= 0.7;
  return TuningParameter;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj PilotAdaptation(metrobj MetrObj, mcmcobj McmcObj) {

  //Set Metropolis objects
  double MetropPsi = MetrObj.MetropPsi;
  double AcceptancePsi = MetrObj.AcceptancePsi;

  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Get acceptance percentages
  double PctPsi = AcceptancePsi / double(PilotAdaptDenominator);

  //Update Tuning Parameter
  MetropPsi = PilotAdaptFunc(MetropPsi, PctPsi);
  MetrObj.MetropPsi = MetropPsi;

  //Zero the acceptance counters
  AcceptancePsi = 0;
  MetrObj.AcceptancePsi = AcceptancePsi;
  return MetrObj;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(metrobj MetrObj) {

  Rcpp::List Out = Rcpp::List::create(Rcpp::Named("AcceptancePsi") = MetrObj.AcceptancePsi,
                                      Rcpp::Named("MetropPsi") = MetrObj.MetropPsi);
  return Out;

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int s, mcmcobj McmcObj) {

  //Set MCMC object
  int NSims = McmcObj.NSims;
  int NBurn = McmcObj.NBurn;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%..  ";

}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples(datobj DatObj, para Para) {

  //Set data object
  int M = DatObj.M;
  int K = DatObj.K;
  int Nu = DatObj.Nu;

  //Set parameter objects
  arma::mat Lambda = Para.Lambda;
  arma::mat BigPhi = Para.BigPhi;
  arma::colvec Sigma2 = Para.Sigma2;
  double Kappa2 = Para.Kappa2;
  arma::colvec Delta = Para.Delta;
  arma::mat Upsilon = Para.Upsilon;
  double Psi = Para.Psi;

  //Save raw samples
  int counter = 0;
  arma::colvec col(M * K + K * Nu + M + 1 + K + ((K + 1) * K) / 2 + 1);
  for (arma::uword i = 0; i < M; i++) {
    for (arma::uword j = 0; j < K; j++) {
      col(counter) = Lambda(i, j);
      counter++;
    }
  }
  for (arma::uword j = 0; j < K; j++) {
    for (arma::uword t = 0; t < Nu; t++) {
      col(counter) = BigPhi(j, t);
      counter++;
    }
  }
  for (arma::uword i = 0; i < M; i++) {
    col(counter) = Sigma2(i);
    counter++;
  }
  col(counter) = Kappa2;
  counter++;
  for (arma::uword j = 0; j < K; j++) {
    col(counter) = Delta(j);
    counter++;
  }
  for (arma::uword i = 0; i < K; i++) {
    for (arma::uword j = 0; j <= i; j++) {
      col(counter) = Upsilon(i, j);
      counter++;
    }
  }
  col(counter) = Psi;
  return col;
}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Add a new star
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }




//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBarInt(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgressInt);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);

  //Add percentage to submited job mode
  Rcpp::Rcout.precision(0);
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%..  ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%..  ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%..  ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%..  ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%..  ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%..  ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%..  ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%..  ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%..  ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!  ";

}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Number of new star in interactive mode
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBar(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  int BarLength = McmcObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rBurn-in progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";

}

