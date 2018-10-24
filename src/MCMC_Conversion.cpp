#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  double ScaleY = DatObj_List["ScaleY"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::colvec YStar = DatObj_List["YStar"];
  arma::colvec YObserved = DatObj_List["YObserved"];
  arma::Mat<int> W = DatObj_List["W"];
  arma::mat TimeDist = DatObj_List["TimeDist"];
  arma::colvec Time = DatObj_List["Time"];
  double Rho = DatObj_List["Rho"];
  int N = DatObj_List["N"];
  int K = DatObj_List["K"];
  int L = DatObj_List["L"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int TempCorInd = DatObj_List["TempCorInd"];
  int FamilyInd = DatObj_List["FamilyInd"];
  arma::mat ICAR = DatObj_List["ICAR"];
  arma::mat ICARInv = DatObj_List["ICARInv"];
  arma::mat EyeNu = DatObj_List["EyeNu"];
  arma::Col<int> SeqL = DatObj_List["SeqL"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::mat EyeKbyNu = DatObj_List["EyeKbyNu"];
  arma::colvec ZeroKbyNu = DatObj_List["ZeroKbyNu"];
  arma::colvec OneNu = DatObj_List["OneNu"];

  //Convert to C++ struct
  datobj DatObj;
  DatObj.ScaleY = ScaleY;
  DatObj.YStarWide = YStarWide;
  DatObj.YStar = YStar;
  DatObj.YObserved = YObserved;
  DatObj.W = W;
  DatObj.TimeDist = TimeDist;
  DatObj.Time = Time;
  DatObj.Rho = Rho;
  DatObj.N = N;
  DatObj.K = K;
  DatObj.L = L;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.TempCorInd = TempCorInd;
  DatObj.FamilyInd = FamilyInd;
  DatObj.ICAR = ICAR;
  DatObj.ICARInv = ICARInv;
  DatObj.EyeNu = EyeNu;
  DatObj.SeqL = SeqL;
  DatObj.EyeM = EyeM;
  DatObj.EyeKbyNu = EyeKbyNu;
  DatObj.ZeroKbyNu = ZeroKbyNu;
  DatObj.OneNu = OneNu;
  return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  double A = HyPara_List["A"];
  double B = HyPara_List["B"];
  double C = HyPara_List["C"];
  double D = HyPara_List["D"];
  double A1 = HyPara_List["A1"];
  double A2 = HyPara_List["A2"];
  double APsi = HyPara_List["APsi"];
  double BPsi = HyPara_List["BPsi"];
  double Gamma = HyPara_List["Gamma"];
  double Beta = HyPara_List["Beta"];
  double Zeta = HyPara_List["Zeta"];
  arma::mat Omega = HyPara_List["Omega"];
  
  //Convert to C++ struct
  hypara HyPara;
  HyPara.A = A;
  HyPara.B = B;
  HyPara.C = C;
  HyPara.D = D;
  HyPara.A1 = A1;
  HyPara.A2 = A2;
  HyPara.APsi = APsi;
  HyPara.BPsi = BPsi;
  HyPara.Gamma = Gamma;
  HyPara.Beta = Beta;
  HyPara.Zeta = Zeta;
  HyPara.Omega = Omega;
  return HyPara;

}



//Function to convert Rcpp::List MetrObj to a custom C++ struct metrobj-----------------------------------------------
metrobj ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  double MetropPsi = MetrObj_List["MetropPsi"];
  int AcceptancePsi = MetrObj_List["AcceptancePsi"];
  double OriginalTuners = MetrObj_List["OriginalTuners"];
  
  //Convert to C++ struct
  metrobj MetrObj;
  MetrObj.MetropPsi = MetropPsi;
  MetrObj.AcceptancePsi = AcceptancePsi;
  MetrObj.OriginalTuners = OriginalTuners;
  return MetrObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {

  //Set objects from List
  arma::colvec Sigma2 = Para_List["Sigma2"];
  double Kappa2 = Para_List["Kappa2"];
  arma::colvec Delta = Para_List["Delta"];
  double Psi = Para_List["Psi"];
  arma::mat Upsilon = Para_List["Upsilon"];
  arma::mat UpsilonInv = Para_List["UpsilonInv"];
  arma::umat Xi = Para_List["Xi"];
  arma::mat Theta = Para_List["Theta"];
  arma::mat Lambda = Para_List["Lambda"];
  arma::colvec Tau = Para_List["Tau"];
  arma::mat BigPhi = Para_List["BigPhi"];
  arma::colvec Eta = Para_List["Eta"];
  arma::cube Alpha = Para_List["Alpha"];
  arma::cube Z = Para_List["Z"];
  arma::mat BigPsi = Para_List["BigPsi"];
  arma::mat Sigma = Para_List["Sigma"];
  arma::mat SigmaInv = Para_List["SigmaInv"];
  arma::mat HPsi = Para_List["HPsi"];
  arma::mat CholHPsi = Para_List["CholHPsi"];
  arma::mat HPsiInv = Para_List["HPsiInv"];
  arma::colvec Mean = Para_List["Mean"];
  arma::cube Weights = Para_List["Weights"];
  arma::cube logWeights = Para_List["logWeights"];

  //Convert to C++ struct
  para Para;
  Para.Sigma2 = Sigma2;
  Para.Kappa2 = Kappa2;
  Para.Delta = Delta;
  Para.Psi = Psi;
  Para.Upsilon = Upsilon;
  Para.UpsilonInv = UpsilonInv;
  Para.Xi = Xi;
  Para.Theta = Theta;
  Para.Lambda = Lambda;
  Para.Tau = Tau;
  Para.BigPhi = BigPhi;
  Para.Eta = Eta;
  Para.Alpha = Alpha;
  Para.Z = Z;
  Para.BigPsi = BigPsi;
  Para.Sigma = Sigma;
  Para.SigmaInv = SigmaInv;
  Para.HPsi = HPsi;
  Para.CholHPsi = CholHPsi;
  Para.HPsiInv = HPsiInv;
  Para.Mean = Mean;
  Para.Weights = Weights;
  Para.logWeights = logWeights;
  return Para;
}



//Function to convert Rcpp::List DatAug to a custom C++ struct dataug-----------------------------------------------------
dataug ConvertDatAug(Rcpp::List DatAug_List) {

  //Set objects from List
  int NBelow = DatAug_List["NBelow"];
  int NAbove = DatAug_List["NAbove"];
  arma::uvec WhichBelow = DatAug_List["WhichBelow"];
  arma::uvec WhichAbove = DatAug_List["WhichAbove"];

  //Convert to C++ struct
  dataug DatAug;
  DatAug.NBelow = NBelow;
  DatAug.NAbove = NAbove;
  DatAug.WhichBelow = WhichBelow;
  DatAug.WhichAbove = WhichAbove;
  return DatAug;
}



//Function to convert Rcpp::List McmcObj to a custom C++ struct mcmcmobj-----------------------------------------------------
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List) {

  //Set objects from List
  int NBurn = McmcObj_List["NBurn"];
  int NSims = McmcObj_List["NSims"];
  int NThin = McmcObj_List["NThin"];
  int NPilot = McmcObj_List["NPilot"];
  int NTotal = McmcObj_List["NTotal"];
  int NKeep = McmcObj_List["NKeep"];
  arma::vec WhichKeep = McmcObj_List["WhichKeep"];
  arma::vec WhichPilotAdapt = McmcObj_List["WhichPilotAdapt"];
  arma::vec WhichBurnInProgress = McmcObj_List["WhichBurnInProgress"];
  arma::vec WhichBurnInProgressInt = McmcObj_List["WhichBurnInProgressInt"];
  arma::vec WhichSamplerProgress = McmcObj_List["WhichSamplerProgress"];
  arma::vec BurnInProgress = McmcObj_List["BurnInProgress"];
  int BarLength = McmcObj_List["BarLength"];
  int PilotAdaptDenominator = McmcObj_List["PilotAdaptDenominator"];

  //Convert to C++ struct
  mcmcobj McmcObj;
  McmcObj.NBurn = NBurn;
  McmcObj.NSims = NSims;
  McmcObj.NThin = NThin;
  McmcObj.NPilot = NPilot;
  McmcObj.NTotal = NTotal;
  McmcObj.NKeep = NKeep;
  McmcObj.WhichKeep = WhichKeep;
  McmcObj.WhichPilotAdapt = WhichPilotAdapt;
  McmcObj.WhichBurnInProgress = WhichBurnInProgress;
  McmcObj.WhichBurnInProgressInt = WhichBurnInProgressInt;
  McmcObj.WhichSamplerProgress = WhichSamplerProgress;
  McmcObj.BurnInProgress = BurnInProgress;
  McmcObj.PilotAdaptDenominator = PilotAdaptDenominator;
  McmcObj.BarLength = BarLength;
  return McmcObj;
}



