//
//  functions used in spBFA package
//

#ifndef __spBFA__
#define __spBFA__

//MCMC Sampler
Rcpp::List bfa_sp_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive);
  
//STRUCT DEFINITIONS
struct datobj {
  double Rho;
  double ScaleY;
  int N;
  int M;
  int Nu;
  int K;
  int L;
  int FamilyInd;
  int TempCorInd;
  arma::colvec YStar;
  arma::colvec YObserved;
  arma::mat YStarWide;
  arma::Mat<int> W;
  arma::colvec Time;
  arma::mat ICAR;
  arma::mat ICARInv;
  arma::mat TimeDist;
  arma::mat EyeNu;
  arma::Col<int> SeqL;
  arma::mat EyeM;
  arma::mat EyeKbyNu;
  arma::colvec ZeroKbyNu;
  arma::colvec OneNu;
};
struct hypara {
  double A;
  double B;
  double C;
  double D;
  double A1;
  double A2;
  double APsi;
  double BPsi;
  double Gamma;
  double Beta;
  double Zeta;
  arma::mat Omega;
};
struct metrobj {
  double MetropPsi;
  int AcceptancePsi;
  double OriginalTuners;
};
struct para {
  arma::colvec Sigma2;
  double Kappa2;
  arma::colvec Delta;
  double Psi;
  arma::mat Upsilon;
  arma::mat UpsilonInv;
  arma::umat Xi;
  arma::mat Theta;
  arma::mat Lambda;
  arma::colvec Tau;
  arma::mat BigPhi;
  arma::colvec Eta;
  arma::cube Alpha;
  arma::cube Z;
  arma::mat BigPsi;
  arma::mat Sigma;
  arma::mat SigmaInv;
  arma::mat HPsi;
  arma::mat CholHPsi;
  arma::mat HPsiInv;
  arma::colvec Mean;
  arma::cube Weights;
  arma::cube logWeights;
};
struct dataug {
  int NBelow;
  int NAbove;
  arma::uvec WhichAbove;
  arma::uvec WhichBelow;
};
struct mcmcobj {
  int NBurn;
  int NSims;
  int NThin;
  int NPilot;
  int NTotal;
  int NKeep;
  arma::vec WhichKeep;
  arma::vec WhichPilotAdapt;
  arma::vec WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress;
  arma::vec BurnInProgress;
  int BarLength;
  int PilotAdaptDenominator;
};

//COVARIANCE FUNCTIONS
arma::mat H(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu);
  
//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
double pnormRcpp(double q);
double lpnormRcpp(double q);
double UpperpnormRcpp(double q);
double lUpperpnormRcpp(double q);
double rigammaRcpp(double Alpha, double Theta);
double rgammaRcpp(double Alpha, double Theta);
arma::mat rwishRcpp(double n, arma::mat const& V);
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
double rtnormRcpp(double mean, double sd, bool Above);
arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper);

//MCMC CONVERSION FUNCTIONS
datobj ConvertDatObj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
metrobj ConvertMetrObj(Rcpp::List MetrObj_List);
para ConvertPara(Rcpp::List Para_List);
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List);
dataug ConvertDatAug(Rcpp::List DatAug_List);

//MCMC SAMPLER FUNCTIONS
para SampleTheta(datobj DatObj, para Para);
para SampleXi(datobj DatObj, para Para);
para SampleZ(datobj DatObj, para Para);
para SampleAlpha(datobj DatObj, para Para);
para SampleKappa2(datobj DatObj, para Para, hypara HyPara);
para SampleDelta(datobj DatObj, para Para, hypara HyPara);
para SampleEta(datobj DatObj, para Para, hypara HyPara);
para SampleUpsilon(datobj DatObj, para Para, hypara HyPara);
std::pair<para, metrobj> SamplePsi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
para SampleSigma2(datobj DatObj, para Para, hypara HyPara);
arma::colvec SampleProbit(datobj DatObj, para Para, dataug DatAug);
arma::colvec SampleTobit(datobj DatObj, para Para, dataug DatAug);
datobj SampleY(datobj DatObj, para Para, dataug DatAug);

//MCMC UTILITY FUNCTIONS
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive);
Rcpp::List OutputMetrObj(metrobj MetrObj);
metrobj PilotAdaptation(metrobj MetrObj, mcmcobj McmcObj);
void SamplerProgress(int s, mcmcobj McmcObj);
arma::colvec StoreSamples(datobj DatObj, para Para);
void UpdateBurnInBar(int s, mcmcobj McmcObj);
void UpdateBurnInBarInt(int s, mcmcobj McmcObj);

//UTILITY FUNCTIONS
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye);
arma::mat CholInv(arma::mat const& Cov);
arma::mat GetLambda(arma::mat const& Theta, arma::umat const& Xi, int K, int M);
arma::cube GetlogWeights(arma::cube const& Alpha, int K, int M, int L);
arma::cube GetWeights(arma::cube const& Alpha, int K, int M, int L);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);

#endif // __spBFA__
