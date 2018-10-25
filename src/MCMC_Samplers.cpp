#define ARMA_DONT_PRINT_ERRORS //So the cholesky warning is suppressed
#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//Function to sample latent probit process using Gibbs sampling step------------------------------------------------
arma::colvec SampleProbit(datobj DatObj, para Para, dataug DatAug) {
  
  //Set data objects
  arma::colvec YStar = DatObj.YStar;
  arma::colvec OneNu = DatObj.OneNu;
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  arma::colvec Sigma2 = Para.Sigma2;
  
  //Set data augmentation objects
  arma::uvec WhichBelow = DatAug.WhichBelow;
  arma::uvec WhichAbove = DatAug.WhichAbove;
  int NBelow = DatAug.NBelow;
  int NAbove = DatAug.NAbove;
  
  //Moments
  arma::colvec Mu = Mean(WhichBelow);
  arma::colvec SDFull = arma::kron(OneNu, arma::sqrt(Sigma2));
  arma::colvec SD = SDFull(WhichBelow);
  
  //Sample latent Variable from full conditional
  for (int i = 0; i < NBelow; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), -100000, 0);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Probit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Probit model. Consider scaling towards zero and re-running.");
    YStar(WhichBelow(i)) = Temp;
  }
  for (int i = 0; i < NAbove; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), false);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), 0, 100000);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Probit sampling step. Most likey cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Probit model. Consider scaling towards zero and re-running.");
    YStar(WhichAbove(i)) = Temp;
  }
  return YStar;
  
}



//Function to sample latent tobit process using Gibbs sampling step------------------------------------------------
arma::colvec SampleTobit(datobj DatObj, para Para, dataug DatAug) {
  
  //Set data objects
  arma::colvec YStar = DatObj.YStar;
  arma::colvec OneNu = DatObj.OneNu;
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  arma::colvec Sigma2 = Para.Sigma2;
  
  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::uvec WhichBelow = DatAug.WhichBelow;
  
  //Moments
  arma::colvec Mu = Mean(WhichBelow);
  arma::colvec SDFull = arma::kron(OneNu, arma::sqrt(Sigma2));
  arma::colvec SD = SDFull(WhichBelow);
  
  //Sample latent Variable from full conditional
  for (int i = 0; i < NBelow; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), true);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), -arma::datum::inf, 0);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit sampling step. Most likely cause for this error is that the data being used is inappropriate (i.e., to far from zero) for a Tobit model. Consider scaling towards zero and re-running.");
    YStar(WhichBelow(i)) = Temp;
  }
  return YStar;
  
}



//Function to sample latent process from its full conditional------------------------------------------------------
datobj SampleY(datobj DatObj, para Para, dataug DatAug) {
  
  //Set data objects
  int FamilyInd = DatObj.FamilyInd;
  int N = DatObj.N;
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  
  //Sample latent process
  arma::vec YStar(N);
  if (FamilyInd == 1) YStar = SampleProbit(DatObj, Para, DatAug);
  if (FamilyInd == 2) YStar = SampleTobit(DatObj, Para, DatAug);
  
  //Save output
  DatObj.YStar = YStar;
  DatObj.YStarWide = arma::reshape(YStar, M, Nu);
  return DatObj;
  
}



//Function to sample sigma2(s_i) using a Gibbs sampler step---------------------------------------------------------------
para SampleSigma2(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  arma::mat YStarWide = DatObj.YStarWide;
  
  //Set parameter objects
  arma::colvec Mean = Para.Mean;
  arma::mat MeanMat = arma::reshape(Mean, M, Nu);
  arma::colvec Sigma2 = Para.Sigma2;
  arma::mat Lambda = Para.Lambda;
  
  //Set hyperparameter objects
  double A = HyPara.A;
  double B = HyPara.B;
  
  //Shape is constant over locations
  double Shape = A + 0.5 * Nu;
  
  //Loop over locations
  for (arma::uword i = 0; i < M; i++) {
    
    //Location specific objects
    arma::rowvec Diff = YStarWide.row(i) - MeanMat.row(i);

    //Calculate rate
    double Resids = arma::as_scalar(Diff * arma::trans(Diff));
    double Rate = B + 0.5 * Resids;
  
    //Sample kappa2
    Sigma2(i) = rigammaRcpp(Shape, Rate);

  //End loop over locations  
  }
  
  //Update dependent parameters 
  arma::mat Sigma = arma::diagmat(Sigma2);
  arma::mat SigmaInv = arma::diagmat(1 / Sigma2);
  arma::mat BigPsi = Lambda * arma::trans(Lambda) + Sigma;
    
  //Update parameters object
  Para.Sigma2 = Sigma2;
  Para.Sigma = Sigma;
  Para.SigmaInv = SigmaInv;
  Para.BigPsi = BigPsi;
  return Para;
  
}



//Function to sample new value of psi using a Metropolis sampler step-----------------------------------------------
std::pair<para, metrobj> SamplePsi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj) {
  
  //Set data objects
  arma::mat TimeDist = DatObj.TimeDist;
  int Nu = DatObj.Nu;
  int TempCorInd = DatObj.TempCorInd;
  arma::mat EyeKbyNu = DatObj.EyeKbyNu;
  arma::colvec ZeroKbyNu = DatObj.ZeroKbyNu;
  arma::mat EyeNu = DatObj.EyeNu;
  
  //Set parameter objects
  double Psi = Para.Psi;
  arma::mat HPsi = Para.HPsi;
  arma::mat CholHPsi = Para.CholHPsi;
  arma::mat Upsilon = Para.Upsilon;
  arma::colvec Eta = Para.Eta;

  //Set hyperparameter objects
  double APsi = HyPara.APsi;
  double BPsi = HyPara.BPsi;
  double Gamma = HyPara.Gamma;
  double Beta = HyPara.Beta;
  
  //Set metropolis objects
  double MetropPsi = sqrt(MetrObj.MetropPsi);
  double AcceptancePsi = MetrObj.AcceptancePsi;
  
  //Transform current state to real line
  double BigDelta = log((Psi - APsi) / (BPsi - Psi));
  
  //Numerical fix for when the propopsal cholesky doesn't exist
  double PsiProposal, BigDeltaProposal;
  arma::mat CholHPsiProposal(Nu, Nu), HPsiProposal(Nu, Nu);
  bool Cholesky = false;
  while (!Cholesky) {
    
    //Sample a new Proposal
    BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropPsi));
    
    //Compute Phi Proposal
    PsiProposal = (BPsi * exp(BigDeltaProposal) + APsi) / (1 + exp(BigDeltaProposal));
    
    //Fix numerical issue where PsiProposal can equal APsi or BPsi
    // arma::vec PsiProposalVec(1), APsiVec(1), BPsiVec(1);
    // PsiProposalVec(0) = PsiProposal;
    // APsiVec(0) = APsi;
    // BPsiVec(0) = BPsi;
    // double TOL = 0.000001;
    // if ((rows_equal(PsiProposalVec, APsiVec, TOL)) || (rows_equal(PsiProposalVec, BPsiVec, TOL))) {
    //   if (rows_equal(PsiProposalVec, APsiVec, TOL)) PsiProposal *= 1.1; //doesn't work when APsi is negative
    //   if (rows_equal(PsiProposalVec, BPsiVec, TOL)) PsiProposal *= 0.99;
    //   BigDeltaProposal = log((PsiProposal - APsi) / (BPsi - PsiProposal));
    // }
    
    //Proposal temporal correlation
    HPsiProposal = H(PsiProposal, TempCorInd, TimeDist, Nu);
    Cholesky = arma::chol(CholHPsiProposal, HPsiProposal);
    
  }
  
  //Eta structure components
  arma::mat CholUpsilon = arma::chol(Upsilon);
  arma::mat RootiEta = arma::solve(arma::trimatu(arma::kron(CholHPsi, CholUpsilon)), EyeKbyNu);
  arma::mat RootiEtaProposal = arma::solve(arma::trimatu(arma::kron(CholHPsiProposal, CholUpsilon)), EyeKbyNu);
  double Component1A = lndMvn(Eta, ZeroKbyNu, RootiEtaProposal);
  double Component1B = lndMvn(Eta, ZeroKbyNu, RootiEta);
  double Component1 = Component1A - Component1B;
  
  //Prior component
  double Component2 = 0; //exponential
  if (TempCorInd == 1) { //ar1
    double Component2A = (Gamma - 1) * log(1 + PsiProposal) + (Beta - 1) * log(1 - PsiProposal);
    double Component2B = (Gamma - 1) * log(1 + Psi) + (Beta - 1) * log(1 - Psi);
    Component2 = Component2A - Component2B;
  }
  
  //Jacobian component 1
  double Component3A = BigDeltaProposal;
  double Component3B = BigDelta;
  double Component3 = Component3A - Component3B;
  
  //Jacobian component 2
  double Component4 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
  
  //Compute log acceptance ratio
  double LogR = Component1 + Component2 + Component3 + Component4;
  
  //Metropolis update
  double RandU = randuRcpp();
  if (log(RandU) < LogR) {
    
    //Keep Count of Acceptances
    AcceptancePsi++;
    MetrObj.AcceptancePsi = AcceptancePsi;
    
    //Update dependent parameters
    arma::mat P = arma::solve(arma::trimatu(CholHPsiProposal), EyeNu);
    
    //Update parameters object
    Para.Psi = PsiProposal;
    Para.HPsi = HPsiProposal;
    Para.CholHPsi = CholHPsiProposal;
    Para.HPsiInv = P * arma::trans(P);

  }
  
  //Return output object
  return std::pair<para, metrobj>(Para, MetrObj);
  
}



//Function to sample Upsilon using a Gibbs sampler step-------------------------------------------------------------------
para SampleUpsilon(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int Nu = DatObj.Nu;

  //Set parameters
  arma::mat BigPhi = Para.BigPhi;
  arma::mat HPsiInv = Para.HPsiInv;
  
  //Set hyperparameter objects
  double Zeta = HyPara.Zeta;
  arma::mat Omega = HyPara.Omega;
  
  //Compute SThetaKappa
  arma::mat SPhiPsi = BigPhi * HPsiInv * arma::trans(BigPhi);
  
  //Sample T
  double n = Zeta + Nu;
  arma::mat V = SPhiPsi + Omega;
  arma::mat UpsilonInv = rwishRcpp(n, CholInv(V));
  arma::mat Upsilon = CholInv(UpsilonInv);
  
  //Update parameters object
  Para.Upsilon = Upsilon;
  Para.UpsilonInv = UpsilonInv;
  return Para;
}



//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
para SampleEta(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::colvec YStar = DatObj.YStar;
  int K = DatObj.K;
  int Nu = DatObj.Nu;
  
  //Set parameters
  arma::mat Lambda = Para.Lambda;
  arma::mat SigmaInv = Para.SigmaInv;
  arma::mat HPsiInv = Para.HPsiInv;
  arma::mat UpsilonInv = Para.UpsilonInv;
  
  //Sample eta
  arma::mat tLambdaSigmaInv = arma::trans(Lambda) * SigmaInv;
  arma::mat CovEta = CholInv(arma::kron(EyeNu, tLambdaSigmaInv * Lambda) + arma::kron(HPsiInv, UpsilonInv));
  arma::colvec MeanEta = CovEta * arma::kron(EyeNu, tLambdaSigmaInv) * YStar;
  arma::colvec Eta = rmvnormRcpp(1, MeanEta, CovEta);
  
  //Update parameters dependent on delta
  arma::mat BigPhi = arma::reshape(Eta, K, Nu);
  arma::mat Mean = arma::kron(EyeNu, Lambda) * Eta;
  
  //Update parameters object
  Para.Eta = Eta;
  Para.BigPhi = BigPhi;
  Para.Mean = Mean;
  return Para;
  
}



//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
para SampleDelta(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int K = DatObj.K;
  int L = DatObj.L;

  //Set parameter objects
  arma::mat Theta = Para.Theta;
  arma::colvec Delta = Para.Delta;
  
  //Set hyperparameter objects
  double A1 = HyPara.A1;
  double A2 = HyPara.A2;
  
  //Loop over K delta precision parameters
  double AH;
  for (arma::uword h = 0; h < K; h++) {
    
    //Remove hth delta
    arma::colvec DeltaMinusH = Delta;
    DeltaMinusH(h) = 1;
    
    //Asign hyperparameter
    if (h == 0) AH = A1;
    if (h > 0) AH = A2;
    
    //Obtain moments
    double Resids = 0;
    for (arma::uword j = h; j < K; j++) {
      arma::colvec ThetaJ = Theta.col(j);
      Resids += arma::as_scalar(arma::trans(ThetaJ) * ThetaJ * arma::prod(DeltaMinusH(arma::span(0, j))));
    }
    double Shape = AH + 0.5 * (K - h) * L;
    double Rate = 1 + 0.5 * Resids;
    
    //Sample Deltah
    Delta(h) = rgammaRcpp(Shape, Rate);

  //End loop over deltas  
  }
  
  //Update parameters object
  Para.Delta = Delta;
  Para.Tau = arma::cumprod(Delta);
  return Para;
  
}



//Function to sample kappa2 using a Gibbs sampler step---------------------------------------------------------------
para SampleKappa2(datobj DatObj, para Para, hypara HyPara) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int L = DatObj.L;
  arma::mat ICAR = DatObj.ICAR;

  //Set parameter objects
  arma::cube Alpha = Para.Alpha;
  
  //Set hyperparameter objects
  double C = HyPara.C;
  double D = HyPara.D;
  
  //Calculate moments
  double Resids = 0;
  for (arma::uword j = 0; j < K; j++) {
    for (arma::uword l = 0; l < L; l++) {
      arma::rowvec AlphaJL = Alpha.slice(j).row(l);
      Resids += arma::as_scalar(AlphaJL * ICAR * arma::trans(AlphaJL));
    }
  }
  double Shape = C + 0.5 * K * L * M;
  double Rate = D + 0.5 * Resids;
  
  //Sample kappa2
  double Kappa2 = rigammaRcpp(Shape, Rate);

  //Update parameters object
  Para.Kappa2 = Kappa2;
  return Para;
  
}



//Function to sample alphajl(s_i)'s using a Gibbs sampler step---------------------------------------------------------------
para SampleAlpha(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int L = DatObj.L;
  arma::mat ICAR = DatObj.ICAR;
  arma::mat EyeM = DatObj.EyeM;
  
  //Set parameter objects
  double Kappa2 = Para.Kappa2;
  arma::cube Z = Para.Z;
  arma::cube Alpha = Para.Alpha;
  
  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Get jth process objects
    arma::mat AlphaJ = Alpha.slice(j); 
    arma::mat ZJ = Z.slice(j);
    
    //Loop over clusters L
    for (arma::uword l = 0; l < L; l++) {
      
      //Sample AlphaJL
      arma::mat CovAlpha = CholInv(EyeM + ICAR / Kappa2);
      arma::colvec MeanAlpha = CovAlpha * arma::trans(ZJ.row(l));
      arma::colvec AlphaJL = rmvnormRcpp(1, MeanAlpha, CovAlpha);
      // AlphaJL = (AlphaJL - arma::mean(AlphaJL)); // center
      AlphaJ.row(l) = arma::trans(AlphaJL);  
        
    //End loop over clusters  
    }
      
    //Update Z
    Alpha.slice(j) = AlphaJ;
    
  //End loop over columns 
  }
  
  //Update Weights
  arma::cube Weights = GetWeights(Alpha, K, M, L);
  arma::cube logWeights = GetlogWeights(Alpha, K, M, L);
  
  //Update parameters object
  Para.Alpha = Alpha;
  Para.logWeights = logWeights;
  Para.Weights = Weights;
  return Para;
  
}



//Function to sample zjl(s_i)'s using a Gibbs sampler step---------------------------------------------------------------
para SampleZ(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int L = DatObj.L;
  
  //Set parameter objects
  arma::cube Alpha = Para.Alpha;
  arma::umat Xi = Para.Xi;
  arma::cube Z = Para.Z;
  
  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Get jth process objects
    arma::mat AlphaJ = Alpha.slice(j);

    //Loop over locations M
    arma::mat ZJ(L, M);
    for (arma::uword i = 0; i < M; i++) {
      
      //Loop over clusters L
      for (arma::uword l = 0; l < L; l++) {
        
        //Sample zjl(s_i)
        arma::uword XiIJ = Xi(i, j);
        double ZIJL;
        if (XiIJ > l) ZIJL = rtnormRcppMSM(AlphaJ(l, i), 1, -arma::datum::inf, 0);
        if (XiIJ == l) ZIJL = rtnormRcppMSM(AlphaJ(l, i), 1, 0, arma::datum::inf);
        ZJ(l, i) = ZIJL;
      
      //End loop over clusters  
      }
    
    //End loop over locations
    }
    
    //Update Z
    Z.slice(j) = ZJ;
  
  //End loop over columns 
  }
  
  //Update parameters object
  Para.Z = Z;
  return Para;
  
}



//Function to sample xi's using a Gibbs sampler step---------------------------------------------------------------
para SampleXi(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int L = DatObj.L;
  arma::Col<int> SeqL = DatObj.SeqL;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeNu = DatObj.EyeNu;
  
  //Set parameter objects
  arma::umat Xi = Para.Xi;
  arma::cube logWeights = Para.logWeights;
  arma::mat Lambda = Para.Lambda;
  arma::mat BigPhi = Para.BigPhi;
  arma::colvec Sigma2 = Para.Sigma2;
  arma::mat Theta = Para.Theta;
  arma::mat Sigma = Para.Sigma;
  arma::colvec Eta = Para.Eta;
  
  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Get jth process objects
    arma::mat logWeightsJ = logWeights.slice(j);
    arma::colvec ThetaJ = Theta.col(j);
    
    //Loop over locations M
    for (arma::uword i = 0; i < M; i++) {
      
      //Get location specific objects
      arma::colvec logWeightsIJ = logWeightsJ.col(i);
      double Sigma2I = arma::as_scalar(Sigma2.row(i));
      arma::rowvec LambdaI = Lambda.row(i);
      
      //Loop over clusters L
      arma::vec logProbsRaw(L);
      for (arma::uword l = 0; l < L; l++) {
        
        //Obtain the un-normalized (raw) probabilities on the log scale
        LambdaI(j) = ThetaJ(l);
        arma::rowvec Resid = YStarWide.row(i) - LambdaI * BigPhi;
        double ResidQ = arma::as_scalar(Resid * arma::trans(Resid));
        double Likelihood = -0.5 * (ResidQ / Sigma2I);
        double logWeightsIJL = logWeightsIJ(l);
        logProbsRaw(l) = logWeightsIJL + Likelihood;
        
      }
      
      //Use log sum exponential trick to get normalized probabilities
      double Max = arma::max(logProbsRaw);
      double Delta = Max + log(arma::sum(arma::exp(logProbsRaw - Max)));
      arma::vec ProbsIJ = arma::exp(logProbsRaw - Delta);
      
      //Sample a new label
      arma::uword XiIJ = arma::as_scalar(sampleRcpp(SeqL, 1, true, ProbsIJ));
      Xi(i, j) = XiIJ;
      
      //Update Lambda
      Lambda(i, j) = Theta(XiIJ, j);
      
    }
  }
  
  //Update parameters object
  Para.Xi = Xi;
  Para.Lambda = Lambda;
  Para.BigPsi = Lambda * arma::trans(Lambda) + Sigma;
  Para.Mean = arma::kron(EyeNu, Lambda) * Eta;
  return Para;
  
}


//Function to sample theta_jl's using a Gibbs sampler step---------------------------------------------------------------
para SampleTheta(datobj DatObj, para Para) {
  
  //Set data objects
  int K = DatObj.K;
  int L = DatObj.L;
  arma::mat YStarWide = DatObj.YStarWide;
  arma::mat EyeNu = DatObj.EyeNu;
  
  //Set parameter objects
  arma::mat BigPhi = Para.BigPhi;
  arma::umat Xi = Para.Xi;
  arma::mat Lambda = Para.Lambda;
  arma::vec Sigma2 = Para.Sigma2;
  arma::vec Tau = Para.Tau;
  arma::mat Theta = Para.Theta;
  arma::colvec Eta = Para.Eta;
  arma::mat Sigma = Para.Sigma;
  
  //Loop over columns K
  for (arma::uword j = 0; j < K; j++) {
    
    //Column specific objects
    arma::colvec EtaJ = arma::trans(BigPhi.row(j));
    arma::uvec XiJ = Xi.col(j);
    double TauJ = arma::as_scalar(Tau.row(j));
    arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
    LambdaMinusJ.shed_col(j);
    BigPhiMinusJ.shed_row(j);
    
    //Loop over clusters L
    for (arma::uword l = 0; l < L; l++) {
      
      //Number of objects in cluster l
      arma::uvec WhichJL = find(XiJ == l);
      int NJL = WhichJL.size();
      
      //Declarations
      double ThetaJL;
      
      //Sample from prior if none belong to cluster l
      if (NJL == 0) ThetaJL = arma::as_scalar(rnormRcpp(1, 0, sqrt(1 / TauJ)));
      
      //Sample from posterior if at least 1 belongs to cluster l
      if (NJL > 0) {
        
        //Cluster l specific objects
        arma::colvec Sigma2InvJL = 1 / Sigma2(WhichJL);
        arma::mat YJL = YStarWide.rows(WhichJL);
        arma::mat LambdaJL = LambdaMinusJ.rows(WhichJL);
        
        //Sample theta
        double ResidsJL = arma::as_scalar(arma::trans((YJL - LambdaJL * BigPhiMinusJ) * EtaJ) * Sigma2InvJL);
        double VarThetaJL = arma::as_scalar(1 / ((arma::trans(EtaJ) * EtaJ) * arma::sum(Sigma2InvJL) + TauJ));
        double MeanThetaJL = VarThetaJL * ResidsJL;
        ThetaJL = arma::as_scalar(rnormRcpp(1, MeanThetaJL, sqrt(VarThetaJL)));
        
      }

      //Update parameters
      Theta(l, j) = ThetaJL;
      arma::colvec LambdaJ = Lambda.col(j);
      arma::vec ThetaJLVec(1);
      ThetaJLVec(0) = ThetaJL;
      LambdaJ(WhichJL) = arma::repmat(ThetaJLVec, NJL, 1);
      Lambda.col(j) = LambdaJ;

    }
  }
  
  //Final updates
  arma::colvec Mean = arma::kron(EyeNu, Lambda) * Eta;
  arma::mat BigPsi = Lambda * arma::trans(Lambda) + Sigma;
  
  //Update parameters object
  Para.Theta = Theta;
  Para.Lambda = Lambda;
  Para.Mean = Mean;
  Para.BigPsi = BigPsi;
  return Para;
  
}