###Function for reading in sampler inputs and creating a list object that contains all relavent data objects--------------------
CreateDatObj <- function(Y, W, Time, K, L, Rho, ScaleY, Family, TemporalStructure) {

  ###Data objects
  YObserved <- Y / ScaleY #scale observed data
  N <- length(YObserved)  #total observations
  M <- dim(W)[1] #number of spatial locations
  Nu <- N / M #number of visits

  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(YObserved, ncol = 1)
  YStarWide <- matrix(YStar, nrow = M, ncol = Nu)

  ###CAR covariance objects
  Dw <- diag(apply(W, 1, sum))
  ICAR <- Dw - Rho * W
  ICARInv <- CholInv(ICAR)
  
  ###Temporal distance matrix
  TimeDist <- abs(outer(Time, Time, "-"))
  
  ###Matrix Objects
  EyeNu <- diag(Nu)
  EyeM <- diag(M)
  SeqL <- matrix(0:(L - 1), ncol = 1)
  EyeKbyNu <- diag(K * Nu)
  ZeroKbyNu <- matrix(0, nrow = K * Nu)
  OneNu <- matrix(1, nrow = Nu)
  
  ###Assign temporal correlation structure
  if (TemporalStructure == "exponential") TempCorInd <- 0
  if (TemporalStructure == "ar1") TempCorInd <- 1
  
  ###Family indicator
  if (Family == "normal") FamilyInd <- 0
  if (Family == "probit") FamilyInd <- 1
  if (Family == "tobit") FamilyInd <- 2

  ###Make parameters global
  DatObj <- list()
  DatObj$YObserved <- YObserved
  DatObj$ScaleY <- ScaleY
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$W <- W
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$K <- K
  DatObj$L <- L
  DatObj$FamilyInd <- FamilyInd
  DatObj$Time <- Time
  DatObj$Rho <- Rho
  DatObj$ICAR <- ICAR
  DatObj$ICARInv <- ICARInv
  DatObj$TempCorInd <- TempCorInd
  DatObj$TimeDist <- TimeDist
  DatObj$EyeNu <- EyeNu
  DatObj$SeqL <- SeqL
  DatObj$EyeM <- EyeM
  DatObj$EyeKbyNu <- EyeKbyNu
  DatObj$ZeroKbyNu <- ZeroKbyNu
  DatObj$OneNu <- OneNu
  return(DatObj)

}




###Function to create Hyperparameter Object------------------------------------------------------------------------------------
CreateHyPara <- function(Hypers, DatObj) {

  ###Set data objects
  K <- DatObj$K
  TempCorInd <- DatObj$TempCorInd
  TimeDist <- DatObj$TimeDist 
  
  ###Which parameters are user defined?
  UserHypers <- names(Hypers)

  ###Set hyperparameters for Sigma2
  if ("Sigma2" %in% UserHypers) {
    A <- Hypers$Sigma2$A
    B <- Hypers$Sigma2$B
  }
  if (!("Sigma2" %in% UserHypers)) {
    A <- 1
    B <- 1
  }
  
  ###Set hyperparameters for Kappa2
  if ("Kappa2" %in% UserHypers) {
    C <- Hypers$Kappa2$C
    D <- Hypers$Kappa2$D
  }
  if (!("Kappa2" %in% UserHypers)) {
    C <- 1
    D <- 1
  }

  ###Set hyperparameters for Delta1
  if ("Delta1" %in% UserHypers) {
    A1 <- Hypers$Delta1$A1
  }
  if (!("Delta1" %in% UserHypers)) {
    A1 <- 1
  }
  
  ###Set hyperparameters for Deltah
  if ("Deltah" %in% UserHypers) {
    A2 <- Hypers$Deltah$A2
  }
  if (!("Deltah" %in% UserHypers)) {
    A2 <- 1
  }

  ###Set hyperparameters for Psi
  if ("Psi" %in% UserHypers) {
    if (TempCorInd == 0) { # exponential
      APsi <- Hypers$Psi$APsi
      BPsi <- Hypers$Psi$BPsi
      Gamma <- 1 #Null values
      Beta <- 1 #Null values
    }
    if (TempCorInd == 1) { # ar1
      APsi <- -1
      BPsi <- 1
      Gamma <- Hypers$Psi$Gamma
      Beta <- Hypers$Psi$Beta
    }
  }
  if (!"Psi" %in% UserHypers) {
    if (TempCorInd == 0) { # exponential
      minDiff <- min(TimeDist[TimeDist > 0])
      maxDiff <- max(TimeDist[TimeDist > 0])
      APsi <- -log(0.95) / maxDiff #longest diff goes up to 95%
      BPsi <- -log(0.01) / minDiff #shortest diff goes down to 1%
      Gamma <- 1 #Null values
      Beta <- 1 #Null values
    }
    if (TempCorInd == 1) { # ar1
      APsi <- -1
      BPsi <- 1
      Gamma <- 1
      Beta <- 1
    }
  }
  
  ###Set hyperparameters for Upsilon
  if ("Upsilon" %in% UserHypers) {
    Zeta <- Hypers$Upsilon$Zeta
    Omega <- Hypers$Upsilon$Omega
  }
  if (!("Upsilon" %in% UserHypers)) {
    Zeta <- K + 1
    Omega <- diag(K)
  }

  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$A <- A
  HyPara$B <- B
  HyPara$C <- C
  HyPara$D <- D
  HyPara$A1 <- A1
  HyPara$A2 <- A2
  HyPara$APsi <- APsi
  HyPara$BPsi <- BPsi
  HyPara$Gamma <- Gamma
  HyPara$Beta <- Beta
  HyPara$Zeta <- Zeta
  HyPara$Omega <- Omega
  return(HyPara)

}



###Function for creating an object containing relevant Metropolis information---------------------------------------------------
CreateMetrObj <- function(Tuning) {

  ###Which parameters are user defined?
  UserTuners <- names(Tuning)

  ###Set tuning parameters for Psi
  if ("Psi" %in% UserTuners) MetropPsi <- Tuning$Psi
  if (!("Psi" %in% UserTuners)) MetropPsi <- 1

  ###Set acceptance rate counters
  AcceptancePsi <- 0

  ###Return metropolis object
  MetrObj <- list()
  MetrObj$MetropPsi <- MetropPsi
  MetrObj$AcceptancePsi <- AcceptancePsi
  MetrObj$OriginalTuners <- MetropPsi
  return(MetrObj)

}



###Function for creating inital parameter object-------------------------------------------------------------------------------
CreatePara <- function(Starting, DatObj, HyPara) {

  ###Set data objects
  K <- DatObj$K
  M <- DatObj$M
  L <- DatObj$L
  Nu <- DatObj$Nu
  TempCorInd <- DatObj$TempCorInd
  TimeDist <- DatObj$TimeDist
  EyeNu <- DatObj$EyeNu
  
  ###Set hyperparameter objects
  APsi <- HyPara$APsi
  BPsi <- HyPara$BPsi
  
  ###Which parameters are user defined?
  UserStarters <- names(Starting)

  ###Set initial values of Sigma2
  if ("Sigma2" %in% UserStarters) Sigma2 <- matrix(Starting$Sigma2, ncol = 1, nrow = M)
  if ((!"Sigma2" %in% UserStarters)) Sigma2 <- matrix(1, ncol = 1, nrow = M)

  ###Set initial values of Kappa2
  if ("Kappa2" %in% UserStarters) Kappa2 <- Starting$Kappa2
  if ((!"Kappa2" %in% UserStarters)) Kappa2 <- 1
  
  ###Set initial values of Delta
  if ("Delta" %in% UserStarters) Delta <- matrix(Starting$Delta, nrow = K, ncol = 1)
  if ((!"Delta" %in% UserStarters)) Delta <- matrix(1, nrow = K, ncol = 1)
  
  ###Set initial values of Psi
  if ("Psi" %in% UserStarters) {
    if (TempCorInd == 0) { #exponential
      Psi <- Starting$Psi
      if ((Psi <= APsi) | (Psi >= BPsi)) stop('Starting: "Psi" must be in (APsi, BPsi)')
    }
    if (TempCorInd == 1) { #ar1
      Psi <- Starting$Psi
      if ((Psi <= -1) | (Psi >= 1)) stop('Starting: "Psi" must be in (-1, 1)')
    }
  }
  if ((!"Psi" %in% UserStarters)) {
    if (TempCorInd == 0) { #exponential
      Psi <- mean(c(APsi, BPsi))
    }
    if (TempCorInd == 1) { #ar1
      Psi <- 0
    }
  }
    
  ###Set initial value of Upsilon
  if ("Upsilon" %in% UserStarters) Upsilon <- matrix(Starting$Upsilon, nrow = K, ncol = K)
  if (!("Upsilon" %in% UserStarters)) Upsilon <- diag(K)

  ###Create atom variances (Tau2) and atoms themselves (Theta)
  Tau <- matrix(cumprod(Delta), nrow = K, ncol = 1)
  Theta <- apply(sqrt(1 / Tau), 1, function(x) rnorm(L, 0, x))
  
  ###Create label parameters
  Xi <- matrix(0, nrow = M, ncol = K)
  Lambda <- matrix(nrow = M, ncol = K)
  for (i in 1:M) {
    for (j in 1:K) {
      Lambda[i, j] <- Theta[Xi[i, j] + 1, j]
    }
  }
  
  ###Factors
  BigPhi <- matrix(0, nrow = K, ncol = Nu)
  Eta <- matrix(as.numeric(BigPhi), ncol = 1)
  
  ###Probit parameters
  Alpha <- Weights <- logWeights <- array(0, dim = c(L, M, K))
  Z <- array(-1, dim = c(L, M, K))
  Z[1, , ] <- 1
  logUpperPhiAlpha <- pnorm(Alpha, lower.tail = FALSE, log.p = TRUE)
  UpperPhiAlpha <- pnorm(Alpha, lower.tail = FALSE)
  for (j in 1:K) {
    for (i in 1:M) {
      for (l in 1:L) {
        if (l == 1) {
          Weights[l, i, j] <- pnorm(Alpha[l, i, j])
          logWeights[l, i, j] <- pnorm(Alpha[l, i, j], log.p = TRUE)
        }
        if (l > 1) {
          Weights[l, i, j] <- pnorm(Alpha[l, i, j]) * prod(UpperPhiAlpha[1:(l - 1) , i, j])
          logWeights[l, i, j] <- pnorm(Alpha[l, i, j], log.p = TRUE) + sum(logUpperPhiAlpha[1:(l - 1), i, j])
        }
      }
    }
  }

  ###Marginal covariance
  Sigma <- diag(as.numeric(Sigma2))
  SigmaInv <- diag(as.numeric(1 / Sigma2))
  BigPsi <- Lambda %*% t(Lambda) + Sigma
  
  ###Temporal parameters
  HPsi <- H(Psi, TempCorInd, TimeDist, Nu)
  CholHPsi <- chol(HPsi)
  HPsiInv <- chol2inv(CholHPsi)
  
  ###Save parameter objects
  Para <- list()
  Para$Sigma2 <- Sigma2
  Para$Kappa2 <- Kappa2
  Para$Delta <- Delta
  Para$Psi <- Psi
  Para$Upsilon <- Upsilon
  Para$UpsilonInv <- CholInv(Upsilon)
  Para$Xi <- Xi
  Para$Theta <- Theta
  Para$Lambda <- Lambda
  Para$Tau <- Tau
  Para$BigPhi <- BigPhi
  Para$Eta <- Eta
  Para$Alpha <- Alpha
  Para$Z <- Z
  Para$BigPsi <- BigPsi
  Para$Sigma <- Sigma
  Para$SigmaInv <- SigmaInv
  Para$HPsi <- HPsi
  Para$CholHPsi <- CholHPsi
  Para$HPsiInv <- HPsiInv
  Para$Mean <- kronecker(EyeNu, Lambda) %*% Eta
  Para$Weights <- Weights
  Para$logWeights <- logWeights
  return(Para)

}



###Function that creates the data augmentation (i.e. Tobit) booleans------------------------------------------------------------
CreateDatAug <- function(DatObj) {

  ###Set data object
  YObserved <- DatObj$YObserved
  FamilyInd <- DatObj$FamilyInd
  YStarWide <- DatObj$YStarWide
  M <- DatObj$M
  Nu <- DatObj$Nu

  ###Initialize Data Augmentation Object
  DatAug <- NULL

  ###Normal objects
  if (FamilyInd == 0) {
    DatAug$NBelow <- 0
    DatAug$NAbove <- 0
    DatAug$WhichBelow <- 0
    DatAug$WhichAbove <- 0
    # DatAug$TobitIndeces <- matrix(c(0,0,0,0),2,2)
    # DatAug$ProbitIndeces <- matrix(c(0,0,0,0),2,2)
  }

  ###Probit objects
  if (FamilyInd == 1) {
    TobitBoolean <- YObserved <= 0
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    TobitBooleanMat <- matrix(TobitBoolean, nrow = M, ncol = Nu)
    YStarBelow <- list()
    for (i in 1 : Nu) YStarBelow[[i]] <- YStarWide[!TobitBooleanMat[,i], i]
    NBelowList <- unlist(lapply(YStarBelow, f<-function(x) M - length(x)))
    TobitIndeces <- which(TobitBooleanMat, arr.ind = TRUE)
    TobitIndeces <- TobitIndeces - 1
    ProbitBoolean <- YObserved > 0
    WhichAbove <- which(ProbitBoolean)
    NAbove <- length(WhichAbove)
    ProbitBooleanMat <- matrix(ProbitBoolean, nrow = M, ncol = Nu)
    YStarAbove <- list()
    for (i in 1 : Nu) YStarAbove[[i]] <- YStarWide[!ProbitBooleanMat[,i], i]
    NAboveList <- unlist(lapply(YStarAbove, f<-function(x) M - length(x)))
    ProbitIndeces <- which(ProbitBooleanMat, arr.ind = TRUE)
    ProbitIndeces <- ProbitIndeces - 1

    ###Save objects
    DatAug <- list()
    DatAug$WhichBelow <- WhichBelow - 1
    DatAug$WhichAbove <- WhichAbove - 1
    DatAug$NBelow <- NBelow
    DatAug$NAbove <- NAbove
    DatAug$TobitIndeces <- TobitIndeces
    DatAug$ProbitIndeces <- ProbitIndeces
  }

  ###Tobit objects
  if (FamilyInd == 2) {
    TobitBoolean <- YObserved <= 0
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    TobitBooleanMat <- matrix(TobitBoolean, nrow = M, ncol = Nu)
    YStarNonZero <- list()
    for (i in 1 : Nu) YStarNonZero[[i]] <- YStarWide[!TobitBooleanMat[,i], i]
    NBelowCount <- unlist(lapply(YStarNonZero, f<-function(x) M - length(x)))
    TobitIndeces <- which(TobitBooleanMat, arr.ind = TRUE) - 1
    # TobitIndeces <- TobitIndeces - 1
    # ZDatAug <- model.matrix(~-1 + as.factor(TobitIndeces[,2]))
    # attributes(ZDatAug) <- NULL
    # ZDatAug <- structure(ZDatAug, class = "matrix", dim = c(NBelow, Nu))
    # WDatAug <- array(FALSE, dim = c(M, M, Nu))
    # for (i in 1:NBelow) {
    #   Visit <- TobitIndeces[i, 2] + 1
    #   Location <- TobitIndeces[i, 1] + 1
    #   WDatAug[ , Location, Visit] <- rep(TRUE, M)
    # }
    # WDatAug <- matrix(which(WDatAug) - 1, ncol = 1)

    ###Save objects
    DatAug <- list()
    DatAug$WhichBelow <- WhichBelow - 1
    DatAug$NBelow <- NBelow
    DatAug$TobitBooleanMat <- TobitBooleanMat
    DatAug$YStarNonZero <- YStarNonZero
    DatAug$NBelowCount <- NBelowCount
    DatAug$TobitIndeces <- TobitIndeces
    DatAug$ProbitIndeces <- matrix(c(0,0,0,0),2,2)
    DatAug$NAbove <- 0
    DatAug$WhichAbove <- 0
    # DatAug$ZDatAug <- ZDatAug
    # DatAug$WDatAug <- WDatAug
  }
  return(DatAug)

}



###Function that creates inputs for MCMC sampler--------------------------------------------------------------------------------
CreateMcmc <- function(MCMC, DatObj) {

  ###Which parameters are user defined?
  UserMCMC <- names(MCMC)

  ###Set MCMC objects
  if ("NBurn" %in% UserMCMC) NBurn <- MCMC$NBurn
  if (!("NBurn" %in% UserMCMC)) NBurn <- 10000
  if ("NSims" %in% UserMCMC) NSims <- MCMC$NSims
  if (!("NSims" %in% UserMCMC)) NSims <- 100000
  if ("NThin" %in% UserMCMC) NThin <- MCMC$NThin
  if (!("NThin" %in% UserMCMC)) NThin <- 10
  if ("NPilot" %in% UserMCMC) NPilot <- MCMC$NPilot
  if (!("NPilot" %in% UserMCMC)) NPilot <- 20

  ###One last check of MCMC user inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!(is.wholenumber(NSims / NThin))) stop('MCMC: "NThin" must be a factor of "NSims"')
  if (!(is.wholenumber(NBurn / NPilot))) stop('MCMC: "NPilot" must be a factor of "NBurn"')

  ###Create MCMC objects
  NTotal <- NBurn + NSims
  WhichKeep <- NBurn + (1:(NSims / NThin)) * NThin
  NKeep <- length(WhichKeep)

  ###Pilot adaptation objects
  WhichPilotAdapt <- (1:NPilot) * NBurn / NPilot
  PilotAdaptDenominator <- WhichPilotAdapt[1]

  ###Burn-in progres bar
  BarLength <- 50 #Burn-in bar length (arbitrary)
  BurnInProgress <- seq(1 / BarLength, 1, 1 / BarLength)
  WhichBurnInProgress <- sapply(BurnInProgress, function(x) tail(which(1 : NBurn <= x * NBurn), 1))

  ###Progress output objects
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichSamplerProgress <- sapply(SamplerProgress, function(x) tail(which(1:NSims <= x * NSims), 1)) + NBurn
  WhichBurnInProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NBurn <= x * NBurn), 1))

  ###Save objects
  MCMC <- list()
  MCMC$NBurn <- NBurn
  MCMC$NSims <- NSims
  MCMC$NThin <- NThin
  MCMC$NPilot <- NPilot
  MCMC$NTotal <- NTotal
  MCMC$WhichKeep <- WhichKeep
  MCMC$NKeep <- NKeep
  MCMC$WhichPilotAdapt <- WhichPilotAdapt
  MCMC$PilotAdaptDenominator <- PilotAdaptDenominator
  MCMC$BurnInProgress <- BurnInProgress
  MCMC$WhichBurnInProgress <- WhichBurnInProgress
  MCMC$WhichBurnInProgressInt <- WhichBurnInProgressInt
  MCMC$BarLength <- BarLength
  MCMC$WhichSamplerProgress <- WhichSamplerProgress
  return(MCMC)

}



###Function that creates a storage object for raw samples-----------------------------------------------------------------------
CreateStorage <- function(DatObj, McmcObj) {

  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu

  ###Set MCMC objects
  NKeep <- McmcObj$NKeep

  ###Create storage object 
  # Lambda: M x K
  # Eta: K x Nu
  # Sigma: M
  # Kappa2: 1
  # Delta: K
  # Upsilon: K x (K + 1) / 2
  # Psi: 1
  Out <- matrix(nrow = (M * K + K * Nu + M + 1 + K + ((K + 1) * K) / 2 + 1), ncol = NKeep)
  return(Out)

}



