###Function for reading in sampler inputs and creating a list object that contains all relavent data objects--------------------
CreateDatObj <- function(Y, Dist, Time, K, L, ScaleY, Family, TemporalStructure, SpatialStructure) {

  ###Data objects
  M <- dim(Y)[1] #number of spatial locations
  O <- dim(Y)[2] #number of spatial observations
  Nu <- dim(Y)[3] #number of visits
  N <- length(Y) #total observations
  YObserved <- Y / ScaleY #scale observed data

  ###Upper bound of L
  if (is.finite(L)) {
    L <- L
    LInf <- 0
  }
  if (is.infinite(L)) {
    L <- M * O
    LInf <- 1
  }
  
  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(as.numeric(YObserved), ncol = 1)
  YStarWide <- matrix(YStar, nrow = M * O, ncol = Nu)
  
  ###Temporal distance matrix
  TimeDist <- abs(outer(Time, Time, "-"))
  
  ###Matrix Objects
  EyeNu <- diag(Nu)
  EyeM <- diag(M)
  EyeO <- diag(O)
  EyeOM <- diag(O * M)
  SeqL <- matrix(0:(L - 1), ncol = 1)
  EyeKbyNu <- diag(K * Nu)
  ZeroKbyNu <- matrix(0, nrow = K * Nu)
  ZeroM <- matrix(0, nrow = M)
  ZeroOM <- matrix(0, nrow = O * M)
  OneNu <- matrix(1, nrow = Nu)
  OneO <- matrix(1, nrow = O)
  
  ###Assign temporal correlation structure
  if (TemporalStructure == "exponential") TempCorInd <- 0
  if (TemporalStructure == "ar1") TempCorInd <- 1

  ###Assign spatial correlation structure
  if (SpatialStructure == "continuous") SpCorInd <- 0
  if (SpatialStructure == "discrete") SpCorInd <- 1
  
  ###Family indicator
  FamilyInd <- numeric(length = O)
  if (length(Family) == O) {
    for (o in 1:O) {
      if (Family[o] == "normal") FamilyInd[o] <- 0
      if (Family[o] == "probit") FamilyInd[o] <- 1
      if (Family[o] == "tobit") FamilyInd[o] <- 2
    }
  }
  if (length(Family) == 1) {
    if (Family == "normal") FamilyInd <- rep(0, O)
    if (Family == "probit") FamilyInd <- rep(1, O)
    if (Family == "tobit") FamilyInd <- rep(2, O)
  }
    
  ###Make parameters global
  DatObj <- list()
  DatObj$YObserved <- YObserved
  DatObj$ScaleY <- ScaleY
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$SpDist <- Dist
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$K <- K
  DatObj$L <- L
  DatObj$O <- O
  DatObj$FamilyInd <- FamilyInd
  DatObj$Time <- Time
  DatObj$TempCorInd <- TempCorInd
  DatObj$SpCorInd <- SpCorInd
  DatObj$TimeDist <- TimeDist
  DatObj$EyeNu <- EyeNu
  DatObj$EyeO <- EyeO
  DatObj$SeqL <- SeqL
  DatObj$EyeM <- EyeM
  DatObj$EyeOM <- EyeOM
  DatObj$EyeKbyNu <- EyeKbyNu
  DatObj$ZeroKbyNu <- ZeroKbyNu
  DatObj$ZeroM <- ZeroM
  DatObj$ZeroOM <- ZeroOM
  DatObj$OneNu <- OneNu
  DatObj$OneO <- OneO
  DatObj$LInf <- LInf
  return(DatObj)

}



###Function to create Hyperparameter Object------------------------------------------------------------------------------------
CreateHyPara <- function(Hypers, DatObj) {

  ###Set data objects
  K <- DatObj$K
  O <- DatObj$O
  TempCorInd <- DatObj$TempCorInd
  SpCorInd <- DatObj$SpCorInd
  TimeDist <- DatObj$TimeDist 
  SpDist <- DatObj$SpDist
  
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
  
  ###Set hyperparameters for Kappa
  if ("Kappa" %in% UserHypers) {
    SmallUpsilon <- Hypers$Kappa$SmallUpsilon
    BigTheta <- Hypers$Kappa$BigTheta
  }
  if (!("Kappa" %in% UserHypers)) {
    SmallUpsilon <- O + 1
    BigTheta <- diag(O)
  }

  ###Set hyperparameters for Rho
  if ("Rho" %in% UserHypers) {
    if (SpCorInd == 0) { # continuous
      ARho <- Hypers$Rho$ARho
      BRho <- Hypers$Rho$BRho
    }
    if (SpCorInd == 1) { # discrete
      ARho <- 0 #null values, because Rho is fixed
      BRho <- 1
    }
  }
  if (!"Rho" %in% UserHypers) {
    if (SpCorInd == 0) { # continuous
      minDiff <- min(SpDist[SpDist > 0])
      maxDiff <- max(SpDist[SpDist > 0])
      ARho <- -log(0.95) / maxDiff #longest diff goes up to 95%
      BRho <- -log(0.01) / minDiff #shortest diff goes down to 1%
    }
    if (SpCorInd == 1) { # discrete
      ARho <- 0 #null values, because Rho is fixed
      BRho <- 1
    }
  }
  
  ###Set hyperparameters for Delta
  if ("Delta" %in% UserHypers) {
    A1 <- Hypers$Delta$A1
    A2 <- Hypers$Delta$A2
  }
  if (!("Delta" %in% UserHypers)) {
    A1 <- 1
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
  HyPara$SmallUpsilon <- SmallUpsilon
  HyPara$BigTheta <- BigTheta
  HyPara$A1 <- A1
  HyPara$A2 <- A2
  HyPara$APsi <- APsi
  HyPara$BPsi <- BPsi
  HyPara$Gamma <- Gamma
  HyPara$Beta <- Beta
  HyPara$Zeta <- Zeta
  HyPara$Omega <- Omega
  HyPara$ARho <- ARho
  HyPara$BRho <- BRho
  return(HyPara)

}



###Function for creating an object containing relevant Metropolis information---------------------------------------------------
CreateMetrObj <- function(Tuning, DatObj) {

  ###Set data objects
  SpCorInd <- DatObj$SpCorInd
  
  ###Which parameters are user defined?
  UserTuners <- names(Tuning)

  ###Set tuning parameters for Psi
  if ("Psi" %in% UserTuners) MetropPsi <- Tuning$Psi
  if (!("Psi" %in% UserTuners)) MetropPsi <- 1

  ###Set tuning parameters for Rho
  if ("Rho" %in% UserTuners) {
    if (SpCorInd == 0) MetropRho <- Tuning$Rho # continuous
    if (SpCorInd == 1) MetropRho <- 1 # null value
  }
  if (!("Rho" %in% UserTuners)) {
    MetropRho <- 1
  }
  
  ###Set acceptance rate counters
  AcceptancePsi <- AcceptanceRho <- 0

  ###Return metropolis object
  MetrObj <- list()
  MetrObj$MetropPsi <- MetropPsi
  MetrObj$AcceptancePsi <- AcceptancePsi
  MetrObj$MetropRho <- MetropRho
  MetrObj$AcceptanceRho <- AcceptanceRho
  MetrObj$OriginalTuners <- c(MetropPsi, MetropRho)
  return(MetrObj)

}



###Function for creating inital parameter object-------------------------------------------------------------------------------
CreatePara <- function(Starting, DatObj, HyPara) {

  ###Set data objects
  K <- DatObj$K
  M <- DatObj$M
  L <- DatObj$L
  Nu <- DatObj$Nu
  O <- DatObj$O
  TempCorInd <- DatObj$TempCorInd
  TimeDist <- DatObj$TimeDist
  EyeNu <- DatObj$EyeNu
  EyeO <- DatObj$EyeO
  SpCorInd <- DatObj$SpCorInd
  SpDist <- DatObj$SpDist
  
  ###Set hyperparameter objects
  APsi <- HyPara$APsi
  BPsi <- HyPara$BPsi
  ARho <- HyPara$ARho
  BRho <- HyPara$BRho
  
  ###Which parameters are user defined?
  UserStarters <- names(Starting)

  ###Set initial values of Sigma2
  if ("Sigma2" %in% UserStarters) Sigma2 <- matrix(Starting$Sigma2, ncol = 1, nrow = M)
  if ((!"Sigma2" %in% UserStarters)) Sigma2 <- matrix(1, ncol = 1, nrow = M)

  ###Set initial values of Kappa
  if ("Kappa" %in% UserStarters) Kappa <- Starting$Kappa
  if ((!"Kappa" %in% UserStarters)) Kappa <- diag(O)
  
  ###Set initial values of Rho
  if ("Rho" %in% UserStarters) {
    if (SpCorInd == 0) { #continuous
      Rho <- Starting$Rho
      if ((Rho <= ARho) | (Rho >= BRho)) stop('Starting: "Rho" must be in (ARho, BRho)')
    }
    if (SpCorInd == 1) { #discrete
      Rho <- Starting$Rho
    }
  }
  if ((!"Rho" %in% UserStarters)) {
    if (SpCorInd == 0) { #continuous
      Rho <- mean(c(ARho, BRho))
    }
    if (SpCorInd == 1) { #discrete
      Rho <- 0.99
    }
  }
  
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
      if ((Psi <= APsi) | (Psi >= BPsi)) stop('Starting: "Psi" must be in (APsi, BPsi)')
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
  Xi <- matrix(0, nrow = M * O, ncol = K)
  Lambda <- matrix(nrow = M * O, ncol = K)
  for (o in 1:O) {
    for (i in 1:M) {
      for (j in 1:K) {
        Lambda[i + M * (o - 1), j] <- Theta[Xi[i + M * (o - 1), j] + 1, j]
      }
    }  
  }

  ###Factors
  BigPhi <- matrix(0, nrow = K, ncol = Nu)
  Eta <- matrix(as.numeric(BigPhi), ncol = 1)
  
  ###Probit parameters
  Alpha <- Weights <- logWeights <- array(0, dim = c(L, M * O, K))
  Z <- array(-1, dim = c(L, M * O, K))
  Z[1, , ] <- 1
  logUpperPhiAlpha <- pnorm(Alpha, lower.tail = FALSE, log.p = TRUE)
  UpperPhiAlpha <- pnorm(Alpha, lower.tail = FALSE)
  for (j in 1:K) {
    for (o in 1:O) {
      for (i in 1:M) {
        Index <- i + M * (o - 1)
        for (l in 1:L) {
          if (l == 1) {
            Weights[l, Index, j] <- pnorm(Alpha[l, Index, j])
            logWeights[l, Index, j] <- pnorm(Alpha[l, Index, j], log.p = TRUE)
          }
          if (l > 1) {
            Weights[l, Index, j] <- pnorm(Alpha[l, Index, j]) * prod(UpperPhiAlpha[1:(l - 1) , Index, j])
            logWeights[l, Index, j] <- pnorm(Alpha[l, Index, j], log.p = TRUE) + sum(logUpperPhiAlpha[1:(l - 1), Index, j])
          }
        }
      } 
    }
  }

  ###Slice sampling latent parameter
  U <- matrix(nrow = M * O, ncol = K)
  for (j in 1:K) {
    for (o in 1:O) {
      for (i in 1:M) {
        Index <- i + M * (o - 1)
        U[Index, j] <- runif(1, 0, Weights[Xi[Index, j] + 1, Index, j])
      }  
    }
  }
  
  ###Upper Bounds for latent sampling
  OneMinusUStar <- 1 - apply(U, 2, min)
  LStarJ <- numeric(K)
  for (j in 1:K) {
    LStarOIJ <- numeric(M * O)
    for (o in 1:O) {
      for (i in 1:M) {
        Index <- i + M * (o - 1)
        LStarOIJ[Index] <- which.max(cumsum(Weights[ , Index, j]) > OneMinusUStar[j]) - 1
      }
    }
    LStarJ[j] <- max(LStarOIJ)
  }
  LStarJ <- matrix(LStarJ, ncol = 1)

  ###Marginal covariance
  Sigma <- diag(as.numeric(Sigma2))
  SigmaInv <- diag(as.numeric(1 / Sigma2))
  BigPsi <- Lambda %*% Upsilon %*% t(Lambda) + kronecker(EyeO, Sigma)
  
  ###Temporal parameters
  HPsi <- H(Psi, TempCorInd, TimeDist, Nu)
  CholHPsi <- chol(HPsi)
  HPsiInv <- chol2inv(CholHPsi)
  
  ###Spatial covariance objects
  CholKappa <- chol(Kappa)
  KappaInv <- chol2inv(CholKappa)
  
  ###Spatial correlation objects
  if (SpCorInd == 1) { #discrete
    Dw <- diag(apply(SpDist, 1, sum))
    SpCovInv <- Dw - Rho * SpDist
    SpCov <- CholInv(SpCovInv)
    CholSpCov <- matrix(0, nrow = M, ncol = M)
  }
  if (SpCorInd == 0) { #discrete
    SpCov <- exp(-Rho * SpDist)
    CholSpCov <- chol(SpCov)
    SpCovInv <- chol2inv(CholSpCov)
  }
  
  ###Save parameter objects
  Para <- list()
  Para$Sigma2 <- Sigma2
  Para$Kappa <- Kappa
  Para$Rho <- Rho
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
  Para$U <- U
  Para$LStarJ <- LStarJ
  Para$SpCov <- SpCov
  Para$CholSpCov <- CholSpCov
  Para$SpCovInv <- SpCovInv
  Para$CholKappa <- CholKappa
  Para$KappaInv <- KappaInv
  return(Para)

}



###Function that creates the data augmentation (i.e. Tobit) booleans------------------------------------------------------------
CreateDatAug <- function(DatObj) {

  ###Set data object
  YObserved <- DatObj$YObserved
  FamilyInd <- DatObj$FamilyInd
  YStarWide <- DatObj$YStarWide
  M <- DatObj$M
  O <- DatObj$O
  Nu <- DatObj$Nu
  N <- DatObj$N

  ###Initialize Data Augmentation Object with Normality
  DatAug <- list()
  DatAug$NBelow <- 0
  DatAug$NAbove <- 0
  DatAug$WhichBelow <- 0
  DatAug$WhichAbove <- 0
  
  ###Lower truncation
  if (any(FamilyInd == 2) | any(FamilyInd == 1)) {
    TobitBoolean <- matrix(FALSE, ncol = 1, nrow = N)
    for (o in 1:O) {
      if (FamilyInd[o] == 2 | FamilyInd[o] == 1) {
        TobitBooleanCube <- array(FALSE, dim = c(M, O, Nu))
        TobitBooleanCube[ , o, ] <- TRUE
        TobitBoolean <- TobitBoolean | matrix(TobitBooleanCube & YObserved <= 0, ncol = 1)
      }
    }
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    DatAug$NBelow <- NBelow
    DatAug$WhichBelow <- WhichBelow - 1
  }
  
  ###Upper truncation
  if (any(FamilyInd == 1)) {
    ProbitBoolean <- matrix(FALSE, ncol = 1, nrow = N)
    for (o in 1:O) {
      if (FamilyInd[o] == 1) {
        ProbitBoolean <- array(FALSE, dim = c(M, O, Nu))
        ProbitBoolean[ , o, ] <- TRUE
        ProbitBoolean <- ProbitBoolean | matrix(ProbitBoolean & YObserved > 0, ncol = 1)
      }
    }
    WhichAbove <- which(ProbitBoolean)
    NAbove <- length(WhichAbove)
    DatAug$WhichAbove <- WhichAbove - 1
    DatAug$NAbove <- NAbove
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
  O <- DatObj$O

  ###Set MCMC objects
  NKeep <- McmcObj$NKeep

  ###Create storage object 
  # Lambda: M * O x K
  # Eta: K x Nu
  # Sigma: M
  # Kappa: O x (O + 1) / 2
  # Rho: 1
  # Delta: K
  # Upsilon: K x (K + 1) / 2
  # Psi: 1
  # Xi: M * O * K
  Out <- matrix(nrow = (M * O * K + K * Nu + M + ((O + 1) * O) / 2  + K + ((K + 1) * K) / 2 + 1 + M * O * K + 1), ncol = NKeep)
  return(Out)

}