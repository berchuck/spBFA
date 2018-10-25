CheckInputs <- function(Y, W, Time, K, L, Starting, Hypers, Tuning, MCMC, Family, TemporalStructure, Rho, ScaleY) {
  
  ###Data dimensions
  N <- length(Y)
  M <- dim(W)[1]
  Nu <- length(Time)

  ###Family
  if (!Family %in% c("normal", "probit", "tobit")) stop('Family: must be one of "normal", "probit" or "tobit"')

  ###Temporal correlation
  if (!TemporalStructure %in% c("ar1", "exponential")) stop('TemporalStructure: must be one of "ar1" or "exponential"')
  
  ###Rho
  if (missing(Rho)) stop("Rho: missing")
  if (!is.scalar(Rho)) stop('Rho must be a scalar')
  if (is.na(Rho)) stop('Rho cannot be NA')
  if (!is.finite(Rho)) stop('Rho cannot be infinite')
  if (!((Rho < 1) & (Rho > 0))) stop('Rho must be in (0, 1)')

  ###ScaleY
  if (missing(ScaleY)) stop("ScaleY: missing")
  if (!is.scalar(ScaleY)) stop('ScaleY must be a scalar')
  if (is.na(ScaleY)) stop('ScaleY cannot be NA')
  if (!is.finite(ScaleY)) stop('ScaleY cannot be infinite')
  if (!(ScaleY > 0)) stop('ScaleY must be positive')

  ###Data checks for Y
  if (!is.numeric(Y)) stop('Y must be a vector')
  if (length(Y) != N) stop(paste0('Y must have length ', N))
  if (any(is.na(Y))) stop("Y may have no missing values")
  if (any(!is.finite(Y))) stop("Y must have strictly finite entries")
  if ((Family == "probit") & ((sum(Y == 1) + sum(Y == 0)) != N)) stop('Y: for "probit" observed data must be binary')
  if ((Family == "tobit") & (any(Y < 0))) stop('Y: for "tobit" observed data must be non-negative')

  ###Data checks for W
  if (!is.matrix(W)) stop('W must be a matrix')
  if (!dim(W)[1] == M) stop(paste0('W must be a ', M ,' x ', M, ' dimensional matrix'))
  if (!dim(W)[2] == M) stop('W must be square')
  if (sum(!((W) == t(W))) > 0) stop('W must be symmetric')
  if (length(table(W)) > 2) stop('W must only contain binaries (i.e. 0\'s or 1\'s)')
  if (any(diag(W) != 0)) stop('W must have only zeros on the diagonal')

  ###Data checks for Time
  if (!is.numeric(Time)) stop('Time must be a vector')
  if (length(Time) != Nu) stop(paste0('Time must have length ', Nu))
  if (any(is.na(Time))) stop("Time may have no missing values")
  if (any(!is.finite(Time))) stop("Time must have strictly finite entries")
  if (is.unsorted(Time)) stop('Time vector is not in increasing order')
  if (!all(Time >= 0)) stop('Time vector has at least one negative point')

  ###Data checks for K
  if (missing(K)) stop("K: missing")
  if (!is.scalar(K)) stop('K must be a scalar')
  if (is.na(K)) stop('K cannot be NA')
  if (!is.finite(K)) stop('K cannot be infinite')
  if (!is.wholenumber(K) | K <= 0) stop('K must be a strictly positive integer')

  ###Data checks for L
  if (missing(L)) stop("L: missing")
  if (!is.scalar(L)) stop('L must be a scalar')
  if (is.na(L)) stop('L cannot be NA')
  if (!is.finite(L)) stop('L cannot be infinite')
  if (!is.wholenumber(L) | L <= 0) stop('L must be a strictly positive integer')
  
  ###Verify dimensions
  if ((N / M) != Nu) stop('Time and Y have contradictory dimensions')

  ###Hypers
  if (!is.null(Hypers)) {
    if (!is.list(Hypers)) stop('Hypers must be a list')
    if (!all(names(Hypers) %in% c("Sigma2", "Kappa2", "Delta", "Psi", "Upsilon"))) stop('Hypers: Can only contain lists with names "Sigma2", "Kappa2" and "Delta", "Psi", or "Upsilon"')

    ###If Sigma2 hyperparameters are provided
    if ("Sigma2" %in% names(Hypers)) {
      if (!is.list(Hypers$Sigma2)) stop('Hypers: "Sigma2" must be a list')
      if (!"A" %in% names(Hypers$Sigma2)) stop('Hypers: "A" value missing')
      if (!is.scalar(Hypers$Sigma2$A)) stop('Hypers: "A" must be a scalar')
      if (is.na(Hypers$Sigma2$A)) stop('Hypers: "A" cannot be NA')
      if (!is.finite(Hypers$Sigma2$A)) stop('Hypers: "A" cannot be infinite')
      if (Hypers$Sigma2$A <= 0) stop('Hypers: "A" must be strictly positive')
      if (!"B" %in% names(Hypers$Sigma2)) stop('Hypers: "B" value missing')
      if (!is.scalar(Hypers$Sigma2$B)) stop('Hypers: "B" must be a scalar')
      if (is.na(Hypers$Sigma2$B)) stop('Hypers: "B" cannot be NA')
      if (!is.finite(Hypers$Sigma2$B)) stop('Hypers: "B" cannot be infinite')
      if (Hypers$Sigma2$B <= 0) stop('Hypers: "B" must be strictly positive')
    }
    
    ###If Kappa2 hyperparameters are provided
    if ("Kappa2" %in% names(Hypers)) {
      if (!is.list(Hypers$Kappa2)) stop('Hypers: "Kappa2" must be a list')
      if (!"C" %in% names(Hypers$Kappa2)) stop('Hypers: "C" value missing')
      if (!is.scalar(Hypers$Kappa2$C)) stop('Hypers: "C" must be a scalar')
      if (is.na(Hypers$Kappa2$C)) stop('Hypers: "C" cannot be NA')
      if (!is.finite(Hypers$Kappa2$C)) stop('Hypers: "C" cannot be infinite')
      if (Hypers$Kappa2$C <= 0) stop('Hypers: "C" must be strictly positive')
      if (!"D" %in% names(Hypers$Kappa2)) stop('Hypers: "D" value missing')
      if (!is.scalar(Hypers$Kappa2$D)) stop('Hypers: "D" must be a scalar')
      if (is.na(Hypers$Kappa2$D)) stop('Hypers: "D" cannot be NA')
      if (!is.finite(Hypers$Kappa2$D)) stop('Hypers: "D" cannot be infinite')
      if (Hypers$Kappa2$D <= 0) stop('Hypers: "D" must be strictly positive')
    }
    
    ###If Delta hyperparameters are provided
    if ("Delta" %in% names(Hypers)) {
      if (!is.list(Hypers$Delta)) stop('Hypers: "Delta" must be a list')
      if (!"A1" %in% names(Hypers$Delta)) stop('Hypers: "A1" value missing')
      if (!is.scalar(Hypers$Delta$A1)) stop('Hypers: "A1" must be a scalar')
      if (is.na(Hypers$Delta$A1)) stop('Hypers: "A1" cannot be NA')
      if (!is.finite(Hypers$Delta$A1)) stop('Hypers: "A1" cannot be infinite')
      if (Hypers$Delta$A1 <= 0) stop('Hypers: "A1" must be strictly positive')
      if (!"A2" %in% names(Hypers$Delta)) stop('Hypers: "A2" value missing')
      if (!is.scalar(Hypers$Delta$A2)) stop('Hypers: "A2" must be a scalar')
      if (is.na(Hypers$Delta$A2)) stop('Hypers: "A2" cannot be NA')
      if (!is.finite(Hypers$Delta$A2)) stop('Hypers: "A2" cannot be infinite')
      if (Hypers$Delta$A2 <= 0) stop('Hypers: "A2" must be strictly positive')
    }
    
    ###If Upsilon hyperparameters are provided
    if ("Upsilon" %in% names(Hypers)) {
      if (!is.list(Hypers$Upsilon)) stop('Hypers: "Upsilon" must be a list')
      if (!"Zeta" %in% names(Hypers$Upsilon)) stop('Hypers: "Zeta" value missing')
      if (!is.scalar(Hypers$Upsilon$Zeta)) stop('Hypers: "Zeta" must be a scalar')
      if (is.na(Hypers$Upsilon$Zeta)) stop('Hypers: "Zeta" cannot be NA')
      if (!is.finite(Hypers$Upsilon$Zeta)) stop('Hypers: "Zeta" cannot be infinite')
      if (Hypers$Upsilon$Zeta < K) stop('Hypers: "Zeta" must be greater than or equal to K')
      if (!"Omega" %in% names(Hypers$Upsilon)) stop('Hypers: "Omega" value missing')
      if (!is.matrix(Hypers$Upsilon$Omega)) stop('Hypers: "Omega" must be a matrix')
      if (!dim(Hypers$Upsilon$Omega)[1] == K) stop('Hypers: "Omega" must be K dimensional')
      if (!all(!is.na(Hypers$Upsilon$Omega))) stop('Hypers: "Omega" cannot have missing values')
      if (!all(is.finite(Hypers$Upsilon$Omega))) stop('Hypers: "Omega" cannot have infinite values')
      if (!dim(Hypers$Upsilon$Omega)[2] == K) stop('Hypers: "Omega" must be square')
      if (sum( !( (Hypers$Upsilon$Omega) == t(Hypers$Upsilon$Omega) ) ) > 0) stop('Hypers: "Omega" must be symmetric')
      if ((det(Hypers$Upsilon$Omega) - 0) < 0.00001) stop('Hypers: "Omega" is close to singular')
    }

    ###If Psi hyperparameters are provided
    if ("Psi" %in% names(Hypers)) {
      if (!is.list(Hypers$Psi)) stop('Hypers: "Psi" must be a list')
      if (TemporalStructure == "exponential") {
        if (!"APsi" %in% names(Hypers$Psi)) stop('Hypers: "APsi" value missing')
        if (!is.scalar(Hypers$Psi$APsi)) stop('Hypers: "APsi" must be a scalar')
        if (is.na(Hypers$Psi$APsi)) stop('Hypers: "APsi" cannot be NA')
        if (!is.finite(Hypers$Psi$APsi)) stop('Hypers: "APsi" cannot be infinite')
        if (!"BPsi" %in% names(Hypers$Psi)) stop('Hypers: "BPsi" value missing')
        if (!is.scalar(Hypers$Psi$BPsi)) stop('Hypers: "BPsi" must be a scalar')
        if (is.na(Hypers$Psi$BPsi)) stop('Hypers: "BPsi" cannot be NA')
        if (!is.finite(Hypers$Psi$BPsi)) stop('Hypers: "BPsi" cannot be infinite')
        if (Hypers$Psi$APsi < 0) stop('Hypers: "APsi" must be non-negative')
        if (Hypers$Psi$BPsi <= 0) stop('Hypers: "BPsi" must be strictly positive')
        if (Hypers$Psi$BPsi < Hypers$Psi$APsi) stop('Hypers: "BPsi" must be greater than "APsi"')
      }
      if (TemporalStructure == "ar1") {
        if (!"APsi" %in% names(Hypers$Psi)) stop('Hypers: "APsi" value missing')
        if (!is.scalar(Hypers$Psi$APsi)) stop('Hypers: "APsi" must be a scalar')
        if (is.na(Hypers$Psi$APsi)) stop('Hypers: "APsi" cannot be NA')
        if (!is.finite(Hypers$Psi$APsi)) stop('Hypers: "APsi" cannot be infinite')
        if (!"BPsi" %in% names(Hypers$Psi)) stop('Hypers: "BPsi" value missing')
        if (!is.scalar(Hypers$Psi$BPsi)) stop('Hypers: "BPsi" must be a scalar')
        if (is.na(Hypers$Psi$BPsi)) stop('Hypers: "BPsi" cannot be NA')
        if (!is.finite(Hypers$Psi$BPsi)) stop('Hypers: "BPsi" cannot be infinite')
        if ((Hypers$Psi$APsi < -1) | (Hypers$Psi$APsi > 1)) stop('Hypers: "APsi" must be in the range (-1, 1)')
        if ((Hypers$Psi$BPsi < -1) | (Hypers$Psi$BPsi > 1)) stop('Hypers: "BPsi" must be in the range (-1, 1)')
        if (Hypers$Psi$BPsi < Hypers$Psi$APsi) stop('Hypers: "BPsi" must be greater than "APsi"')
        if (!"Beta" %in% names(Hypers$Psi)) stop('Hypers: "Beta" value missing')
        if (!is.scalar(Hypers$Psi$Beta)) stop('Hypers: "Beta" must be a scalar')
        if (is.na(Hypers$Psi$Beta)) stop('Hypers: "Beta" cannot be NA')
        if (!is.finite(Hypers$Psi$Beta)) stop('Hypers: "Beta" cannot be infinite')
        if (!"Gamma" %in% names(Hypers$Psi)) stop('Hypers: "Gamma" value missing')
        if (!is.scalar(Hypers$Psi$Gamma)) stop('Hypers: "Gamma" must be a scalar')
        if (is.na(Hypers$Psi$Gamma)) stop('Hypers: "Gamma" cannot be NA')
        if (!is.finite(Hypers$Psi$Gamma)) stop('Hypers: "Gamma" cannot be infinite')
        if (Hypers$Psi$Beta <= 0) stop('Hypers: "Beta" must be strictly positive')
        if (Hypers$Psi$Gamma <= 0) stop('Hypers: "Gamma" must be strictly positive')
      }
    }

  ###End Hyperparameters
  }

  ###Starting Values
  if (!is.null(Starting)) {
    if (!is.list(Starting)) stop('Starting must be a list')
    if (!all(names(Starting) %in% c("Sigma2", "Kappa2", "Delta", "Psi", "Upsilon"))) stop('Starting: Can only contain objects with names "Sigma2", "Kappa2", "Delta", "Psi", and "Upsilon"')

    ###If Delta starting values is provided
    if ("Delta" %in% names(Starting)) {
      if (!is.numeric(Starting$Delta)) stop('Starting: "Delta" must be a vector')
      if (length(Starting$Delta) != K) stop('Starting: "Delta" must be length K')
      if (!all(!is.na(Starting$Delta))) stop('Starting: "Delta" cannot have missing values')
      if (!all(is.finite(Starting$Delta))) stop('Starting: "Delta" cannot have infinite values')
    }

    ###If Sigma2 starting values is provided
    if ("Sigma2" %in% names(Starting)) {
      if (!is.scalar(Starting$Sigma2)) stop('Starting: "Sigma2" must be a scalar')
      if (is.na(Starting$Sigma2)) stop('Starting: "Sigma2" cannot be NA')
      if (!is.finite(Starting$Sigma2)) stop('Starting: "Sigma2" cannot be infinite')
      if (Starting$Sigma2 <= 0) stop('Starting: "Sigma2" must be strictly positive')
    }
    
    ###If Kappa2 starting values is provided
    if ("Kappa2" %in% names(Starting)) {
      if (!is.scalar(Starting$Kappa2)) stop('Starting: "Kappa2" must be a scalar')
      if (is.na(Starting$Kappa2)) stop('Starting: "Kappa2" cannot be NA')
      if (!is.finite(Starting$Kappa2)) stop('Starting: "Kappa2" cannot be infinite')
      if (Starting$Kappa2 <= 0) stop('Starting: "Kappa2" must be strictly positive')
    }

    ###If Upsilon starting values is provided
    if ("Upsilon" %in% names(Starting)) {
      if (!is.matrix(Starting$Upsilon)) stop('Starting: "Upsilon" must be a matrix')
      if (!dim(Starting$Upsilon)[1] == K) stop('Starting: "Upsilon" must be K dimensional')
      if (!dim(Starting$Upsilon)[2] == K) stop('Starting: "Upsilon" must be square')
      if (!all(!is.na(Starting$Upsilon))) stop('Starting: "Upsilon" cannot have missing values')
      if (!all(is.finite(Starting$Upsilon))) stop('Starting: "Upsilon" cannot have infinite values')
      if (sum( !( (Starting$Upsilon) == t(Starting$Upsilon) ) ) > 0) stop('Starting: "Upsilon" must be symmetric')
      if ((det(Starting$Upsilon) - 0) < 0.0000000001) stop('Starting: "Upsilon" is close to singular')
    }

    ###If Psi starting values is provided
    if ("Psi" %in% names(Starting)) {
      if (!is.scalar(Starting$Psi)) stop('Starting: "Psi" must be a scalar')
      if (is.na(Starting$Psi)) stop('Starting: "Psi" cannot be NA')
      if (!is.finite(Starting$Psi)) stop('Starting: "Psi" cannot be infinite')
      # I make sure that Psi is in (APsi, BPsi) in CreatePara();
    }

  ###End Starting Values
  }

  ###Tuning Values
  if (!is.null(Tuning)) {
    if (!is.list(Tuning)) stop('Tuning must be a list')
    if (!all(names(Tuning) %in% c("Psi"))) stop('Tuning: Can only contain objects with names "Psi"')

    ###If Psi tuning value is provided
    if ("Psi" %in% names(Tuning)) {
      if (!is.scalar(Tuning$Psi)) stop('Tuning: "Psi" must be a scalar')
      if (is.na(Tuning$Psi)) stop('Tuning: "Psi" cannot be NA')
      if (!is.finite(Tuning$Psi)) stop('Tuning: "Psi" cannot be infinite')
      if (Tuning$Psi < 0) stop('Tuning: "Psi" must be non-negative')
    }

  ###End Tuning Values
  }

  ###MCMC Values
  if (!is.null(MCMC)) {
    if (!is.list(MCMC)) stop('MCMC must be a list')
    if (!all(names(MCMC) %in% c("NBurn", "NSims", "NThin", "NPilot"))) stop('MCMC: Can only contain objects with names "NBurn", "NSims", "NThin" and "NPilot"')

    ###If NBurn is provided
    if ("NBurn" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NBurn)) stop('MCMC: "NBurn" must be a scalar')
      if (is.na(MCMC$NBurn)) stop('MCMC: "NBurn" cannot be NA')
      if (!is.finite(MCMC$NBurn)) stop('MCMC: "NBurn" cannot be infinite')
      if (!is.wholenumber(MCMC$NBurn) | MCMC$NBurn < 0) stop('MCMC: "NBurn" must be a non-negative integer')
      if (MCMC$NBurn < 100) stop('MCMC: "NBurn" must be at least 100')
    }

    ###If NSims is provided
    if ("NSims" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NSims)) stop('MCMC: "NSims" must be a scalar')
      if (is.na(MCMC$NSims)) stop('MCMC: "NSims" cannot be NA')
      if (!is.finite(MCMC$NSims)) stop('MCMC: "NSims" cannot be infinite')
      if (!is.wholenumber(MCMC$NSims) | MCMC$NSims <= 0) stop('MCMC: "NSims" must be a positive integer')
      if (MCMC$NSims < 100) stop('MCMC: "NSims" must be at least 100')
    }

    ###If NThin is provided
    if ("NThin" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NThin)) stop('MCMC: "NThin" must be a scalar')
      if (is.na(MCMC$NThin)) stop('MCMC: "NThin" cannot be NA')
      if (!is.finite(MCMC$NThin)) stop('MCMC: "NThin" cannot be infinite')
      if (!is.wholenumber(MCMC$NThin) | MCMC$NThin <= 0) stop('MCMC: "NThin" must be a positive integer')
      # if (!is.wholenumber(MCMC$NSims / MCMC$NThin)) stop('MCMC: "NThin" must be a factor of "NSims"') enforced in CreateMCMC();
    }

    ###If NPilot is provided
    if ("NPilot" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NPilot)) stop('MCMC: "NPilot" must be a scalar')
      if (is.na(MCMC$NPilot)) stop('MCMC: "NPilot" cannot be NA')
      if (!is.finite(MCMC$NPilot)) stop('MCMC: "NPilot" cannot be infinite')
      if (!is.wholenumber(MCMC$NPilot) | MCMC$NPilot < 0) stop('MCMC: "NPilot" must be a positive integer')
      # if (!is.wholenumber(MCMC$NBurn / MCMC$NPilot)) stop('MCMC: "NPilot" must be a factor of "NBurn"') enforced in CreateMCMC();
    }

  ###End MCMC Values
  }

}

###Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
