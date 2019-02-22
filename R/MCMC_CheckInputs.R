CheckInputs <- function(Y, Dist, Time, K, L, Trials, Starting, Hypers, Tuning, MCMC, Family, TemporalStructure, SpatialStructure) {
  
  ###Data dimensions
  N <- length(as.numeric(Y))
  M <- dim(Y)[1]
  O <- dim(Y)[2]
  Nu <- dim(Y)[3]

  ###Family
  if ((length(Family) != O) & (length(Family) != 1)) stop(paste0('Family: must have 1 or O = ', O, ' entries'))
  bool <- logical(length = length(Family))
  for (i in 1:length(Family)) bool[i] <- !(Family[i] %in% c("normal", "probit", "tobit", "binomial"))
  if (any(bool)) stop('Family: All entries must be one of "normal", "probit", "tobit", or "binomial"')
  
  ###Temporal correlation
  if (!TemporalStructure %in% c("ar1", "exponential")) stop('TemporalStructure: must be one of "ar1" or "exponential"')

  ###Spatial correlation
  if (!SpatialStructure %in% c("discrete", "continuous")) stop('SpatialStructure: must be one of "discrete" or "continuous"')
  
  ###Data checks for Trials
  if (is.null(Trials) & any(Family %in% "binomial")) stop('Trials must be specified for Family "binomial"')
  if (!is.null(Trials)) {
    if (!any(Family %in% "binomial")) stop('Trials can not be specified without Family = "binomial"')
    if (!is.array(Trials)) stop('Trials must be an array')
    # if (length(Trials) != N) stop(paste0('Trials must have exactly ', N, 'values')) # I check dimensions later
    if (any(is.na(Trials))) stop("Trials may have no missing values")
    if (any(!is.finite(Trials))) stop("Trials must have strictly finite entries")
    if (!isTRUE(all(Trials == floor(Trials)))) stop("Trials must have integers only")
    if (any(Trials < 1)) stop("Trials must contain positive integers only")
  }
  
  ###Data checks for Y
  if (!is.array(Y)) stop('Y must be an array')
  if (length(Y) != N) stop(paste0('Y must have exactly ', N, 'values'))
  if (any(is.na(Y))) stop("Y may have no missing values")
  if (any(!is.finite(Y))) stop("Y must have strictly finite entries")
  if (any(Family == "probit")) {
    for (o in 1:O) {
      if (Family[o] == "probit") {
        if ((sum(Y[ , o, ] == 1) + sum(Y[ , o, ] == 0)) != (Nu * M)) stop('Y: for "probit" observed data must be binary')
      }
    }
  } 
  if (any(Family == "tobit")) {
    for (o in 1:O) {
      if (Family[o] == "tobit") {
        if (any(Y[ , o, ] < 0)) stop('Y: for "tobit" observed data must be non-negative')
      }
    }
  } 
  if (any(Family == "binomial")) {
    count <- 1
    for (o in 1:O) {
      if (Family[o] == "binomial") {
        if (any(Y[ , o, ] < 0)) stop('Y: for "binomial" observed data must be non-negative')
        if (!isTRUE(all(Y[, o, ] == floor(Y[ , o, ])))) stop('Y: for "binomial" observed data must be non-negative integers')
        if (any(Y[, o, ] > Trials[, count, ])) stop('Y: for "binomial" observed data must be less than the corresponding number of Trials')
      }
    }
  } 
  
  ###Data checks for Dist
  if (!is.matrix(Dist)) stop('Dist must be a matrix')
  if (!dim(Dist)[1] == M) stop(paste0('Dist must be a ', M ,' x ', M, ' dimensional matrix'))
  if (!dim(Dist)[2] == M) stop('Dist must be square')
  if (sum(!((Dist) == t(Dist))) > 0) stop('Dist must be symmetric')
  if (any(diag(Dist) != 0)) stop('Dist must have only zeros on the diagonal')
  if (!all(!is.na(Dist))) stop('Dist cannot have missing values')
  if (!all(is.finite(Dist))) stop('Dist cannot have infinite values')
  if (SpatialStructure == "discrete") if (length(table(Dist)) > 2) stop('Dist must only contain binaries (i.e. 0\'s or 1\'s) for "discrete" space')
  
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
  if (!is.scalar(L) & !is.infinite(L)) stop('L must be a scalar or Inf')
  if (is.na(L)) stop('L cannot be NA')
  if (is.finite(L)) if (!is.wholenumber(L) | L <= 0) stop('L as a scalar must be a strictly positive integer')
  
  ###Hypers
  if (!is.null(Hypers)) {
    if (!is.list(Hypers)) stop('Hypers must be a list')
    if (!all(names(Hypers) %in% c("Sigma2", "Kappa", "Rho", "Delta", "Psi", "Upsilon"))) stop('Hypers: Can only contain lists with names "Sigma2", "Kappa", "Rho", "Delta", "Psi", or "Upsilon"')

    ###If Sigma2 hyperparameters are provided
    if ("Sigma2" %in% names(Hypers)) {
      if (any(Family %in% c("normal", "probit", "tobit"))) {
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
      } else stop('Hypers: "Sigma2" cannot be included for "binomial" likelihood')
    }
    
    ###If Kappa hyperparameters are provided
    if ("Kappa" %in% names(Hypers)) {
      if (!is.list(Hypers$Kappa)) stop('Hypers: "Kappa" must be a list')
      if (!"SmallUpsilon" %in% names(Hypers$Kappa)) stop('Hypers: "Kappa" value missing for SmallUpsilon')
      if (!is.scalar(Hypers$Kappa$SmallUpsilon)) stop('Hypers: "SmallUpsilon" must be a scalar')
      if (is.na(Hypers$Kappa$SmallUpsilon)) stop('Hypers: "SmallUpsilon" cannot be NA')
      if (!is.finite(Hypers$Kappa$SmallUpsilon)) stop('Hypers: "Kappa" cannot be infinite')
      if (Hypers$Kappa$SmallUpsilon <= 0) stop('Hypers: "SmallUpsilon" must be strictly positive')
      if (!"BigTheta" %in% names(Hypers$Kappa)) stop('Hypers: "BigTheta" value missing')
      if (!is.matrix(Hypers$Kappa$BigTheta)) stop('Hypers: "BigTheta" must be a matrix')
      if (!dim(Hypers$Kappa$BigTheta)[1] == O) stop('Hypers: "BigTheta" must be O dimensional')
      if (!all(!is.na(Hypers$Kappa$BigTheta))) stop('Hypers: "BigTheta" cannot have missing values')
      if (!all(is.finite(Hypers$Kappa$BigTheta))) stop('Hypers: "BigTheta" cannot have infinite values')
      if (!dim(Hypers$Kappa$BigTheta)[2] == O) stop('Hypers: "BigTheta" must be square')
      if (sum( !( (Hypers$Kappa$BigTheta) == t(Hypers$Kappa$BigTheta) ) ) > 0) stop('Hypers: "BigTheta" must be symmetric')
      if ((det(Hypers$Kappa$BigTheta) - 0) < 0.00001) stop('Hypers: "BigTheta" is close to singular')
      
    }

    ###If Rho hyperparameters are provided
    if ("Rho" %in% names(Hypers)) {
      if (SpatialStructure == "discrete") if (!is.null(Hypers$Rho)) stop('Hypers: When SpatialStructure = "discrete", "Rho" must be missing or NULL')
      if (SpatialStructure == "continuous") {
        if (!is.list(Hypers$Rho)) stop('Hypers: "Rho" must be a list')
        if (!"ARho" %in% names(Hypers$Rho)) stop('Hypers: "ARho" value missing')
        if (!is.scalar(Hypers$Rho$ARho)) stop('Hypers: "ARho" must be a scalar')
        if (is.na(Hypers$Rho$ARho)) stop('Hypers: "ARho" cannot be NA')
        if (!is.finite(Hypers$Rho$ARho)) stop('Hypers: "ARho" cannot be infinite')
        if (Hypers$Rho$ARho <= 0) stop('Hypers: "ARho" must be strictly positive')
        if (!"BRho" %in% names(Hypers$Rho)) stop('Hypers: "BRho" value missing')
        if (!is.scalar(Hypers$Rho$BRho)) stop('Hypers: "BRho" must be a scalar')
        if (is.na(Hypers$Rho$BRho)) stop('Hypers: "BRho" cannot be NA')
        if (!is.finite(Hypers$Rho$BRho)) stop('Hypers: "BRho" cannot be infinite')
        if (Hypers$Rho$BRho <= 0) stop('Hypers: "BRho" must be strictly positive')
        if (Hypers$Rho$BRho < Hypers$Rho$ARho) stop('Hypers: "BRho" must be greater than "ARho"')
      }
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
    if (!all(names(Starting) %in% c("Sigma2", "Kappa", "Rho", "Delta", "Psi", "Upsilon"))) stop('Starting: Can only contain objects with names "Sigma2", "Kappa", "Rho", "Delta", "Psi", and "Upsilon"')

    ###If Delta starting values is provided
    if ("Delta" %in% names(Starting)) {
      if (is.vector(Starting$Delta)) {
        if (!is.numeric(Starting$Delta)) stop('Starting: "Delta" must be a vector')
        if (length(Starting$Delta) != K) stop('Starting: "Delta" must be length K')
        if (!all(!is.na(Starting$Delta))) stop('Starting: "Delta" cannot have missing values')
        if (!all(is.finite(Starting$Delta))) stop('Starting: "Delta" cannot have infinite values')
        if (any(Starting$Delta < 0)) stop('Starting: "Delta" cannot have non-negative values')
      }
    }
    if (is.scalar(Starting$Delta)) {
      if (is.na(Starting$Delta)) stop('Starting: "Delta" cannot be NA')
      if (!is.finite(Starting$Delta)) stop('Starting: "Delta" cannot be infinite')
      if (Starting$Delta < 0) stop('Starting: "Delta" must be non-negative')
    }
    if ((!is.vector(Starting$Delta)) & (!is.scalar(Starting$Delta))) stop('Starting: "Delta" must be a scalar or a vector')
    
    ###If Sigma2 starting values is provided
    if ("Sigma2" %in% names(Starting)) {
      if (any(Family %in% c("normal", "probit", "tobit"))) {
        if (is.matrix(Starting$Sigma2)) {
          if (!dim(Starting$Sigma2)[1] == M) stop(paste0('Starting: "Sigma2" must have', M, 'rows'))
          if (!dim(Starting$Sigma2)[2] == (O - C)) stop(paste0('Starting: "Sigma2" must have', O - C, 'columns'))
          if (!all(!is.na(Starting$Sigma2))) stop('Starting: "Sigma2" cannot have missing values')
          if (!all(is.finite(Starting$Sigma2))) stop('Starting: "Sigma2" cannot have infinite values')
          if (Starting$Sigma2 <= 0) stop('Starting: "Sigma2" must be strictly positive')
        }
        if (is.scalar(Starting$Sigma2)) {
            if (is.na(Starting$Delta)) stop('Starting: "Sigma2" cannot be NA')
            if (!is.finite(Sigma2$Sigma2)) stop('Starting: "Sigma2" cannot be infinite')
            if (Starting$Sigma2 < 0) stop('Starting: "Sigma2" must be non-negative')
          }
        if ((!is.matrix(Starting$Sigma2)) & (!is.scalar(Starting$Sigma2))) stop('Starting: "Sigma2" must be a scalar or a matrix')
      } else stop('Starting: "Sigma2" does not get included for "binomial" likelihood')
    }
    
    ###If Kappa starting values is provided
    if ("Kappa" %in% names(Starting)) {
      if (!is.matrix(Starting$Kappa)) stop('Starting: "Kappa" must be a matrix')
      if (!dim(Starting$Kappa)[1] == O) stop('Starting: "Kappa" must be O dimensional')
      if (!dim(Starting$Kappa)[2] == O) stop('Starting: "Kappa" must be square')
      if (!all(!is.na(Starting$Kappa))) stop('Starting: "Kappa" cannot have missing values')
      if (!all(is.finite(Starting$Kappa))) stop('Starting: "Kappa" cannot have infinite values')
      if (sum( !( (Starting$Kappa) == t(Starting$Kappa) ) ) > 0) stop('Starting: "Kappa" must be symmetric')
      if ((det(Starting$Kappa) - 0) < 0.0000000001) stop('Starting: "Kappa" is close to singular')
    }

    ###If Rho starting values is provided
    if ("Rho" %in% names(Starting)) {
      if (!is.scalar(Starting$Rho)) stop('Starting: "Rho" must be a scalar')
      if (is.na(Starting$Rho)) stop('Starting: "Rho" cannot be NA')
      if (!is.finite(Starting$Rho)) stop('Starting: "Rho" cannot be infinite')
      if (SpatialStructure == "discrete") {
        if ((Starting$Rho <= 0) | (Starting$Rho >= 1)) stop('Starting: "Rho" must be in (0, 1) for discrete spatial process')
      }
      if (SpatialStructure == "continuous") {
        # I make sure that Rho is in (ARho, BRho) in CreatePara();
      }
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
    if (!all(names(Tuning) %in% c("Psi", "Rho"))) stop('Tuning: Can only contain objects with names "Psi", "Rho"')

    ###If Psi tuning value is provided
    if ("Psi" %in% names(Tuning)) {
      if (!is.scalar(Tuning$Psi)) stop('Tuning: "Psi" must be a scalar')
      if (is.na(Tuning$Psi)) stop('Tuning: "Psi" cannot be NA')
      if (!is.finite(Tuning$Psi)) stop('Tuning: "Psi" cannot be infinite')
      if (Tuning$Psi < 0) stop('Tuning: "Psi" must be non-negative')
    }

    ###If Rho tuning value is provided
    if ("Rho" %in% names(Tuning)) {
      if (SpatialStructure == "discrete") {
        if (!is.null(Tuning$Rho)) stop('Tuning: No tuning is needed for "Rho" when "SpatialStructure" is discrete')
      }
      if (SpatialStructure == "continuous") {
        if (!is.scalar(Tuning$Rho)) stop('Tuning: "Rho" must be a scalar')
        if (is.na(Tuning$Rho)) stop('Tuning: "Rho" cannot be NA')
        if (!is.finite(Tuning$Rho)) stop('Tuning: "Rho" cannot be infinite')
        if (Tuning$Rho < 0) stop('Tuning: "Rho" must be non-negative')
      }
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
