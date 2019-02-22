###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  C <- DatObj$C
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(O * M * K)]
  Eta <- RawSamples[, (O * M * K + 1):(O * M * K + K * Nu)]
  if (C == O) Sigma2 <- NULL
  if (C != O) Sigma2 <- RawSamples[, (O * M * K + K * Nu + 1):(O * M * K + K * Nu + M * (O - C))]
  Kappa <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2)]
  Delta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K)]
  Tau <- t(apply(Delta, 1, cumprod))
  Upsilon <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2)]
  Psi <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1), drop = FALSE]
  Xi <- RawSamples[, (O * M * K + K * Nu + M  * (O - C)+ (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + O * M * K)]
  Rho <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + O * M * K + 1), drop = FALSE]
  LambdaInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Lambda) <- paste0("Lambda_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  colnames(Eta) <- as.character((apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Eta", 1:Nu, "_"), x))))
  if (C != O) Sigma2Ind <- expand.grid(which(DatObj$FamilyInd != 3), 1:M)
  if (C != O) colnames(Sigma2) <- paste0("Sigma2_", Sigma2Ind[, 1], "_", Sigma2Ind[, 2])
  KappaInd <- which(lower.tri(apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x)), diag = TRUE), arr.ind = TRUE)
  colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]]
  colnames(Delta) <- paste0("Delta", 1:K)
  colnames(Tau) <- paste0("Tau", 1:K)
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  colnames(Psi) <- "Psi"
  colnames(Xi) <- paste0("Xi_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  colnames(Rho) <- "Rho"
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa = Kappa, Delta = Delta, Tau = Tau, Upsilon = Upsilon, Psi = Psi, Xi = Xi, Rho = Rho)
  return(Out)
}



###Function for creating a data object that contains objects needed for ModelFit-----------------------------------
OutputDatObj <- function(DatObj) {

  ###Collect needed objects
  # DatObjOut <- list(M = DatObj$M,
  #                   Nu = DatObj$Nu,
  #                   AdjacentEdgesBoolean = DatObj$AdjacentEdgesBoolean,
  #                   W = DatObj$W,
  #                   EyeM = DatObj$EyeM,
  #                   EyeNu = DatObj$EyeNu,
  #                   OneM = DatObj$OneM,
  #                   OneN = DatObj$OneN,
  #                   OneNu = DatObj$OneNu,
  #                   YStarWide = DatObj$YStarWide,
  #                   Rho = DatObj$Rho,
  #                   FamilyInd = DatObj$FamilyInd,
  #                   ScaleY = DatObj$ScaleY,
  #                   YObserved = DatObj$YObserved,
  #                   ScaleDM = DatObj$ScaleDM,
  #                   Time = DatObj$Time,
  #                   TimeVec = DatObj$TimeVec,
  #                   YObserved = DatObj$YObserved,
  #                   tNu = DatObj$tNu,
  #                   t1 = DatObj$t1,
  #                   XThetaInd = DatObj$XThetaInd,
  #                   N = DatObj$N,
  #                   EyeN = DatObj$EyeN)
  DatObjOut <- DatObj
  return(DatObjOut)

}



###Function for creating a data augmentation object that contains objects needed for ModelFit----------------------
OutputDatAug <- function(DatAug) {

  ###Collect needed objects
  # DatAugOut <- list(NBelow = DatAug$NBelow,
  #                   NBelowList = DatAug$NBelowList,
  #                   TobitBooleanMat = DatAug$TobitBooleanMat,
  #                   YStarNonZeroList = DatAug$YStarNonZero)
  DatAugOut <- DatAug
  return(DatAugOut)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  SpCorInd <- DatObj$SpCorInd

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropPsi <- MetropRcpp$MetropPsi
  AcceptancePsi <- MetropRcpp$AcceptancePsi
  OriginalTuners <- MetrObj$OriginalTuners[1]
  AcceptancePct <- AcceptancePsi / NSims
  MetrSummary <- cbind(AcceptancePct, MetropPsi, OriginalTuners)
  rownames(MetrSummary) <- "Psi"
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  
  ###Add Rho
  if (SpCorInd == 0) {
    MetropRho <- MetropRcpp$MetropRho
    AcceptanceRho <- MetropRcpp$AcceptanceRho
    OriginalTuners <- MetrObj$OriginalTuners
    AcceptancePct <- c(AcceptancePsi, AcceptanceRho) / NSims
    MetrSummary <- cbind(AcceptancePct, c(MetropPsi, MetropRho), OriginalTuners)
    rownames(MetrSummary) <- c("Psi", "Rho")
    colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  }

  ###Summarize and output
  return(MetrSummary)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.spBFA
#'
#' \code{is.spBFA} is a general test of an object being interpretable as a
#' \code{\link{spBFA}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{spBFA}} class is defined as the regression object that
#'  results from the \code{\link{spBFA}} regression function.
#'
#' @export
is.spBFA <- function(x) {
  identical(attributes(x)$class, "spBFA")
}