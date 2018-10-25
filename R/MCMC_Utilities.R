###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(M * K)]
  Eta <- RawSamples[, (M * K + 1):(M * K + K * Nu)]
  Sigma2 <- RawSamples[, (M * K + K * Nu + 1):(M * K + K * Nu + M)]
  Kappa2 <- RawSamples[, M * K + K * Nu + M + 1, drop = FALSE]
  Delta <- RawSamples[, (M * K + K * Nu + M + 1 + 1):(M * K + K * Nu + M + 1 + K)]
  Tau <- t(apply(Delta, 1, cumprod))
  Upsilon <- RawSamples[, (M * K + K * Nu + M + 1 + K + 1):(M * K + K * Nu + M + 1 + K + (K * (K + 1)) / 2)]
  Psi <- RawSamples[, (M * K + K * Nu + M + 1 + K + (K * (K + 1)) / 2 + 1), drop = FALSE]
  colnames(Lambda) <- as.character(t(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Lambda", 1:M, "_"), x))))
  colnames(Eta) <- as.character((apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Eta", 1:Nu, "_"), x))))
  colnames(Sigma2) <- paste0("Sigma2_", 1:M)
  colnames(Kappa2) <- "Kappa2"
  colnames(Delta) <- paste0("Delta", 1:K)
  colnames(Tau) <- paste0("Tau", 1:K)
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Lambda", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Lambda", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  colnames(Psi) <- "Psi"
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa2 = Kappa2, Delta = Delta, Tau = Tau, Upsilon = Upsilon, Psi = Psi)
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
  M <- DatObj$M

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropPsi <- MetropRcpp$MetropPsi
  AcceptancePsi <- MetropRcpp$AcceptancePsi
  OriginalTuners <- MetrObj$OriginalTuners

  ###Summarize and output
  AcceptancePct <- AcceptancePsi / NSims
  MetrSummary <- cbind(AcceptancePct, MetropPsi, OriginalTuners)
  rownames(MetrSummary) <- "Psi"
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
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