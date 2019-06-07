###Function to get model fit diagnostics given a spBFA object
#'
#' predict.spBFA
#'
#' Predicts future observations from the \code{\link{spBFA}} model.
#'
#' @param object a \code{\link{spBFA}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.spBFA} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged observed outcomes along with the
#'  observational level parameters (\code{mu}, \code{tau}, and \code{alpha}).
#'
#' @return \code{predict.spBFA} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{MuTauAlpha}}{A \code{list} containing three matrices, \code{mu},
#'   \code{tau} and \code{alpha}. Each matrix is dimension \code{NKeep x s}, where
#'   \code{s} is the number of new time points. Each matrix contains posterior
#'   samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{s} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x s}, where \code{s}
#'   is the number of new time points. Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spBFA function
predict.spBFA <- function(object, NewTimes, NewTrials = NULL, type = "temporal", ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.spBFA(object)) stop('"object" must be of class spBFA')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')
  if (!is.character(type)) stop('"type" must be a character string')
  if (!(type %in% c("temporal", "spatial"))) stop('"type" must be one of "spatial" or "temporal"')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O

  ###Create updated distance matrix
  TimeFixed <- DatObj$Time
  Time <- sort(c(TimeFixed, NewTimes))
  TimeDist <- abs(outer(Time, Time, "-" ))
  NNewVisits <- length(NewTimes)
  NewVisits <- OriginalVisits <- NULL
  for (i in 1:NNewVisits) NewVisits <- c(NewVisits, which(NewTimes[i] == Time) - 1)
  for (i in 1:Nu) OriginalVisits <- c(OriginalVisits, which(TimeFixed[i] == Time) - 1)

  ###Update DatObj
  DatObj$NewVisits <- NewVisits
  DatObj$OriginalVisits <- OriginalVisits
  DatObj$TimeDist <- TimeDist
  DatObj$NNewVisits <- NNewVisits
  DatObj$EyeK <- diag(DatObj$K)
  
  ###Create Trials object
  if (DatObj$C == 0) {
    DatObj$Trials <- array(0, dim = c(M, O, Nu))
  }
  if ((DatObj$C > 0) & is.null(NewTrials)) {
    Trials <- array(dim = c(M, DatObj$C, NNewVisits))
    for (n in 1:NNewVisits) Trials[, , n] <- DatObj$Trials[, , Nu]
  }
  if ((DatObj$C > 0) & !is.null(NewTrials)) {
    Trials <- NewTrials
    if (!is.array(Trials)) stop('Trials must be an array')
    if (dim(Trials)[1] != M) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (dim(Trials)[1] != C) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (dim(Trials)[1] != NNewVisits) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (any(is.na(Trials))) stop("Trials may have no missing values")
    if (any(!is.finite(Trials))) stop("Trials must have strictly finite entries")
    if (!isTRUE(all(Trials == floor(Trials)))) stop("Trials must have integers only")
    if (any(Trials < 1)) stop("Trials must contain positive integers only")
  }
  
  ###Set mcmc object
  NKeep <- dim(object$rho)[1]

  ###Create parameter object
  Para <- list()
  Para$Psi <- object$psi
  Para$Upsilon <- object$upsilon
  Para$Lambda <- object$lambda
  Para$Eta <- object$eta
  if (is.null(object$sigma2)) Para$Sigma2 <- matrix(1)
  if (!is.null(object$sigma2)) Para$Sigma2 <- object$sigma2
  
  ###Obtain samples of eta using Bayesian krigging
  EtaKrig <- EtaKrigging(DatObj, Para, NKeep)

  ###Obtain samples of observed Y
  YKrig <- YKrigging(DatObj, Para, EtaKrig, NKeep)

  ###Format theta samples for output
  EtaOut <- list()
  Eta <- array(t(EtaKrig), dim = c(K, NKeep, NNewVisits))
  for (n in 1:NNewVisits) EtaOut[[n]] <- t(Eta[, , n])
  for (n in 1:NNewVisits) colnames(EtaOut[[n]]) <- paste0("Eta", 1:K, "_", Nu + n)
  names(EtaOut) <- paste0("Eta", Nu + 1:NNewVisits)
  
  ###Format Y samples for output
  YOut <- list()
  YInd <- expand.grid(1:M, 1:O)
  for (n in 1:NNewVisits) YOut[[n]] <- t(YKrig[, n, ])
  for (n in 1:NNewVisits) colnames(YOut[[n]]) <- paste0("Y_", Nu + n, "_", YInd[, 2], "_", YInd[, 1])
  names(YOut) <- paste0("Y", Nu + 1:NNewVisits)  
    
  ###Return formated samples
  return(list(Eta = EtaOut, Y = YOut))

}
