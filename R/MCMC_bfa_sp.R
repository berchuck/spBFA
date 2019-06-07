#' Spatial factor analysis using a Bayesian hierarchical model.
#'
#' \code{bfa_sp} is a Markov chain Monte Carlo (MCMC) sampler for a spatial factor analysis model. The spatial component is 
#' introduced using a Probit stick-breaking process prior on the factor loadings. The model is implemented using a Bayesian hierarchical framework.
#'
#' @param Y A \code{M x O x Nu} dimensional array containing the observed outcome data.
#'  Here, \code{M} represents the number of spatial locations, \code{O} the number of different observation types
#'  and \code{Nu} the number of temporal visits. The observations in \code{Y} must be first
#'  ordered spatially, second by observation type and then temporally. This means that each slice, Y[ , , t] is an \code{M x O} matrix
#'  with each column having observations from each spatial location \code{M} of a spatial observation type.
#'
#' @param Dist A \code{M x M} dimensional distance matrix. For a \code{discrete} spatial process the matrix contains binary adjacencies that dictate the
#'  spatial neigborhood structure and for \code{continuous} spatial processes the matrix should be a continuous distance matrix (e.g., Euclidean).
#'
#' @param Time A \code{Nu} dimensional vector containing the observed time points
#'  in increasing order.
#'
#' @param K A scalar that indicates the dimension (i.e., quantity) of latent factors.
#'  
#' @param L The number of latent clusters. If finite, a scalar indicating the number of clusters for each column of the factor loadings matrix. By default \code{L} is set at \code{Inf}
#'  so that the Probit stick-breaking process becomes an infinite mixture model.
#'  
#' @param Trials A \code{M x C x Nu} dimensional array containing the number of trials for each of the observations in Y. The second dimension 
#'  now is \code{C} dimensional, which corresponds to the number of \code{binomial} random variable types. If there is no count data, \code{Trials} should be left missing.
#'  
#' @param Starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon} or
#'  \code{Psi} containing appropriate objects. \code{Delta} must either be a \code{K} dimensional
#'  vector or a scalar (the scalar populates the entire vector). \code{Sigma2} must be either a \code{M x (O - C)} matrix or a scalar.
#'  \code{Kappa} must be a \code{O x O} dimensional matrix, \code{Rho} a scalar, \code{Upsilon} a \code{K x K} matrix, and \code{Psi} a scalar.
#'
#' @param Hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in (Berchuck et al.).
#'  Otherwise a \code{list} must be provided with names \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon} or
#'  \code{Psi} containing further hyperparameter information. These objects are themselves
#'  \code{lists} and may be constructed as follows.
#'
#'  \code{Delta} is a \code{list} with two objects, \code{A1} and \code{A2}. These values represent the prior shape 
#'  parameters for the multiplicative Gamma shrinkage prior.
#'  
#'  \code{Sigma2} is a \code{list} with two objects, \code{A} and \code{B}. These values represent the shape and scale for the variance parameters.
#'  
#'  \code{Kappa} is a \code{list} with two objects,
#'  \code{SmallUpsilon} and \code{BigTheta}. \code{SmallUpsilon} represents the degrees of freedom parameter for the
#'  inverse-Wishart hyperprior and must be a real number scalar, while \code{BigTheta} represents
#'  the scale matrix and must be a \code{O x O} dimensional positive definite matrix.
#'  
#'  \code{Rho} is a \code{list} with two objects, \code{ARho} and \code{BRho}. \code{ARho}
#'  represents the lower bound for the uniform hyperprior, while \code{BRho} represents
#'  the upper bound. The bounds must be specified carefully. This is only specified for continuous spatial processes.
#'  
#'  \code{Upsilon} is a \code{list} with two objects,
#'  \code{Zeta} and \code{Omega}. \code{Zeta} represents the degrees of freedom parameter for the
#'  inverse-Wishart hyperprior and must be a real number scalar, while \code{Omega} represents
#'  the scale matrix and must be a \code{K x K} dimensional positive definite matrix.
#'
#'  \code{Psi} is a \code{list} with two objects, dependent on if the temporal kernel is \code{exponential} or \code{ar1}.
#'  For \code{exponential}, the two objects are \code{APsi} and \code{BPsi}. \code{APsi}
#'  represents the lower bound for the uniform hyperprior, while \code{BPsi} represents
#'  the upper bound. The bounds must be specified carefully. For \code{ar1}, the two objets are \code{Beta} and \code{Gamma}, which are the 
#'  two shape parameters of a Beta distribution shifted to have domain in (-1, 1). 
#'  
#' @param Tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the MCMC Metropolis steps. If \code{NULL} is not chosen then all
#'  of the tuning values must be specified.
#'
#'  When \code{NULL} is chosen then default tuning values are automatically generated to
#'  \code{1}. Otherwise a \code{list} must be provided with names \code{Psi}, 
#'  or \code{Rho}. Each of these entries must be scalars containing tuning variances for their corresponding Metropolis updates.
#'  \code{Rho} is only specified for continuous spatial processes.
#'
#' @param MCMC Either \code{NULL} or a \code{list} containing input values to be used
#'  for implementing the MCMC sampler. If \code{NULL} is not chosen then all
#'  of the MCMC input values must be specified.
#'
#'  \code{NBurn}: The number of sampler scans included in the burn-in phase. (default =
#'  \code{10,000})
#'
#'  \code{NSims}: The number of post-burn-in scans for which to perform the
#'   sampler. (default = \code{10,000})
#'
#'  \code{NThin}: Value such that during the post-burn-in phase, only every
#'  \code{NThin}-th scan is recorded for use in posterior inference (For return values
#'  we define, NKeep = NSims / NThin (default = \code{1}).
#'
#'  \code{NPilot}: The number of times during the burn-in phase that pilot adaptation
#'  is performed (default = \code{20})
#'
#' @param Family Character string indicating the distribution of the observed data. Options
#'  include: \code{"normal"}, \code{"probit"}, \code{"tobit"}, and \code{"binomial"}. \code{Family} must have either \code{O} or
#'  \code{1} dimensions (the one populates the rest). Any combination of likelihoods can be used.
#'
#' @param TemporalStructure Character string indicating the temporal kernel. Options include:
#'  \code{"exponential"} and \code{"ar1"}.
#'
#' @param SpatialStructure Character string indicating the type of spatial process. Options include:
#'  \code{"continuous"} (i.e., Gaussian process with exponential kernel) and \code{"discrete"} (i.e., proper CAR).
#'
#' @param Seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'
#' @details Details of the underlying statistical model proposed by
#'  Berchuck et al. 2019. are forthcoming.
#'
#' @return \code{bfa_sp} returns a list containing the following objects
#'
#'   \describe{
#'
#'   \item{\code{lambda}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for factor loadings matrix \code{lambda}.
#'   The labels for each column are Lambda_O_M_K.}
#'
#'   \item{\code{eta}}{\code{NKeep x (Nu x K)} \code{matrix} of posterior samples for the latent factors \code{eta}.
#'   The labels for each column are Eta_Nu_K.}
#'
#'   \item{\code{sigma2}}{\code{NKeep x (M * (O - C))} \code{matrix} of posterior samples for the variances \code{sigma2}.
#'   The labels for each column are Sigma2_O_M.}
#'
#'   \item{\code{kappa}}{\code{NKeep x ((O * (O + 1)) / 2)} \code{matrix} of posterior samples for \code{kappa}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Kappa3_2} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{delta}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{delta}.}
#'
#'   \item{\code{tau}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{tau}.}
#'
#'   \item{\code{upsilon}}{\code{NKeep x ((K * (K + 1)) / 2)} \code{matrix} of posterior samples for \code{Upsilon}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Upsilon3_2} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{psi}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{psi}.}
#'
#'   \item{\code{xi}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for factor loadings cluster labels \code{xi}.
#'   The labels for each column are Xi_O_M_K.}
#'   
#'   \item{\code{rho}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{rho}.}
#'
#'   \item{\code{metropolis}}{\code{2 (or 1) x 3} \code{matrix} of metropolis
#'   acceptance rates, updated tuners, and original tuners that result from the pilot
#'   adaptation.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{bfa_sp} functions
#'   and should be ignored by the user.}
#'
#'   \item{\code{dataug}}{A \code{list} of data augmentation objects that are used in future
#'   \code{bfa_sp} functions and should be ignored by the user.}
#'
#'   }
#'
# @author Samuel I. Berchuck
#' @references Reference for Berchuck et al. 2019 is forthcoming.
#' @export
bfa_sp <- function(Y, Dist, Time, K, L = Inf, Trials = NULL,
                   Starting = NULL, Hypers = NULL, Tuning = NULL, MCMC = NULL, 
                   Family = "normal", TemporalStructure = "exponential", SpatialStructure = "discrete", Seed = 54) {
  
  ###Function Inputs
  # Y = Data
  # Dist = W
  # Time = Time
  # Trials = Trials
  # Starting = Starting
  # Hypers = Hypers
  # Tuning = Tuning
  # MCMC = MCMC
  # Family = "binomial"
  # TemporalStructure = "exponential"
  # SpatialStructure = "discrete"
  # Seed = 54
  # K = K
  # L = Inf
  
  ###Check for missing objects
  if (missing(Y)) stop("Y: missing")
  if (missing(Dist)) stop("Dist: missing")
  if (missing(Time)) stop("Time: missing")
  if (missing(K)) stop("K: missing")

  ###Check model inputs
  CheckInputs(Y, Dist, Time, K, L, Trials, Starting, Hypers, Tuning, MCMC, Family, TemporalStructure, SpatialStructure)

  ####Set seed for reproducibility
  set.seed(Seed)

  ###Check to see if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObj(Y, Dist, Time, Trials, K, L, Family, TemporalStructure, SpatialStructure)
  HyPara <- CreateHyPara(Hypers, DatObj) 
  MetrObj <- CreateMetrObj(Tuning, DatObj)
  Para <- CreatePara(Starting, DatObj, HyPara)
  DatAug <- CreateDatAug(DatObj)
  McmcObj <- CreateMcmc(MCMC, DatObj)
  RawSamples <- CreateStorage(DatObj, McmcObj)

  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run MCMC sampler in Rcpp
  RegObj <- bfa_sp_Rcpp(DatObj, HyPara, MetrObj, Para, DatAug, McmcObj, RawSamples, Interactive)

  ###Set regression objects
  RawSamples <- RegObj$rawsamples
  MetropRcpp <- RegObj$metropolis

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Collect output to be returned
  DatObjOut <- OutputDatObj(DatObj)
  DatAugOut <- OutputDatAug(DatAug)
  Metropolis <- SummarizeMetropolis(DatObj, MetrObj, MetropRcpp, McmcObj)
  Samples <- FormatSamples(DatObj, RawSamples)

  ###Return spBFA object
  spBFA <- list(lambda = Samples$Lambda,
                eta = Samples$Eta,
                sigma2 = Samples$Sigma2,
                kappa = Samples$Kappa,
                delta = Samples$Delta,
                tau = Samples$Tau,
                upsilon = Samples$Upsilon,
                psi = Samples$Psi,
                xi = Samples$Xi,
                rho = Samples$Rho,
                metropolis = Metropolis,
                datobj = DatObjOut,
                dataug = DatAugOut,
                runtime = paste0("Model runtime: ",round(RunTime, 2), " ", attr(RunTime, "units")))
  spBFA <- structure(spBFA, class = "spBFA")
  return(spBFA)

###End sampler
}
