rm(list = ls()) #Start with a clean working directory

###Load Library (set the path that you save the package. You end with the package in the file path)
path.package <- "/Users/sam/Documents/Postdoc/Software/spBFA/"
devtools::load_all(path.package, export_all = TRUE)
devtools::document(path.package)
library(spBFA)

# devtools::build_win(pkg = path.package)
# devtools::use_vignette("spCP-example", pkg = path.package)
# devtools::use_build_ignore("spCP-example.pdf", pkg = path.package)
# devtools::use_build_ignore(files = "cran-comments.md", pkg = path.package)
# devtools::use_build_ignore("NEWS.md", pkg = path.package)
# devtools::use_build_ignore("README.md", pkg = path.package)
# devtools::build_win(pkg = path.package)
# devtools::spell_check()
# devtools::release(path.package) # for submitting package for the first time
# devtools::submit_cran(path.package) # submit to CRAN without going through checks

###Format data for MCMC sampler
library(womblR)
VFSeries <- VFSeries[order(VFSeries$Visit), ] #sort by visit
VFSeries <- VFSeries[!VFSeries$Location %in% c(26, 35), ] #remove blind spot locations
Y <- VFSeries$DLS #assign observed outcome data
Time <- unique(VFSeries$Time) / 365 #time since first visit
Nu <- length(Time)
K <- 6
L <- 50

###Create original adjacency matrix
W <- HFAII_Queen[-c(26, 35), -c(26, 35)] #Visual field adjacency matrix

###Bounds for temporal tuning parameter
TimeDist <- abs(outer(Time, Time, "-"))
minDiff <- min(TimeDist[TimeDist > 0])
maxDiff <- max(TimeDist[TimeDist > 0])
APsi <- -log(0.95) / maxDiff #longest diff goes up to 95%
BPsi <- -log(0.01) / minDiff #shortest diff goes down to 1%

###Initial values
Starting <- list(Sigma2 = 1,
                 Kappa2 = 1,
                 Delta = rep(1, K),
                 Psi = 1,
                 Upsilon = diag(K))

###Hyperparameters
Hypers <- list(Sigma2 = list(A = 0.001, B = 0.001),
               Kappa2 = list(C = 0.001, D = 0.001),
               Delta = list(A1 = 1, A2 = 2),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))

###Metropolis tuners
Tuning <- list(Psi = 1)

###MCMC inputs
MCMC <- list(NBurn = 1000, NSims = 2500, NThin = 1, NPilot = 10)

###Fit sampler
reg.bfa_sp <- bfa_sp(Y = Y, W = W, Time = Time, K = K, L = L,
                     Starting = Starting, Hypers = Hypers, Tuning = Tuning, MCMC = MCMC,
                     ScaleY = 1, Rho = 0.99, Family = "tobit", Seed = 54)

###Posterior checks
library(coda)
traceplot(as.mcmc(reg.bfa_sp$psi))
traceplot(as.mcmc(reg.bfa_sp$lambda[, 100]))
traceplot(as.mcmc(reg.bfa_sp$eta[, 6]))

traceplot(as.mcmc(reg.bfa_sp$sigma2[, 2]))
traceplot(as.mcmc(reg.bfa_sp$kappa2))
par(mfcol = c(2, 3))
traceplot(as.mcmc(1 / reg.bfa_sp$tau))



