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
O <- 2
M <- 52 
Nu <- dim(VFSeries)[1] / M
YWide <- matrix(Y, nrow = M, ncol = Nu)
YWide2 <- YWide + pmax(matrix(sample(-10:10, M * Nu, replace = TRUE), nrow = M, ncol = Nu), 0)
Data <- array(dim = c(M, O, Nu))
Data[ , 1, ] <- YWide
Data[ , 2, ] <- YWide2
Time <- unique(VFSeries$Time) / 365 #time since first visit
K <- 5
L <- Inf
W <- HFAII_Queen[-c(26, 35), -c(26, 35)] #Visual field adjacency matrix
Trials <- array(dim = c(M, 1, Nu))
Trials[ , , ] <- 50

###Load simulated dataset
# load("/Volumes/Macintosh HD/Users/sam/Box Sync/Postdoc/Projects/SFA/Simulations/Simulation1/Data/SimData.RData")
# Time <- SimData[[3]]
# Y <- matrix(SimData[[1]][, 1], ncol = 1)
# W <- SimData[[2]]
# K <- 10

###Center data
# YWide <- matrix(Y, nrow = 100, ncol = 21)
# YCenter <- matrix(scale(YWide, scale = FALSE), ncol = 1)

###Bounds for temporal tuning parameter
TimeDist <- abs(outer(Time, Time, "-"))
minDiff <- min(TimeDist[TimeDist > 0])
maxDiff <- max(TimeDist[TimeDist > 0])
APsi <- -log(0.95) / maxDiff #longest diff goes up to 95%
BPsi <- -log(0.01) / minDiff #shortest diff goes down to 1%
# APsi <- -1
# BPsi <- 1
# Beta <- Gamma <- 1

###Bounds for spatial tuning parameter
# SpDist <- as.matrix(dist(rnorm(100, 10, 1)))
# minDiff <- min(SpDist[SpDist > 0])
# maxDiff <- max(SpDist[SpDist > 0])
# ARho <- -log(0.95) / maxDiff #longest diff goes up to 95%
# BRho <- -log(0.05) / minDiff #shortest diff goes down to 1%

###Initial values
Starting <- list(
                  # Sigma2 = 1,
                 Kappa = diag(O),
                 Rho = 0.99,
                 Delta = 2 * (1:K),
                 Psi = 1,
                 Upsilon = diag(K))

###Hyperparameters
Hypers <- list(
              # Sigma2 = list(A = 0.001, B = 0.001),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Delta = list(A1 = 1, A2 = 20),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))

###Metropolis tuners
Tuning <- list(Psi = 1)

###MCMC inputs
# MCMC <- list(NBurn = 10000, NSims = 250000, NThin = 25, NPilot = 10)
MCMC <- list(NBurn = 100, NSims = 250, NThin = 1, NPilot = 2)

###Fit sampler
reg.bfa_sp <- bfa_sp(Y = Data, Dist = W, Time = Time, K = K, L = Inf, Trials = Trials,
                     Starting = Starting, Hypers = Hypers, Tuning = Tuning, MCMC = MCMC,
                     TemporalStructure = "exponential", SpatialStructure = "discrete",
                     Family = c("binomial", "normal"))



save(reg.bfa_sp, file = "/Users/sam/Desktop/out.RData")

# # load("/Users/sam/Desktop/Sim102.RData")
# load("/Users/sam/Desktop/Sim133.RData")
# 
# ###Check a prediction
# # s <- 5000
# # K <- 4
# # Nu <- 21
# # M <- 100
# # NSims <- 10000
# s <- 500
# K <- 10
# Nu <- 21
# M <- 100
# NSims <- 10000
# Lambda <- matrix(reg.bfa_sp$lambda[s, ], nrow = M, ncol = K, byrow = TRUE)
# Eta <- matrix(as.numeric(matrix(reg.bfa_sp$eta[s, ], nrow = K, ncol = Nu, byrow = TRUE)), ncol = 1)
# Mean <- matrix(kronecker(diag(Nu), Lambda) %*% Eta, nrow = M, ncol = Nu)
# par(mfcol = c(1, 2))
# Yt <- matrix(reg.bfa_sp$datobj$YObserved, nrow = M, ncol = Nu)[, 1]
# plot.sim(Mean[, 1], zlim = c(min(Yt), max(Yt)))
# # plot.sim(matrix(YCenter, nrow = M, ncol = Nu)[, 1])
# plot.sim(Yt)
# 
# 
# 
# ###Check posterior covariances
# BigPsi <- matrix(0, nrow = M, ncol = M)
# Mean <- matrix(0, nrow = M, ncol = Nu)
# for (s in 1:NSims) {
#   Lambda <- matrix(reg.bfa_sp$lambda[s, ], nrow = M, ncol = K, byrow = TRUE)
#   Eta <- matrix(as.numeric(matrix(reg.bfa_sp$eta[s, ], nrow = K, ncol = Nu, byrow = TRUE)), ncol = 1)
#   Mean <- Mean + matrix(kronecker(diag(Nu), Lambda) %*% Eta, nrow = M, ncol = Nu)
#   Upsilon <- matrix(0, nrow = K, ncol = K)
#   Upsilon[upper.tri(Upsilon, diag = TRUE)] <- reg.bfa_sp$upsilon[s, ]
#   Upsilon[lower.tri(Upsilon)] <- t(Upsilon)[lower.tri(Upsilon)]
#   Sigma <- diag(reg.bfa_sp$sigma2[s, ])
#   BigPsi <- BigPsi + cov2cor(Lambda %*% Upsilon %*% t(Lambda) + Sigma)
#   # image(cov2cor(BigPsi)[index, index])
#   # exp(-reg.bfa_sp$psi[s] * TimeDist)[1, 1]
# }
# par(mfcol = c(1, 1))
# image(BigPsi[index, rev(index)] / NSims)
# image(BigPsi / NSims)
# 
# 
# par(mfcol = c(1, 2))
# plot.sim((Mean / NSims)[, 1], zlim = c(min(Yt), max(Yt)))
# plot.sim(matrix(YCenter, nrow = 100, ncol = 21)[, 1])
# 
# 
# 
# ###Check heatmap probabilities
# out <- array(0, dim = c(M, M, K))
# for (s in 1:NSims) {
#   Xis <- matrix(reg.bfa_sp$xi[s, ], nrow = M, ncol = K, byrow = TRUE)
#   # Lambda <- matrix(reg.bfa_sp$lambda[s, ], nrow = M, ncol = K, byrow = TRUE)
#   # Eta <- matrix(as.numeric(matrix(reg.bfa_sp$eta[s, ], nrow = K, ncol = Nu, byrow = TRUE)), ncol = 1)
#   # matrix(Eta, nrow = K, ncol = Nu)
#   for (j in 1:K) {
#     mat <- out[ , , j]
#     mat <- mat + 1 * (as.matrix(dist(Xis[, j])) == 0)
#     out[ , , j] <- mat
#   }
# }
# out <- out / NSims
# 
# index <- Indices[order(Indices[, 3]), 4]
# 
# par(mfcol = c(2, 5))
# for (i in 1:10) image(out[index, rev(index), i])
# 
# image(out[index, rev(index), 1])
# image(out[index, rev(index), 2])
# image(out[index, rev(index), 3])
# image(out[index, rev(index), 4])
# 
# 
# image(out[, , 1])
# 
# ###Posterior checks
# library(coda)
# par(mfcol = c(1, 1))
# traceplot(as.mcmc(reg.bfa_sp$psi))
# traceplot(as.mcmc(reg.bfa_sp$lambda[, 5]))
# traceplot(as.mcmc(reg.bfa_sp$eta[, 1]))
# 
# traceplot(as.mcmc(reg.bfa_sp$xi[, 2]))
# 
# 
# traceplot(as.mcmc(reg.bfa_sp$sigma2[, 2]))
# traceplot(as.mcmc(reg.bfa_sp$kappa2))
# par(mfcol = c(2, 2))
# traceplot(as.mcmc(1 / reg.bfa_sp$tau))
# 
# 
# install.packages("RcppArmadillo", lib = "/home/sib2/R")
# install.packages("Cairo", lib = "/home/sib2/R")
# 
# install.packages("spBFA_1.0.tar.gz", lib = "/home/sib2/R", type = "source", repo = NULL)
# .libPaths( c( .libPaths(), "/home/sib2/R") )
# install.packages("/home/sib2/Packages/spBFA_1.0.tar.gz", lib = "/home/sib2/R", type = "source", repo = NULL)
# install.packages("womblR_1.0.3.tar.gz", lib = "/home/sib2/R", type = "source", repo = NULL)
