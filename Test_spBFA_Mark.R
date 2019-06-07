rm(list = ls()) #Start with a clean working directory

###Load Library (set the path that you save the package. You end with the package in the file path)
path.package <- "/Users/sam/Documents/Postdoc/Software/spBFA/"
devtools::load_all(path.package, export_all = TRUE)
devtools::document(path.package)
library(spBFA)
library(womblR)

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
analdata <- read.csv("/Users/sam/Box Sync/Postdoc/Projects/SFA/Mark/Data/malaria_data_for_Sam.csv")
W <- read.csv("/Users/sam/Box Sync/Postdoc/Projects/SFA/Mark/Data/adj_matrix.csv", header = FALSE)[, -1]
W <- as.matrix(W)
M <- 51
O <- 2
# Nu <- 52 * 1
Nu <- 10
Data <- Trials <- array(dim = c(M, O, Nu))
Data1 <- matrix(analdata$nfalciparum[analdata$year %in% 2015 & analdata$epiweek %in% 1:10], nrow = M, ncol = Nu)
Data2 <- matrix(analdata$nvivax[analdata$year %in% 2015 & analdata$epiweek %in% 1:10], nrow = M, ncol = Nu)
Trial <- matrix(analdata$population[analdata$year %in% 2015 & analdata$epiweek %in% 1:10], nrow = M, ncol = Nu)
Data[ , 1, ] <- Data1
Data[ , 2, ] <- Data2
Trials[ , 1, ] <- Trials[ , 2, ] <- Trial
Time <- seq(0, 1, length.out = Nu)
TimeDist <- as.matrix(dist(Time))
BPsi <- log(0.025) / -min(TimeDist[TimeDist > 0])
APsi <- log(0.975) / -max(TimeDist)
K <- 10
Starting <- list(Kappa = diag(O),
                 Delta = 2 * (1:K),
                 Psi = (APsi + BPsi) / 2,
                 Upsilon = diag(K))
Hypers <- list(Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Delta = list(A1 = 1, A2 = 20),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
Tuning <- list(Psi = 1)
MCMC <- list(NBurn = 100, NSims = 100, NThin = 1, NPilot = 2)
reg.bfa_sp <- bfa_sp(Y = Data, Trials = Trials, Dist = W, Time = Time, K = K, Starting = Starting, Hypers = Hypers, Tuning = Tuning, MCMC = MCMC, Family = "binomial")

pred <- predict(reg.bfa_sp, NewTimes = 0.999, type = "temporal")

###Check posterior covariances
Mean <- matrix(0, nrow = M * O, ncol = Nu)
for (s in 1:dim(reg.bfa_sp$rho)[1]) {
  Lambda <- matrix(reg.bfa_sp$lambda[s, ], nrow = M * O, ncol = K, byrow = TRUE)
  Eta <- matrix(reg.bfa_sp$eta[s, ], ncol = 1)
  Mean <- Mean + matrix(kronecker(diag(Nu), Lambda) %*% Eta, nrow = M * O, ncol = Nu)
}
Mean <- Mean / dim(reg.bfa_sp$rho)[1]




time <- 10
Yt <- matrix(reg.bfa_sp$datobj$YObserved, nrow = M * O, ncol = Nu)[1:51 , time]
Trialst <- matrix(reg.bfa_sp$datobj$Trials, nrow = M * O, ncol = Nu)[1:51, time]
Predt <- Trialst * exp(Mean[1:51, time]) / (1 + exp(Mean[1:51, time]))

Predt2 <- apply(pred$Y$Y11, 2, mean)[1:51]

quantities_of_interest <- rnorm(51)
### libraries
library(tidyverse)
library(sf)

### read in polygon
polygon <- st_read(dsn = "/Users/sam/Box Sync/Postdoc/Projects/SFA/Mark/Data/distritoloreto/distritoloreto.shp")

### plot polygon
gg <- ggplot(data = polygon) +
  geom_sf(aes(fill = Yt)) +
  theme_bw() +
  labs(fill = paste0("True: ", time))
gg

gg <- ggplot(data = polygon) +
  geom_sf(aes(fill = Predt)) +
  theme_bw() +
  labs(fill = paste0("Pred: ", time))
gg

gg <- ggplot(data = polygon) +
  geom_sf(aes(fill = Predt2)) +
  theme_bw() +
  labs(fill = paste0("Pred Krig: ", time))
gg

