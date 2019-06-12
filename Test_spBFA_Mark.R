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
M <- length(unique(analdata$ubigeo))
Types <- c("nfalciparum", "nvivax")
O <- length(Types)
analdata2015 <- analdata[(analdata$year == 2015) & (analdata$epiweek %in% 1:10), ]
analdata2015$season <- 1 * (as.numeric(format(as.Date(analdata2015$date, format = "%m/%d/%Y"), "%m")) %in% c(2, 3, 4, 5, 6, 7))


Nu <- dim(analdata2015)[1] / M
N <- M * O * Nu
Space <- unique(analdata2015$ubigeo)
Time <- unique(analdata2015$epiweek)
dat <- matrix(nrow = N, ncol = 8)
for (t in 1:Nu) {
  for (o in 1:O) {
    for (i in 1:M) {
      Index <- i + (o - 1) * M + M * O * (t - 1)
      dat[Index, 1] <- analdata2015[analdata2015$ubigeo == Space[i] & analdata2015$epiweek == Time[t], Types[o]]
      dat[Index, 2] <- analdata2015[analdata2015$ubigeo == Space[i] & analdata2015$epiweek == Time[t], "population"]
      dat[Index, 3] <- analdata2015[analdata2015$ubigeo == Space[i] & analdata2015$epiweek == Time[t], "rainMM"]
      dat[Index, 4] <- analdata2015[analdata2015$ubigeo == Space[i] & analdata2015$epiweek == Time[t], "TmeanC"]
      dat[Index, 5] <- analdata2015[analdata2015$ubigeo == Space[i] & analdata2015$epiweek == Time[t], "season"]
      dat[Index, 6] <- t
      dat[Index, 7] <- o
      dat[Index, 8] <- i
    }
  }
}
dat <- data.frame(dat)
colnames(dat) <- c("malaria", "population", "rain", "temp", "season", "time", "type", "location")

Time <- unique(dat$time / Nu) 
TimeDist <- as.matrix(dist(Time))
BPsi <- log(0.025) / -min(TimeDist[TimeDist > 0])
APsi <- log(0.975) / -max(TimeDist)
K <- round(log(M * O))
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
reg.bfa_spTEST <- bfa_sp(malaria ~ rain + temp + season, data = dat, family = "binomial", trials = "population", dist = W, time = Time, K = K, starting = Starting, hypers = Hypers, tuning = Tuning, mcmc = MCMC)

reg.bfa_sp <- bfa_sp(malaria ~ rain + temp + season, data = dat, family = "binomial", trials = "population", dist = W, time = Time, K = K, starting = Starting, hypers = Hypers, tuning = Tuning, mcmc = MCMC)

pred <- predict(reg.bfa_sp, NewTimes = 0.1923, type = "temporal", X = X[reg.bfa_sp$datobj$Indeces == (Nu - 1), ], NewTrials = array(reg.bfa_sp$datobj$Trials, dim = c(M, O, 1)))

###Check posterior covariances
X <- model.matrix(malaria ~ rain + temp + season, data = dat)
Mean <- matrix(0, nrow = M * O, ncol = Nu)
for (s in 1:dim(reg.bfa_sp$rho)[1]) {
  Lambda <- matrix(reg.bfa_sp$lambda[s, ], nrow = M * O, ncol = K, byrow = TRUE)
  Eta <- matrix(reg.bfa_sp$eta[s, ], ncol = 1)
  Beta <- matrix(reg.bfa_sp$beta[s, ], ncol = 1)
  Mean <- Mean + matrix(kronecker(diag(Nu), Lambda) %*% Eta + X %*% Beta, nrow = M * O, ncol = Nu)
}
Mean <- Mean / dim(reg.bfa_sp$rho)[1]


time <- 10
Yt <- matrix(reg.bfa_sp$datobj$YObserved, nrow = M * O, ncol = Nu)[52:102, time]
Trialst <- matrix(reg.bfa_sp$datobj$Trials, nrow = M * O, ncol = Nu)[52:102, time]
Predt <- Trialst * exp(Mean[52:102, time]) / (1 + exp(Mean[52:102, time]))

Predt2 <- apply(pred$Y$Y11, 2, mean)[52:102]

Yt <- matrix(reg.bfa_sp$datobj$YObserved, nrow = M * O, ncol = Nu)[1:51, time]
Trialst <- matrix(reg.bfa_sp$datobj$Trials, nrow = M * O, ncol = Nu)[1:51, time]
Predt <- Trialst * exp(Mean[1:51, time]) / (1 + exp(Mean[1:51, time]))

Predt2 <- apply(pred$Y$Y11, 2, mean)[1:51]



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

