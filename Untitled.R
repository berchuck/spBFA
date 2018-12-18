###Load my library
library(spBFA)
library(womblR)

###Read in data
analdata <- read.csv("/Volumes/Macintosh HD/Users/sam/Box Sync/Postdoc/Projects/SFA/DataApplication/Data/analdata.csv")[, -1]
# analdata <- read.csv("/home/sib2/Projects/SFA/DataApplication/Data/analdata.csv")[, -1]
ID <- unique(analdata$eyeid)
patdata <- analdata[analdata$eyeid == ID[9], ]

###Format data
VF <- t(patdata[, 26:77])
VF[VF < 0] <- 0
RNFL <- t(patdata[, 78:129])
M <- dim(VF)[1]
Nu <- dim(VF)[2]
O <- 2
Y <- array(dim = c(M, O, Nu))
Y[ , 1, ] <- VF
Y[ , 2, ] <- RNFL
Time <- patdata$sapfollowup

###Spatial distance
blind_spot <- c(26, 35)
W <- HFAII_Queen[-blind_spot, -blind_spot]

###Bounds for psi
TimeDist <- as.matrix(dist(Time))
BPsi <- log(0.025) / -min(TimeDist[TimeDist > 0])
APsi <- log(0.975) / -max(TimeDist)

###Initial values
K <- 25
Starting <- list(Sigma2 = 1,
                 Kappa = diag(O),
                 Delta = rep(1, K),
                 Psi = (APsi + BPsi) / 8,
                 Upsilon = diag(K))

###Hyperparameters
Hypers <- list(Sigma2 = list(A = 0.001, B = 0.001),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Delta = list(A1 = 1, A2 = 2),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))

###Metropolis tuners
Tuning <- list(Psi = 1)

###MCMC inputs
MCMC <- list(NBurn = 50000, NSims = 250000, NThin = 25, NPilot = 100)

###Fit Bayesian factor analysis model
reg.bfa_sp <- bfa_sp(Y = Y, Dist = W, Time = Time, K = K, Starting = Starting, Hypers = Hypers, Tuning = Tuning, MCMC = MCMC, ScaleY = 10, Family = c("tobit", "normal"))
