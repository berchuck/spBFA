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
load("/Users/sam/Box Sync/Postdoc/Projects/SFA/Glaucoma/Data/analdata.RData")
TD <- as.numeric(analdata[44, 7:58])
M <- length(TD)
O <- 1
Nu <- 1
Y <- array(dim = c(M, O, Nu))
Y[ , 1, 1] <- TD
Y <- Y / 10
dat <- data.frame(TD = Y)

Time <- 0

blind_spot <- c(26, 35)
W <- HFAII_Queen[-blind_spot, -blind_spot]
TimeDist <- as.matrix(dist(Time))
BPsi <- 1
APsi <- 0
K <- 10
Starting <- list(Sigma2 = 1,
                 Kappa = diag(O),
                 Delta = 2 * (1:K),
                 Psi = (APsi + BPsi) / 2,
                 Upsilon = diag(K))
Hypers <- list(Sigma2 = list(A = 0.001, B = 0.001),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Delta = list(A1 = 1, A2 = 20),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))
Tuning <- list(Psi = 1)
MCMC <- list(NBurn = 1000, NSims = 1000, NThin = 1, NPilot = 10)
reg.bfa_sp <- bfa_sp(TD ~ 0, data = dat, dist = W, time = Time, K = K, starting = Starting, hypers = Hypers, tuning = Tuning, mcmc = MCMC,
                     gamma.shrinkage = FALSE, include.space = TRUE, clustering = FALSE)

Yt <- matrix(reg.bfa_sp$datobj$YObserved, nrow = M * O, ncol = Nu)


makeSquare <- function(x, bg = 0) {
  square <- matrix(NA, nrow = 9, ncol = 8)
  x <- as.numeric(x)
  x <- c(x[1:25], bg, x[26:33], bg, x[34:52])
  layout.matrix <- matrix(c(bg, bg, bg,  1,  2,  3,  4, bg, bg,
                            bg, bg,  5,  6,  7,  8,  9, 10, bg,
                            bg, 11, 12, 13, 14, 15, 16, 17, 18,
                            19, 20, 21, 22, 23, 24, 25, 26, 27,
                            28, 29, 30, 31, 32, 33, 34, 35, 36,
                            bg, 37, 38, 39, 40, 41, 42, 43, 44,
                            bg, bg, 45, 46, 47, 48, 49, 50, bg,
                            bg, bg, bg, 51, 52, 53, 54, bg, bg), nrow = 8, ncol = 9, byrow = TRUE)
  square[t(layout.matrix != bg)] <- x
  square[is.na(square)] <- bg
  square <- t(square)
  square <- cbind(rep(bg, 12), rep(bg, 12), rbind(rep(bg, 9), rep(bg, 9), square, rep(bg, 9), rep(bg, 9)), rep(bg, 12))
  return(square)
}


library(reshape2)
library(ggplot2)
library(scales)

for (i in 1:dim(analdata)[1]){

TD <- analdata[44, 7:58]
  
pdf(paste0("/Users/sam/Desktop/figs/", i, ".pdf"), height = 4, width = 5)
TD <- makeSquare(TD, bg = -37)
TD[TD == -37] <- NA
data <- melt(t(TD))
p <- ggplot(data, aes(Var1, Var2, fill = value))
p <- p + geom_tile()
p <- p + scale_fill_gradient2(low = "black", high = "red", oob = squish,
                              na.value = "black", limits = c(-40, 10))
p <- p + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank())
print(p)
dev.off()

}


Mean <- apply(diagnostics(reg.bfa_sp, keepPPD = TRUE)$PPD, 2, mean)


pdf(paste0("/Users/sam/Desktop/Mean.pdf"), height = 4, width = 5)
TD <- makeSquare(Mean * 10, bg = -37)
TD[TD == -37] <- NA
data <- melt(t(TD))
p <- ggplot(data, aes(Var1, Var2, fill = value))
p <- p + geom_tile()
p <- p + scale_fill_gradient2(low = "black", high = "red", oob = squish,
                              na.value = "black", limits = c(-40, 10))
p <- p + theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank())
print(p)
dev.off()
