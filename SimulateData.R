###Start with a clean space
rm(list = ls())

###Set seed
set.seed(54)

###Set simulation parameters
Grid <- 10
M <- Grid^2
Nu <- 21
GroupMeans <- c(-5, -2, 2, 5)
K <- length(GroupMeans)
Time <- seq(0, 1, length.out = Nu)

###Create indices
Indices <- cbind(rep(1:10, 10), rep(1:10, each = 10), 0)
Indices <- Indices[order(Indices[, 1]), ]
Indices[Indices[, 2] <= 5, 3] <- c(rep(3, M / 4), rep(4, M / 4))
Indices[Indices[, 2] > 5, 3] <- c(rep(1, M / 4), rep(2, M / 4))
Indices <- Indices[order(Indices[, 2], decreasing = TRUE), ]

###Visualize grid
plot(1:11, 1:11, xlab = "", ylab = "", type = "n", xaxt = "n", yaxt = "n", bty = "n", asp = 1)
counter <- 0
for (i in 1:M) {
  counter <- counter + 1
  points(Indices[counter, 1], Indices[counter, 2], pch = Indices[counter, 3])
}

###Create adjacency matrix
library(hierarchicalDS)
W <- square_adj(Grid)

###Set true lambda parameters
Lambda <- matrix(0, nrow = M, ncol = K)
for (i in 1:M) Lambda[i, Indices[i, 3]] <- 1

###Temporal components / factors
Psi <- 0.9
Upsilon <- diag(K)
BigPhi <- matrix(nrow = K, ncol = Nu)  
BigPhi[, 1] <- GroupMeans
for (t in 2:Nu) BigPhi[, t] <- mvtnorm::rmvnorm(1, Psi * BigPhi[, t - 1], Upsilon)  
Eta <- matrix(as.numeric(BigPhi), ncol = 1)

###Variances
Sigma2 <- as.numeric(exp(mvtnorm::rmvnorm(1, rep(0, M), diag(M))), ncol = 1)

###Simulated moments
JointMean <- kronecker(diag(Nu), Lambda) %*% Eta
JointVar <- rep(Sigma2, Nu)

###Simulations
NSims <- 200
Sims <- matrix(nrow = Nu * M, ncol = NSims)
for (i in 1:NSims) Sims[, i] <- matrix(rnorm(Nu * M, JointMean, JointVar), ncol = 1)

###Save data
SimData <- list(Sims, W, Time)


###Visualize simulated data
plot.sim <- function(Y) {
  zlim <- c(min(Y), max(Y))
  bins <- 100
  labs <- levels(cut(zlim, bins))
  labs <- cbind(lower = as.numeric(sub("\\((.+),.*","\\1", labs)), upper = as.numeric(sub("[^,]*,([^]]*)\\]","\\1", labs)))
  legvals <- as.numeric(c(labs[1, 1], labs[ , 2]))
  legvals[1] <- -Inf
  legvals[length(legvals)] <- Inf
  colbr <- colorRampPalette(c("yellow", "orange", "red"))
  colpal <- colbr(bins)
  cuts <- cut(Y[!is.na(Y)], breaks = legvals)
  cuts <- colpal[as.numeric(cuts)]
  plot(1:11, 1:11, xlab = "", ylab = "", type = "n", xaxt = "n", yaxt = "n", bty = "n", asp = 1)
  counter <- 0
  for (i in 1:M) {
    counter <- counter + 1
    # points(Indices[counter, 1], Indices[counter, 2], pch = Indices[counter, 3])
    x <- Indices[counter, 1]
    y <- Indices[counter, 2]
    rect(x, y, x + 1, y + 1, col = cuts[i], border = cuts[i])
  }
  NColors <- length(colpal)
  Vertical <- seq(4, 8, length.out = NColors)
  for (i in 1:NColors) segments(12, Vertical[i], 12.75, Vertical[i], col = colpal[i], lwd = 1.5)
  minx <- zlim[1]
  maxx <- zlim[2]
  LegendPV <- seq(minx, maxx, length.out = 5)
  segments(12.75, 4, 12.75, 8, lwd = 1.5)
  segments(12, 4, 12, 8, lwd = 1.5)
  segments(12, 8, 12.75, 8, lwd = 1.5)
  segments(12, 4, 12.75, 4, lwd = 1.5)
  format0 <- function(x, legend.round) format(round(x,legend.round),nsmall=legend.round)
  for (i in 1 : length(LegendPV)) {
    text(13.75, (4:8)[i], format0(LegendPV[i], 0))
    segments(12.75, (4:8)[i], 13, (4:8)[i], lwd = 1.5)
  }
  text(12.5, 8.5, "Y")
}
Y <- matrix(Sims[, 1], nrow = M, ncol = Nu)
plot.sim(Y[, 1])


