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
plot(1:10, 1:10, xlab = "", ylab = "", type = "n", xaxt = "n", yaxt = "n", bty = "n", asp = 1)
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





