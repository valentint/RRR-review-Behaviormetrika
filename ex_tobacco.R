##  Example: Tobacco dataset
##  
##  Anderson, R. L. and Bancroft, T. A. (1952). Statistical Theory in Research. McGraw-Hill, New York, NY. 
##      https://dn790002.ca.archive.org/0/items/StatisticalTheoryInResearch/StatisticalTheoryInResearch.pdf.
##  p.205
##
##  Izenman, A. J. (2008). Modern Multivariate Statistical Techniques. Springer, New York, NY. 6.3.3
##

library(xtable)

here::i_am("ex_tobacco.R")
library(here)

source("ccsrrr.R")


##  Read the data
df <- read.csv("tobacco.csv")
df <- df[ , -1]
head(df)
dim(df)
colnames(df) <- c(paste0("y", 1:3), paste0("x", 1:6)) 
tobacco <- df
save(tobacco, file="tobacco.rda", version=2)

X <- as.matrix(tobacco[, 4:9])
Y <- as.matrix(tobacco[,1:3]) 
n <- nrow(X)

## Scale/whiten the data
X <- scale(X, center = TRUE, scale = FALSE)/sqrt(n-1) # whitened data
Y <- scale(Y, center = TRUE, scale = FALSE)/sqrt(n-1) # whitened data

##
## Prepare Table 1
##
##======================================================
W <- svd(X)
s <- sum(W$d > 1e-6)
U <- W$u[, 1:s, drop=FALSE]
V <- W$v[, 1:s, drop=FALSE] 
D <- W$d[1:s]

XX <- t(X) %*% X
XY <- t(X) %*% Y
XXinv <- solve(XX)
XXinvsq <- V %*% diag(1/D) %*% t(V)

## Set the required rank
r <- 1

##=====================================================
##  Eq. (6)
## (X'X)^(-1/2)(X'Y ) = UDV
WW <- svd(XXinvsq %*% XY)
Ur <- WW$u[, 1:r, drop=FALSE]
Vr <- WW$v[, 1:r, drop=FALSE] 
Dr <- WW$d[1:r]

(C_6 <- XXinvsq %*% Ur %*% Dr %*% t(Vr))
(rf_6 <- norm(Y - X %*% C_6, "F"))

##=====================================================
## Eq. (10)
## (X'Y ) = UDV
WW <- svd(t(X) %*% Y)
Ur <- WW$u[, 1:r, drop=FALSE]
Vr <- WW$v[, 1:r, drop=FALSE] 
Dr <- WW$d[1:r]

(C_10 <- XXinv %*% Ur %*% Dr %*% t(Vr))
(rf_10 <- norm(Y - X %*% C_10, "F"))

## Table 1
(C <- cbind(C_6,C_10))
xtable(C)

##
##  Prepare Table 2
##
##==============================================================
##  Chen, L. and Huang, J. Z. (2012). Sparse reduced-rank regression for simultaneous
##  dimension reduction and variable selection. Journal of the American
##  Statistical Association, 107:1533--1545.

library(rrpack)

(rfit1 <- srrr(Y, X, nrank=1))
summary(rfit1)
C1 <- coef(rfit1) 
dimnames(C1) <- list(1:6, 1:3) 
C1
(rf1 <- norm(Y - X %*% C1, "F"))

##===============================================================
## Chen, K., Chan, K.-S., and Stenseth, N. C. (2012). Reduced rank stochastic
##  regression with a sparse singular value decomposition. Journal of the Royal
##  Statistical Society: Series B, 74:203--221.
##
##  rssvd{rrpack}
##

library(rrpack)

(rfit2 <- rssvd(Y, X, nrank=1))
summary(rfit2)
C2 <- coef(rfit2)
dimnames(C2) <- list(1:6, 1:3) 
C2
(rf2 <- norm(Y - X %*% C2, "F"))

##==============================================================
##  Mishra, A., Dey, D. K., and Chen, K. (2017). Sequential co-sparse factor
##  regression. Journal of Computational and Graphical Statistics, 26:814{825.
##  https://doi.org/10.1080/10618600.2017.1340891.

library(secure)

control <- secure.control(spU=100/p, spV=1)
rfit3 <- secure.path(Y, X, nrank=1, control=control)
C3 <- rfit3$C.est
dimnames(C3) <- list(1:6, 1:3) 
C3
(rf3 <- norm(Y - X %*% C3, "F"))


##  Table 2
(C <- cbind(C1, C2, C3))
xtable(C)

##=============================================================
##
##  Prepare Table 3
##
## Simulate the solution of Chen and Huang, 2012
rfit4 <- ccsrrr(Y, X, nrank=1, nonzR=5, nonzC=0)
(C4 <- rfit4$C)
rfit4$f

## Simulate the solution of Chen et al., 2012
rfit5 <- ccsrrr(Y, X, nrank=1, nonzR=5, nonzC=2)
(C5 <- rfit5$C)
rfit5$f

## The default solution with rank 1
rfit6 <- ccsrrr(Y, X, nrank=1)
(C6 <- rfit6$C)
rfit6$f

##  Table 3
(C <- cbind(C4, C5, C6))
xtable(C)


