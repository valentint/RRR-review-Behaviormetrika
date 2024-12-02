library(microbenchmark)
library(rrpack)
library(secure)

here::i_am("timing-microbenchmark.R")
library(here)

source("ccsrrr.R")

n <- 100; p <- 30; p0 <- 10; q <- 10; q0 <- 10; nrank = 3; rho_X = 0; rho_E = 0.5; # 1d

res1 <- rrr.sim2(n, p, p0, q, q0, nrank, s2n = 1, sigma = NULL, rho_X, rho_E)
         
X <- res1$X
Y <- res1$Y

X <- scale(X, center = TRUE, scale = FALSE)/sqrt(n-1) # whitened data
Y <- scale(Y, center = TRUE, scale = FALSE)/sqrt(n-1) # whitened data
XY <- t(X) %*% Y

nbr <- 2

fit1 <- ccsrrr(Y, X, XY, nrank, nbr) 
fit2 <- secure.path(Y, X, nrank, nlambda = 300,  control=secure.control(elnetAlpha=1))

res <- microbenchmark(csrrr=ccsrrr(Y, X, nrank, nbr), 
                      secure=secure.path(Y, X, nrank, nlambda = 300,  control=secure.control(elnetAlpha=1)), times=100L) 
res
plot(res)

