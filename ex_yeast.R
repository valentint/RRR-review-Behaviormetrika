##
##  Yeast cell cycle data
##
##  We use all the n=542 (out of 1790) genes, to examine the association between 
##  the RNA levels along the q=18 time points and the binding information 
##  of p=106 (out of 113) transcription factors (TF). The values n=542 and p=106
##  result after removing the missing data from the original data set of n=1790 and p=113.
##

library(rrpack) 
library(spls)
library(remMap)
library(secure)
library(MRCE)
library(SiER)

here::i_am("yeast.R")
library(here)

source("ccsrrr.R")

data(yeast)
X <- yeast$x
Y <- yeast$y

n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)

##  The 21 known yeast TFs related to cell cycle process.
##      Note that "CBF1" and "GCN4" were not selected by gSCAD.

TF_gSCAD <- c("ACE2", "SWI4", "SWI5", "SWI6", "MBP1", "STB1", "FKH1", 
              "FKH2", "NDD1", "MCM1", "ABF1", "BAS1", "CBF1", "GCN4",
              "GCR1", "GCR2", "LEU3", "MET31", "REB1", "SKN7", "STE12")

tf <- colnames(X)
TF <- substr(tf, 1, nchar(tf)-4)
colnames(X) <- TF
time <- as.numeric(substr(colnames(Y), 6, 9))

##============================================================
## SPLS with eta=0.7 and 8 hidden components
(cv <- cv.spls(X, Y, fold=10, K=c(1:10), eta = seq(0.1,0.9,0.1), plot.it=FALSE)) 
(f <- spls(X, X, K=7, eta=0.7))
Abar <- f$betahat
(TF_spls <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_spls)
(TF_confirmed_spls <- TF_spls[which(TF_spls %in% TF_gSCAD)])
length(TF_confirmed_spls)

##  spls estimates the rank as 7 (by 10-fold cross validation) and selects 36 TFs 
##  out of 106 from which 6 are in the 21 experimentally confirmed ones
##  [1] "CBF1" "MBP1" "MCM1" "STB1" "SWI4" "SWI6"

##===========================================================
##  ccsrrr
##
ff <- cv.ccsrrr(Y, X, nrank=1:10, nbr=3) 
ff$nrank
ff <- ccsrrr(Y, X, nrank=3, nbr=3) 
Abar <- ff$C
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 
Abar[abs(Abar) < 1e-6] <- 0
(TF_ccsrrr <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_ccsrrr)
(TF_confirmed_ccsrrr <- TF_ccsrrr[which(TF_ccsrrr %in% TF_gSCAD)])
length(TF_confirmed_ccsrrr)

##  ccsrrr estimates the rank as 2 or 3. With nrank=3 selects 53 TFs out of 106 from which 
##  15 are in the 21 experimentally confirmed ones
##   [1] "ACE2"  "FKH2"  "GCR1"  "GCR2"  "LEU3"  "MBP1"  "MCM1"  "MET31" "NDD1"  
##       "REB1"  "STB1"  "STE12" "SWI4"  "SWI5"  "SWI6" 

##===========================================================
##  secure
##
control <- secure.control(spU=n/p, spV=1)
cv.fit <- secure.path(Y, X, nrank=4, control=control) 
Abar <- cv.fit$C.est
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 
Abar[abs(Abar) < 1e-6] <- 0

(TF_secure <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_secure)
(TF_confirmed_secure <- TF_secure[which(TF_secure %in% TF_gSCAD)])
length(TF_confirmed_secure)

##  secure estimates the rank as 4 and selects 64 TFs out of 106 from which 
##  15 are in the 21 experimentally confirmed ones
##  [1] "ACE2"  "BAS1"  "FKH2"  "GCN4"  "GCR1"  "MBP1"  "MCM1"  "MET31" "NDD1"  "REB1"  "STB1"  "STE12"
##  [13] "SWI4"  "SWI5"  "SWI6" 

##===========================================================
##  srrr
##
cv.fit <- cv.srrr(Y, X, nfold=10, nrank=4, method = "adglasso")
#cv.fit <- srrr(Y, X, nrank=4, method = "adglasso")
Abar <- cv.fit$coef
Abar[abs(Abar) < 1e-6] <- 0
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 

(TF_srrr <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_srrr)
(TF_confirmed_srrr <- TF_srrr[which(TF_srrr %in% TF_gSCAD)])
length(TF_confirmed_srrr)

##  srrr estimates the rank as 4 and selects 50 TFs out of 106 from which 
##  13 are in the 21 experimentally confirmed ones

##===========================================================
##  remmap
##

lamL1.v <- exp(seq(log(1.5), log(50), length=5))
lamL2.v <- seq(0, 5, length=5)
cv.fit <- remMap.CV(X=X, Y=Y, lamL1.v, lamL2.v, C.m=NULL, fold=5, seed=12345)
pick <- which.min(as.vector(cv.fit$ols.cv))
lamL1.pick <- cv.fit$l.index[1, pick] 
lamL2.pick <- cv.fit$l.index[2, pick]
result <- remMap(X, Y, lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
Abar <- result$phi
Abar[abs(Abar) < 1e-6] <- 0
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 

(TF_remmap <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_remmap)
TF_confirmed_remmap <- TF_remmap[which(TF_remmap %in% TF_gSCAD)]
length(TF_confirmed_remmap)

##  remmap with rank 4 selects 58 TFs out of 106 from which 
##  all 16 are in the 21 experimentally confirmed ones

##===========================================================
##  rssvd
##
cv.fit <- rssvd(Y=Y, X=X, nrank=18)     # nrank=min(n, p, q)
Abar <- coef(cv.fit)
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 
Abar[abs(Abar) < 1e-6] <- 0

(TF_rssvd <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_rssvd)
(TF_confirmed_rssvd <- TF_rssvd[which(TF_rssvd %in% TF_gSCAD)])
length(TF_confirmed_rssvd)

##  RSSVD estimates the rank as 4 and selects 55 TFs out of 106 from which 
##  14 are in the 21 experimentally confirmed ones
##   [1] "ACE2"  "BAS1"  "FKH2"  "GCN4"  "MBP1"  "MCM1"  "MET31" "NDD1"  "REB1" 
##  [10] "STB1"  "STE12" "SWI4"  "SWI5"  "SWI6" 

##===========================================================
##  MRCE
##
lam1.vec <- rev(10^seq(from=-2, to=0, by=0.5))
lam2.vec <- rev(10^seq(from=-2, to=0, by=0.5))

##lam1.vec=rev(10^seq(from=-0.5, to=2, by=0.5))
##lam2.vec=rev(10^seq(from=-2.5, to=0, by=0.5))

cv.fit <- mrce(Y=Y, X=X, lam1.vec=lam1.vec, lam2.vec=lam2.vec, method="cv", standardize=TRUE)

Abar <- cv.fit$Bhat
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 
Abar[abs(Abar) < 1e-6] <- 0

(TF_mrce <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_mrce)
(TF_confirmed_mrce <- TF_mrce[which(TF_mrce %in% TF_gSCAD)])
length(TF_confirmed_mrce)

##  MRCE estimates the rank as 4 and selects 101 TFs out of 106 from which 
##  21 are in the 21 experimentally confirmed ones

##===========================================================
##  sofar
##
cv.fit <- sofar(Y, X, ic.type="GIC", nrank=10, control=list(methodA="adlasso", methodB="adlasso"))

## pred y sofar
P <- cv.fit$U
Q <- cv.fit$V
d <- cv.fit$D
Abar <- if(length(d) == 1) P %*% t(Q) * d  else  P %*% diag(d) %*% t(Q)
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 
Abar[abs(Abar) < 1e-6] <- 0

(TF_sofar <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_sofar)
(TF_confirmed_sofar <- TF_sofar[which(TF_sofar %in% TF_gSCAD)])
length(TF_confirmed_sofar)

##  SOFAR estimates the rank as 4 and selects 50 TFs out of 106 from which 
##  14 are in the 21 experimentally confirmed ones

##==============================================================================
##
##  Table 1
##
library(xtable)
method <- c("SPLS", "SRRR", "SECURE", "REMMAP", "RSSVD", "SOFAR", "CCSRRR")
nn <- matrix(NA, nrow=length(method), ncol=3)
colnames(nn) <- c("Rank", "TFs_selected", "TFs_confirmed") 
dfnn <- data.frame(Method=method, nn)
dfnn[dfnn$Method=="SPLS", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(7, length(TF_spls), length(TF_confirmed_spls))
dfnn[dfnn$Method=="SRRR", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(4, length(TF_srrr), length(TF_confirmed_srrr))
dfnn[dfnn$Method=="SECURE", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(4, length(TF_secure), length(TF_confirmed_secure))
dfnn[dfnn$Method=="REMMAP", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(4, length(TF_remmap), length(TF_confirmed_remmap))
dfnn[dfnn$Method=="RSSVD", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(4, length(TF_rssvd), length(TF_confirmed_rssvd))
dfnn[dfnn$Method=="SOFAR", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(4, length(TF_sofar), length(TF_confirmed_sofar))
dfnn[dfnn$Method=="CCSRRR", c("Rank", "TFs_selected", "TFs_confirmed")] <- c(3, length(TF_ccsrrr), length(TF_confirmed_ccsrrr))
dfnn

print(xtable(dfnn), include.rownames=FALSE)

##===========================================================
##  Figure 1
##
##  The next Figure shows the estimated effects (rows in \hat{C}) of the experimentally 
##  confirmed TFs identified by ccsrrr. As expected, the effects are mostly periodic and 
##  the two cycles are clearly seen (except in "STE12").

ff <- ccsrrr(Y, X, nrank=3, nbr=3) 
Abar <- ff$C
colnames(Abar) <-  paste0("alpha", time)
rownames(Abar) <- TF 
Abar[abs(Abar) < 1e-6] <- 0
(TF_ccsrrr <- TF[which(rowSums(abs(Abar)) != 0)])
length(TF_ccsrrr)
length(which(TF_ccsrrr %in% TF_gSCAD))
(selected <- TF_ccsrrr[which(TF_ccsrrr %in% TF_gSCAD)])

library(reshape2)
library(ggplot2)
a <- cbind.data.frame(time=time, t(Abar[selected,]))
along <- melt(a, id.vars=c("time"))

p <- ggplot(along, aes(time, value)) +
  geom_line(show.legend = FALSE) +
  geom_hline(yintercept=0, linewidth=0.1, linetype="dotted") +
  facet_wrap(~variable, ncol=5, scales = "free_y") +
  ylab("Beta(t)") + xlab("Time") +
  theme_light() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) 
p

ggsave(filename=here::here("Output", "yeast-ccsrrr.png"), device=png, width=10, height=7)

cairo_pdf(filename=here::here("Output", "yeast-ccsrrr.pdf"), width=10, height=7)
print(p)
dev.off()

