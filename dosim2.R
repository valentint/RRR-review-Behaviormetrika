##  Simulation framework for sparse redused rank rgression following [1] and [2]:
##
##  [1] Bunea, F., She, Y., and Wegkamp, M. H. (2012). Joint variable and rank selection 
##      for parsimonious estimation of high-dimensional matrices. The Annals of 
##      Statistics, 40(5):2359–2388.
##  [2] Hilafu, H., Safo, S. E., and Haine, L. (2020). Sparse reduced-rank regression 
##      for integrating omics data. BMC Bioinformatics, 21(283):1–17.
##
##  Methods: 
##      spls, SRRR, SARRS, SeCURE, RSSVD, remMap, MRCE, SIER, SOFAR, CCSRRR (ours)
##
##  Cases:
##  1. n > p = q (n = 100, p = 25, q = 25, s = 15, r = 5),      b = 0.2, 0.4,   rho = 0.1, 0.5, 0.9
##  2. q < n < p (n = 30, p = 100, q = 10, s = 15, r = 2),      b = 0.5, 1,     rho = 0.1, 0.5, 0.9
##  3. n < p = q (n = 30, p = 100, q = 100, s = 15, r = 2),     b = 0.5, 1,     rho = 0.1, 0.5, 0.9
##  4. n < p < q (n = 30, p = 100, q = 1000, s = 15, r = 5),    b = 0.5, 1,     rho = 0.1, 0.5, 0.9
##
##  Meassures:
##  Delta = ||C - \hat{C}||^2/(pq) = mean((C-\hat{C})^2) = mean((A-Abar)^2)
##  MSPE = ||Yt - \hat{Yt}||^2/(ntq) = mean((Yt-\hat{Yt})^2) + mean((y_test - pred_y)^2)
##  TPR: true positive rate, the ratio of truly important variables that the method selects as important. 
##  FPR: false positive rate, the ratio of unimportant predictors that the method selects as important.
##      TPR values close to one and FPR values close to zero indicate a better variable selection performance.
##

## Create the predictor covariance matrix
getMyCov <- function(p, rho) {
    fun <- function(i, j) rho^(abs(i-j))
    rows <- 1:p
    cols <- 1:p

    outer(rows, cols, FUN=fun)
}

## Create coefficient matrix A
getMyCoef <- function(p, q, r, s, b) {
    b0 <- matrix(rnorm(s*r, mean=0, sd=1), nrow=s, ncol=r)
    b1 <- matrix(rnorm(r*q, mean=0, sd=1), nrow=r, ncol=q)
    a1 <- b*b0 %*% b1

    rbind(a1, matrix(0, nrow=p-s, ncol=q))
}

getMyData <- function(n=100, p=25, q=25, s=15, r=5, ntest=0, rho=0.1, b=0.2, distr=c("norm", "t", "U")){

    distr <- match.arg(distr)
    N <- n + ntest

    Covmatrix <- getMyCov(p, rho)
    A <- getMyCoef(p, q, r, s, b)

    x <- mvrnorm(N, rep(0,p), Sigma=Covmatrix)
    if(distr == "norm")
        y <- x %*% A + matrix(rnorm(N*q, mean=0, sd=1), nrow=N, ncol=q)
    else if(distr == "t")
         y <- x %*% A + sqrt(3/5) * matrix(rt(N*q, df=5), nrow=N, ncol=q)
    else if(distr == "U")
        y <- x %*% A + matrix(runif(N*q, min=-1, max=1), nrow=N, ncol=q) + 
                       matrix(runif(N*q, min=-1, max=1), nrow=N, ncol=q) + 
                       matrix(runif(N*q, min=-1, max=1), nrow=N, ncol=q)        
    else
        stop(paste("Unknown error distribution: ", distr))
            
    ## Create training and test set
    train_index <- 1:n
    test_index <- setdiff(1:nrow(x), train_index)
    x_train <- x[train_index,]
    y_train <- y[train_index,]
    x_test <- x[test_index,]
    y_test <- y[test_index,]
    
    if(ntest > 0)
        list(x_train=x_train, y_train=y_train, x_test=x_test, y_test=y_test, A=A)
    else
        list(X=x_train, Y=y_train, A=A)
}

getMySummary <- function(output, f1=mean, f2=sd) {
    xf1 <- deparse(substitute(f1))
    xf2 <- deparse(substitute(f2))
    ddply(output, c("rho", "b"), .fun=function(x) {
        cbind.data.frame(x=c(xf1, xf2), 
        rbind(apply(x[, 3:ncol(x)], 2, f1), apply(x[,3:ncol(x)], 2, f2)))})
}

library(conflicted)     # avoid masking of functions
library(MASS)
library(plyr)

library(spls)
library(SiER)
library(remMap)
library(gglasso)
library(rrpack)
library(MRCE)
library(secure)

here::i_am("dosim2.R")
library(here)

source("ccsrrr.R")

##### Set case and method:
set.seed(12345)

case <- 1
rep <- 50
ntest <- 10000

method <- "ccsrrr"          # c("spls", "srrr", "sarrs", "secure", "remmap", "rssvd", "mrce", "sier", "sofar", "ccsrrr")
distr <- "norm"             # c("norm", "t", "U")

switch(case, 
"1" = {
    ##  Case 1: ................
    nsamples <- 100             # number of observations n
    n <- ntest + nsamples       # number of train plus test (10000) observations
    p <- 25                     # dimension of predictors
    q <- 25                     # dimesnion of response
    s <- 15                     # number of important variables (the first s variables in the predictor matrix X are the important ones)
    r <- 5                      # rank
    cat("\nRunning Case 1, Distribution=", distr, " and method = ", method, " ... \n")    
},
"2" = {
    ##  Case 2: ................
    nsamples <- 30              # number of observations n
    n <- ntest + nsamples       # number of train plus test (10000) observations
    p <- 100                    # dimension of predictors
    q <- 10                     # dimesnion of response
    s <- 15                     # number of important variables (the first s variables in the predictor matrix X are the important ones)
    r <- 2                      # rank
    cat("\nRunning Case 2, Distribution=", distr, " and method = ", method, " ... \n")    
},
"3" = {
    ##  Case 3:.................
    nsamples <- 30              # number of observations n
    n <- ntest + nsamples       # number of train plus test (10000) observations
    p <- 100                    # dimension of predictors
    q <- 100                    # dimesnion of response
    s <- 15                     # number of important variables (the first s variables in the predictor matrix X are the important ones)
    r <- 2                      # rank
    cat("\nRunning Case 3, Distribution=", distr, " and method = ", method, " ... \n")    
},
"4" = {
##  Case 4:.................
    nsamples <- 30              # number of observations n
    n <- ntest + nsamples       # number of train plus test (10000) observations
    p <- 100                    # dimension of predictors
    q <- 1000                   # dimesnion of response
    s <- 15                     # number of important variables (the first s variables in the predictor matrix X are the important ones)
    r <- 5                      # rank
    cat("\nRunning Case 4, Distribution=", distr, " and method = ", method, " ... \n")    
},
{
   stop("Invalid case defined!!!")
}
)

rhooo <- c(0.1, 0.5, 0.9)   # rho to define the predictor covariance matrix)
boo   <- c(0.2, 0.4)        # signal-to-noise ratio
if(nsamples < p) 
    boo <- c(0.5, 1)
output <- NULL


start_time_0 <- proc.time()
for(rhoo in 1:length(rhooo)) {
  rho <- rhooo[rhoo]
  
  ## Create the predictor covariance matrix
  Covmatrix <- getMyCov(p, rho)
  
  for(bo in 1:2) {
    b <- boo[bo]  
    
    ## Create coefficient matrix A
    A <- getMyCoef(p, q, r, s, b)
    
    ## Simulation output initialize
    Ys <- XAs <- As <- Nnz <- TPRs <- FPRs <- Times <- vector(mode="numeric", length=rep)
    
    for(i in 1:rep) {
        cat("\nrho=", rho, "b=", b, "i=", i, "...\n")

        # Generate x and y
        x <- mvrnorm(n, rep(0,p), Sigma=Covmatrix)
        if(distr == "norm")
            y <- x %*% A + matrix(rnorm(n*q, mean=0, sd=1), nrow=n, ncol=q)
        else if(distr == "t")
             y <- x %*% A + sqrt(3/5) * matrix(rt(n*q, df=5), nrow=n, ncol=q)
        else if(distr == "U")
            y <- x %*% A + matrix(runif(n*q, min=-1, max=1), nrow=n, ncol=q) + 
                           matrix(runif(n*q, min=-1, max=1), nrow=n, ncol=q) + 
                           matrix(runif(n*q, min=-1, max=1), nrow=n, ncol=q)        
        else
            stop(paste("Unknown error distribution: ", distr))
            
        ## Create training and test set
        train_index <- 1:nsamples
        test_index <- setdiff(1:nrow(x), train_index)
        x_train <- x[train_index,]
        y_train <- y[train_index,]
        x_test <- x[test_index,]
        y_test <- y[test_index,]

        ##data <- getMyData(n, p, q, s, r, ntest, rho=rho, b=b, distr=distr)
        ##A <- data$A
        ##x_train <- data$x_train
        ##y_train <- data$y_train
        ##x_test <- data$x_test
        ##y_test <- data$y_test

        start_time <- proc.time()

        if(method == "spls") {
            ## Fit CV spls
            tt <- system.time(
                {   cv.fit <- cv.spls(x=x_train, y=y_train, fold=5, K=c(1:10), eta = seq(0.1,0.9,0.1), scale.x=FALSE, plot.it=FALSE);
                    fit <- spls(x=x_train, y=y_train, K=cv.fit$K.opt, eta=cv.fit$eta.opt, kappa=0.5, select="pls2", fit="simpls", scale.x=FALSE)
                })
            
            ## Predict spls
            Abar <- coef(fit)
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- predict(fit, newx=x_test, type="fit")
            XA <- x_test %*% A
        }else if(method == "sier") {
            ## fit cv.SiER
            tt <- system.time(cv.fit <- cv.SiER(x_train, y_train, K.cv=5, upper.comp=r))
            
            ## pred y sier
            pred_y <- pred.SiER(cv.fit, x_test)
            Abar <- getcoef.SiER(cv.fit)$beta
            Abar[abs(Abar) < 1e-6] <- 0
            XA <- x_test%*%A
        }else if(method == "remmap") {
            
            ## fit remMap.CV
            tt <- system.time({
                lamL1.v <- exp(seq(log(1.5), log(50), length=5))
                lamL2.v <- seq(0, 5, length=5)
                cv.fit <- remMap.CV(X=x_train, Y=y_train, lamL1.v, lamL2.v, C.m=NULL, fold=5, seed=12345)
                pick <- which.min(as.vector(cv.fit$ols.cv))
                lamL1.pick <- cv.fit$l.index[1, pick] 
                lamL2.pick <- cv.fit$l.index[2, pick]
                result <- remMap(x_train, y_train, lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL)
            })
                        
            ## pred y remMap
            Abar <- result$phi
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A
        }else if(method == "sarrs") {
           
            ## create initial V(0)
            tt <- system.time({
                V_0 <- svd(y_train, nu=r, nv=r)$v
                ##sss = svd(A,nv=r)
                ##V_0 = sss$v
                
                ## group penalized regression.
                ybar <- as.vector(y_train %*% V_0)
                xbar <- kronecker(diag(1, r), x_train)
                group <- rep(1:p, r)
                
                cv.fit <- cv.gglasso(xbar, ybar, group=group, pred.loss="L2", nfolds=5)
                lambda <- cv.fit$lambda.1se
                fit <- gglasso(xbar, ybar, group=group, lambda=lambda)
            })
                        
            ## pred y ssr
            B_1 <- matrix(fit$beta, p, r, byrow=FALSE)
            Abar <- B_1 %*% t(V_0)
            
            ## thresholding step
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A
        }else if(method == "srrr") {
            ## fit cv.srrr
            tt <- system.time(cv.fit <- cv.srrr(y_train, x_train, nfold=5, nrank=r, method="adglasso"))
            
            ## pred y srrr
            Abar <- cv.fit$coef
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A

        }else if(method == "secure") {

            ## fit cv.secure
            ##  secure model with adaptive lasso penalty, without orthogonality constraint
            control <- secure.control(elnetAlpha=1, spU=3.3, spV=1)
            ##  This is faster but worse sparseness performance (larger MSPE)
            if(case == 4)
                control <- secure.control(elnetAlpha=1, spU=0.25, spV=1)
            tt <- system.time(cv.fit <- secure.path(y_train, x_train, nrank=r+2, nlambda=100, control=control)
            )
            
            ## pred y secure
            Abar <- cv.fit$C.est
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A
            
        }else if(method == "mrce") {
            
            ## fit mrce.CV           
            tt <- system.time({
                ## These lambda vectors work with n>pbut not with n <= p. Use theese below
                ##  (recommended by Adam)...
                ##  lam1.vec <- rev(10^seq(from=-2, to=0, by=0.5))
                ##  lam2.vec <- rev(10^seq(from=-2, to=0, by=0.5))
                
                ## These lambdas work for n < p
                lam1.vec=rev(10^seq(from=-0.5, to=2, by=0.5))
                lam2.vec=rev(10^seq(from=-2.5, to=0, by=0.5))
                cv.fit <- mrce(Y=y_train, X=x_train, lam1.vec=lam1.vec, lam2.vec=lam2.vec, method="cv")
            })
                        
            ## pred y mrce
            Abar <- cv.fit$Bhat
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A
            
        }else if(method == "rssvd") {
            ## fit rssvd.CV
            tt <- system.time(cv.fit <- rssvd(Y=y_train, X=x_train, nrank=min(n, p, q)))
            
            ## pred y rssvd
            Abar <- coef(cv.fit)
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A

        }else if(method == "ccsrrr") {
            
            ## fit csrrr
            ## Set nbr (number of breaks) depending on n and p
            nbr <- if(nsamples > p) 2 else 3
            tt <- system.time(cv.fit <- ccsrrr(Y=y_train, X=x_train, nrank=r, nbr=nbr))
            
            ## pred y csrrr
            Abar <- cv.fit$C
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A

        }else if(method == "sofar") {
            ## fit cv.sofar
            tt <- system.time(cv.fit <- sofar(y_train, x_train, ic.type="GIC", nrank=r+2, control=list(methodA="adlasso", methodB="adlasso")))
            
            ## pred y sofar
            P <- cv.fit$U
            Q <- cv.fit$V
            d <- cv.fit$D
            
            Abar <- if(length(d) == 1) P %*% t(Q) * d  else  P %*% diag(d) %*% t(Q)
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- x_test %*% Abar
            XA <- x_test %*% A
        }else
            stop(paste("Unknown method:", method))

        end_time <- proc.time()
        
        ## Store result for each simulation
        Ys[i] <- mean((y_test - pred_y)^2)
        XAs[i] <- mean((XA - pred_y)^2)
        As[i] <- mean((A - Abar)^2)
        Nnz[i] <- length(which(Abar != 0))/(p*q)
        TPRs[i] <- mean(apply(Abar[1:s,] != 0, 1, max))
        FPRs[i] <- mean(apply(Abar[(s+1):p,] != 0, 1, max))
        Times[i] <- tt[3]   #(end_time-start_time)[3]
        cat("\nrho=", rho, "b=", b, "i=", i, "elapsed time: ", Times[i], "\n")
    }
    
    results <- cbind.data.frame(rho=rho, b=b, Ys, XAs, As, Nnz, TPRs, FPRs, Times)
    output <- rbind.data.frame(output, results)
    }
}

end_time_0 <- proc.time()
cat("\nTotal ealpsed time: ", (end_time_0-start_time_0)[3], "\n")

head(output)
dim(output)

### Show summary of the result 
rbind(apply(output[1:rep,], 2, mean), apply(output[1:rep,], 2, sd))
rbind(apply(output[(rep+1):(2*rep),], 2, mean), apply(output[(rep+1):(2*rep),], 2, sd))
rbind(apply(output[(2*rep+1): (3*rep),], 2, mean), apply(output[(2*rep+1): (3*rep),], 2, sd))
rbind(apply(output[(3*rep+1):(4*rep),], 2, mean), apply(output[(3*rep+1):(4*rep),], 2, sd))
rbind(apply(output[(4*rep+1):(5*rep),], 2, mean), apply(output[(4*rep+1):(5*rep),], 2, sd))
rbind(apply(output[(5*rep+1):(6*rep),], 2, mean), apply(output[(5*rep+1):(6*rep),], 2, sd))

getMySummary(output)

fname <- paste0(method, "_", distr, "_n", nsamples, "_p", p, "_q", q, ".csv")
write.csv(output, file=here::here("Results", fname), row.names=FALSE)

