#'
#' Cardinality constrained sparse redused rank regression
#'
#' Row- and column-sparse reduced-rank regresssion for a prespecified rank
#'
#'
#' @param Y response matrix
#' @param X predictor matrix
#' @param nrank prespecified rank
#' @param nonzC requested number of nonzero columns. If missing, heuristic will be applied, see parameter \code{nbr}. 
#'  If \code{nonzC==0} or \code{nonzC==ncol(Y)} no column sparseness will be applied.
#' @param nonzR requested number of nonzero rows. If missing, heuristic will be applied, see parameter \code{nbr}. 
#'  If \code{nonzR==0} or \code{nonzR==ncol(X)} no row sparseness will be applied.
#' @param nbr How many bins of the histogram of absolute values to ignore in choosing the threshold
#' @param maxiter maximal number of iterations, default is \code{maxiter=10000}.
#' @param tol tolerance to use for the stoping criterion, by default \code{tol=1e-6}.
#' @return A list of fitting results
#'
#' @examples
#'  data(tobacco)
#'  X <- as.matrix(tobacco[, 4:9])
#'  Y <- as.matrix(tobacco[,1:3]) 
#'  n <- nrow(X)
#'  fit1 <- csrrr(Y, X, nrank=1)
#'  fit1$C
#'  fit1$f

ccsrrr <- function(Y, X, nrank, nonzC, nonzR, nbr=2, maxiter=10000, tol=1e-6) 
{

    ## ALS algorithm to find Q and P(=PD) 
    #+ no need to find D explicitly
    
    n <- nrow(X)    # cases
    p <- ncol(X)    # independent vars
    q <- ncol(Y)    # dependent variables
    r <- nrank    
    stopifnot(nrow(X) == nrow(Y))
    
    ## Scale/whiten the variables (dependent and independent)
    X <- scale(X, center = TRUE, scale = FALSE)/sqrt(n-1)   # whitened data
    Y <- scale(Y, center = TRUE, scale = FALSE)/sqrt(n-1)   # whitened data
    
    XX <- t(X) %*% X
    XY <- t(X) %*% Y
    
    ## Initialization with AB method
    W <- svd(X)
    s <- sum(W$d > tol)
    if(r>s) 
        r <- s

    U <- W$u[, 1:s, drop=FALSE]
    V <- W$v[, 1:s, drop=FALSE] 
    D <- W$d[1:s]

    UY <- t(U) %*% Y
    PD <- V %*% diag(1/D)
    Xinv <- PD %*% UY

    WW <- V %*% UY
    WW <- svd(WW)
    Q <- WW$v[,1:r, drop=FALSE]		
    P <- Xinv %*% Q

    C0 <- P %*% t(Q)
    f_in <- norm(Y - X %*% C0, "F")
    
 
    ## patern of COLUMN sparceness
    if(missing(nonzC)){
        cs <- rowSums(abs(Q))
        rw <- hist(cs, plot = "FALSE")
        th <- rw$breaks[nbr] 
        sparseC <- as.numeric(cs>th)
    }else if(nonzC > 0 && nonzC < q) {
        sparseC <- rep(0, q)
        ind <- order(colSums(C0*C0), decreasing=TRUE) 
        sparseC[ind[1:nonzC]] <- 1
    }else
        sparseC <- rep(1, q)    
 
    ## patern of ROW sparceness 
    if(missing(nonzR)) {
        rs <- rowSums(abs(P))
        rw <- hist(rs, plot = "FALSE")
        th <- rw$breaks[nbr] 
        sparseR <- as.numeric(rs > th)
    }else if(nonzR > 0 && nonzR < p) {
        sparseR <- rep(0, p)
        ind <- order(rowSums(C0*C0), decreasing=TRUE) 
        sparseR[ind[1:nonzR]] <- 1
    }else
        sparseR <- rep(1, p)    
    
    ## Alternating Least Squares algorithm    
    f0 <- norm(Y, 'f')
    Frec <- f0
    iter <- 0
    
    ## Main loop for ALS algorithm
    while(iter <= maxiter) {
    
        ## update PD in P
        if(iter != 0) 
            P <- Xinv %*% Q
        if(length(which(sparseR == 0)) > 0) { 
            Ps <- P * matrix(rep(sparseR, r), p, r)
            P <- Ps
        }
    	
        ## update Q
        W <- svd(t(P) %*% XY)
        Q <- W$v %*% t(W$u)
        if(length(which(sparseC == 0)) > 0) { 
            Qs <- Q * matrix(rep(sparseC, r), q, r)
            Q <- Qs
        } 

        C <- P %*% t(Q)

        f <- norm( Y - X %*% C,'f')   
        Frec <- rbind(Frec, f)
        iter <- iter +1
    
        # Stop criterion
        if(abs(f0-f) <= tol || f0 > f) 
            break
              
        f0 <- f
    }

    ret <- list(P=P, Q=Q, C=C, C0=C0, iter=iter, Frec=Frec, f=f, f_in=f_in, nrank=nrank)
    ret
}

cv.ccsrrr <- function(Y, X, nrank=1:10, kfold=5, nbr=2, ...){
    ##
    ##  From https://stackoverflow.com/questions/7402313/generate-sets-for-cross-validation
    ##
    f_K_fold <- function(Nobs,K=5){
        rs <- runif(Nobs)
        id <- seq(Nobs)[order(rs)]
        k <- as.integer(Nobs*seq(1,K-1)/K)
        k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
        k[,1] <- k[,1]+1
        l <- lapply(seq.int(K),function(x,k,d) 
                    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
                         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
       return(l)
    } 

    n <- nrow(X)    # cases
    p <- ncol(X)    # independent vars
    q <- ncol(Y)    # dependent variables
    r <- nrank    
    stopifnot(nrow(X) == nrow(Y))
    
    mspe <- vector(mode="numeric", length=length(nrank))
    for(j in 1:length(nrank)){
        r <- nrank[j]
    
        ll <- f_K_fold(n, K=kfold)
        mspex <- vector(mode="numeric", length=kfold)
        for(i in 1:kfold) {
            f <- ccsrrr(Y=Y[ll[[i]]$train,], X=X[ll[[i]]$train,], nrank=r, nbr=nbr, ...)
            Abar <- f$C
            Abar[abs(Abar) < 1e-6] <- 0
            pred_y <- X[ll[[i]]$test,] %*% Abar
            mspex[i] <- mean((Y[ll[[i]]$test,] - pred_y)^2)       
        }
        ##  print(mspex)
        mspe[j] <- mean(mspex)
    }
    
    ##  print(mspe)
    
    erank <- nrank[which.min(mspe)]    
    ccsrrr(Y, X, nrank=erank, nbr=nbr, ...)
}
