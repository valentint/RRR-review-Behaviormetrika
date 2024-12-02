library(plyr)
library(ggplot2)
library(xtable)

getMySummary2 <- function(output, f1=mean, f2=sd) {
    xf1 <- deparse(substitute(f1))
    xf2 <- deparse(substitute(f2))
    o1 <- ddply(output, c("rho", "b"), .fun=function(x) apply(x[, 3:ncol(x)], 2, f1))
    o2 <- ddply(output, c("rho", "b"), .fun=function(x) apply(x[, 3:ncol(x)], 2, f2))
    o3 <- o1[, 1:2]
    for( i in 3:ncol(o1)){
        o3 <- cbind(o3, round(o1[,i, drop=FALSE], 4), round(o2[,i, drop=FALSE], 4))
        colnames(o3)[ncol(o3)] <- "sd"
    }
    o3
}

here::i_am("getResults2.R")
library(here)

distr <- "norm"     # "t", "U"
case <- 4

switch(case, 
"1" = {
    ##  Case 1: ................
    nsamples <- 100             # number of observations n
    p <- 25                     # dimension of predictors
    q <- 25                     # dimesnion of response
    meth <- c("srrr", "sarrs", "sier", "secure", "mrce", "rssvd", "remmap", "spls", "sofar", "ccsrrr")
    cat("\nDo Case 1, Distribution=", distr, "... \n")    
},
"2" = {
    ##  Case 2: ................
    nsamples <- 30              # number of observations n
    p <- 100                    # dimension of predictors
    q <- 10                     # dimesnion of response
    meth <- c("srrr", "sarrs", "sier", "secure", "mrce", "rssvd", "remmap", "spls", "sofar", "ccsrrr")
    cat("\nDo Case 2, Distribution=", distr, "... \n")    
},
"3" = {
    ##  Case 3:.................
    nsamples <- 30              # number of observations n
    p <- 100                    # dimension of predictors
    q <- 100                    # dimesnion of response
    
    ## Remve MRCE, because it is too slow, more than 600 hours expected
    meth <- c("srrr", "sarrs", "sier", "secure", "rssvd", "spls", "sofar", "ccsrrr")
    cat("\nDo Case 3, Distribution=", distr, "... \n")    
},
"4" = {
##  Case 4:.................
    nsamples <- 30              # number of observations n
    p <- 100                    # dimension of predictors
    q <- 1000                   # dimesnion of response
    
    ## Remve MRCE, because it is too slow
    meth <- c("srrr", "sarrs", "sier", "secure", "spls", "sofar", "ccsrrr")
    cat("\nDo Case 4, Distribution=", distr, "... \n")    
},
{
   stop("Invalid case defined!!!")
}
)

otab <- NULL
ooo <- NULL

for(i in 1:length(meth)){
    imeth <- toupper(meth[i])
    fname <- paste0(imeth, "_", distr, "_n", nsamples, "_p", p, "_q", q, ".csv")
    cat("\nReading file", fname, "...\n")
    o <- read.csv(here::here("Results", fname))
    oo <- cbind.data.frame(o[,1:2], method=imeth, o[,3:ncol(o)])
    ooo <- rbind(ooo, oo)
    
    oo <- getMySummary2(o)
    oo1 <- oo[,1:2]
    oo2 <- oo[,-c(1,2,5,6)]
    oo <- cbind(oo1, method=imeth, oo2)
    otab <- rbind(otab, oo)
} 

getMySummary2(o)
otab <- otab[order(otab$rho, otab$b),]
colnames(otab) <- c("rho", "b", "method", "MSPE", "SD", "Delta", "SD", "Nnz", "SD", "TPR", "SD", "FPR", "SD",  "Time", "SD")
otab

print(xtable(otab, digits=3), include.rownames=FALSE)

## MSPE
p <- ggplot(ooo) +
  geom_boxplot(aes(method, Ys), show.legend = FALSE) +
  facet_grid(vars(b), vars(rho), scales = "free_y") +
  ylab("MSPE") + xlab(element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) 
p

cairo_pdf(filename=here::here("Output", paste0(distr, "_", "case", case, "-MSPE.pdf")), width=14, height=7)
print(p)
dev.off()


## Delta
p <- ggplot(ooo) +
  geom_boxplot(aes(method, As), show.legend = FALSE) +
  facet_grid(vars(b), vars(rho), scales = "free_y") +
  ylab("Delta") + xlab(element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) 
p

cairo_pdf(filename=here::here("Output", paste0(distr, "_", "case", case, "-Delta.pdf")), width=14, height=7)
print(p)
dev.off()

## TPR
p <- ggplot(ooo) +
  geom_boxplot(aes(method, TPRs), show.legend = FALSE) +
  facet_grid(vars(b), vars(rho), scales = "free_y") +
  ylab("TPR") + xlab(element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) 
p

cairo_pdf(filename=here::here("Output", paste0(distr, "_", "case", case, "-TPR.pdf")), width=14, height=7)
print(p)
dev.off()

## FPR
p <- ggplot(ooo) +
  geom_boxplot(aes(method, FPRs), show.legend = FALSE) +
  facet_grid(vars(b), vars(rho), scales = "free_y") +
  ylab("FPR") + xlab(element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) 
p

cairo_pdf(filename=here::here("Output", paste0(distr, "_", "case", case, "-FPR.pdf")), width=14, height=7)
print(p)
dev.off()

## Time
##  Extract only secure and csrrr2
ooo1 <- ooo[ooo$method %in% c("SECURE", "CCSRRR"),]
p <- ggplot(ooo1) +
  geom_boxplot(aes(method, Times), show.legend = FALSE) +
  facet_grid(vars(b), vars(rho), scales = "free_y") +
  ylab("Time (sec.)") + xlab(element_blank()) +
  theme_light() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) 
p

cairo_pdf(filename=here::here("Output", paste0(distr, "_", "case", case, "-Time.pdf")), width=14, height=7)
print(p)
dev.off()

