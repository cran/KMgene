#' Optimal KM for Quantitative Traits in Longitudinal GWAS Data (fit null model)
#'
#' This function (LKMO) is used to perform optimal KM analysis (Yan et al., 2016) for quantitative traits in GWAS longitudinal data. \cr
#' # It considers random intercept and random time
#'
#' @name LKMO_Null_Model
#' @aliases LKMO_Null_Model
#' @param phenotype A vector of quantitative trait in the analysis (class: vector). The order should match the vector yid. No missing.
#' @param yid A vector of id (class: vector). Although it doesn't have to be sorted, observations from the same subject have to be connected with each other. The repeated id numbers indicate multiple time points for one subject. Make sure it is not a factor. No missing.
#' @param time A vector of time points (class: vector). The order should match the vector yid. No missing.
#' @param covariates A matrix of covariates (class: data.frame). The order of rows should match the vector yid. Default NULL. No missing.
#' @import MASS
#' @import nlme
#' @importFrom stats as.formula
#' @return output: object as input for LKMO
#' @export
LKMO_Null_Model <- function(phenotype, time, yid, covariates=NULL){

# Regular checks
n<-length(unique(yid))
m<-length(phenotype)
if(!is.null(covariates)) {
  if(nrow(covariates)!=m) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

## block diagonal matrix function ##
bdiag <- function(x){
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
}

order <- seq(1, length(as.vector(yid)))
yid_order <- cbind.data.frame(yid, order)

intercept <- rep(1,length(yid))
design <- as.data.frame(cbind(yid, intercept, time))
y <- as.matrix(phenotype)

design$yid <- factor(design$yid, levels=unique(design$yid))
Z_pre=split(as.data.frame(design), design$yid)


if(is.null(covariates)){
   X <- cbind(intercept, time)
   X <- as.matrix(X)
   exprs<-paste("y ~ time")}
else if(!is.null(covariates)){
   X <- cbind(intercept, time, as.matrix(covariates))
   X <- as.matrix(X)
   exprs<-paste("y ~ time +", paste(names(covariates),collapse=" + "))}

 data_pre <- data.frame(yid=yid,y=y, time=time)
 if(is.null(covariates)) data <- data_pre
 else if(!is.null(covariates)) data <- cbind(data_pre, covariates)

#class(data$yid) <- "character"
mix_model <- lme(as.formula(exprs), random= ~ time | yid, data=data, method="REML", control=lmeControl(opt = "optim"))
error <- as.numeric(VarCorr(mix_model)[3])
var_int <- as.numeric(VarCorr(mix_model)[1])
var_time <- as.numeric(VarCorr(mix_model)[2])
sd_int <- as.numeric(VarCorr(mix_model)[4])
sd_time <- as.numeric(VarCorr(mix_model)[5])
corr <- as.numeric(VarCorr(mix_model)[8])
covariance <- corr*sd_int*sd_time

K_est = rbind(c(var_int, covariance), c(covariance, var_time))
Z_prime <- V_prime <- V_prime_inv <- list()
for (j in 1:length(Z_pre)){
 Z_prime[[j]] <- as.matrix(Z_pre[[j]][,2:3])
 class(Z_prime[[j]]) <- "numeric"
 V_prime[[j]] <- Z_prime[[j]]%*%K_est%*%t(Z_prime[[j]]) + error*diag(1, dim(Z_prime[[j]])[1])
 V_prime_inv[[j]] <- solve(V_prime[[j]])
}
V_prime <- V_prime[!sapply(V_prime,is.null)] #exclude null list elements
V_prime_inv <- V_prime_inv[!sapply(V_prime_inv,is.null)]
V_est <- bdiag(V_prime)
V_est_inv <- bdiag(V_prime_inv)

# The following beta calculation methods give same results
# beta = solve(t(X)%*%V_est_inv%*%X)%*%t(X)%*%V_est_inv%*%y

beta = summary(mix_model)$tTable[,1]
res = y - X%*%beta
P = V_est_inv - V_est_inv %*% X %*% solve(t(X) %*% V_est_inv %*% X) %*% t(X) %*% V_est_inv
svdob<-eigen(P, symmetric=TRUE)
Phalf<-svdob$vectors%*%diag(sqrt(round(svdob$values, 6)))%*%t(svdob$vectors)

list("res"=res, "V_est"=V_est, "V_est_inv"=V_est_inv, "n"=n, "X"=X, "yid"=yid, "yid_order"=yid_order, "P"=P, "Phalf"=Phalf)
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

