#' KM for Multiple Quantitative Traits in GWAS Data (fit null model)
#'
#' This function (MKM) is used to perform KM analysis (Tzeng et al. 2012) for multiple (two) quantitative traits in GWAS data. \cr
#' # It considers variances of trait 1 and 2, and covariance between trait 1 and 2
#'
#' @name MKM_Null_Model
#' @aliases MKM_Null_Model
#' @param phenotype A vector of quantitative trait in the analysis (class: vector). The order should match the vector yid. No missing.
#' @param yid A vector of id (class: vector). Although it doesn't have to be sorted, observations from the same subject have to be connected with each other. The repeated id numbers indicate multiple traits for one subject. Make sure it is not a factor. No missing.
#' @param trait A vector of multivariate traits (class: vector). The order should match the vector yid. No missing.
#' @param covariates A matrix of covariates (class: data.frame). The order of rows should match the vector yid. Default NULL. No missing.
#' @param eq.cov.effect Whether assume equal covariates effects on different traits (Default=False).
#' @import MASS
#' @import nlme
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @return output: object as input for MKM
#' @export
MKM_Null_Model <- function(phenotype, trait, yid, covariates=NULL, eq.cov.effect=F){

# Regular checks
n<-length(unique(yid))
m<-length(phenotype)
if(!is.null(covariates)) {
  if(nrow(covariates)!=m) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

## block diagonal matrix function ##
bdiag2 <- function(x){
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
dummy <- model.matrix( ~ trait - 1)
design <- as.data.frame(cbind(yid, intercept))
y <- as.matrix(phenotype)

design$yid <- factor(design$yid, levels=unique(design$yid))
Z_pre=split(as.data.frame(design), design$yid)


# if equal effects ###################
 ran <- paste("random=~-1+", paste(colnames(dummy),collapse=" + "), "|yid", sep="")
if(eq.cov.effect){
 if(is.null(covariates)){
   X <- intercept
   X <- as.matrix(X)
   data <- data.frame(y=y, yid=yid, dummy)
   exprs<-paste("y ~ 1")}
 else if(!is.null(covariates)){
   X <- cbind(intercept, as.matrix(covariates))
   X <- as.matrix(X)
   data <- data.frame(y=y, yid=yid, dummy, covariates)
   exprs<-paste("y ~ ", paste(names(covariates),collapse=" + "), sep="")}
}

if(!eq.cov.effect){
 if(is.null(covariates)){
   X <- cbind(intercept,dummy[,1])
   X <- as.matrix(X)
   data <- data.frame(y=y, yid=yid, dummy)
   exprs<-paste("y ~ ", colnames(dummy)[1], sep="")}
 else if(!is.null(covariates)){
   covariates1=as.data.frame(covariates)*dummy[,1]
   names(covariates1)=paste(names(covariates), colnames(dummy)[1], sep="_")
   covariates2=as.data.frame(covariates)*dummy[,2]
   names(covariates2)=paste(names(covariates), colnames(dummy)[2], sep="_")
   covariates_all=cbind(covariates1,covariates2)
   X <- cbind(intercept, dummy[,1], as.matrix(covariates_all))
   X <- as.matrix(X)
   data <- data.frame(y=y, yid=yid, dummy, covariates_all)
   exprs<-paste("y ~ ", colnames(dummy)[1], " + ", paste(names(covariates_all),collapse=" + "), sep="")}
}

 mix_model <- lme(as.formula(exprs), as.formula(ran), data=data, method="REML", control=lmeControl(opt = "optim"))
 beta = summary(mix_model)$tTable[,1]
 sigma2 <- as.numeric(VarCorr(mix_model)[3])
 var_ph1 <- as.numeric(VarCorr(mix_model)[1])
 var_ph2 <- as.numeric(VarCorr(mix_model)[2])
 sd_ph1 <- as.numeric(VarCorr(mix_model)[4])
 sd_ph2 <- as.numeric(VarCorr(mix_model)[5])
 corr <- as.numeric(VarCorr(mix_model)[8])
 covariance <- corr*sd_ph1*sd_ph2

 K_est = rbind(c(var_ph1, covariance), c(covariance, var_ph2))
 Z_prime <- V_prime <- V_prime_inv <- list()
 for (j in 1:length(Z_pre)){
  Z_prime[[j]] <- as.matrix(Z_pre[[j]][,2])
  class(Z_prime[[j]]) <- "numeric"
  V_prime[[j]] <- K_est + sigma2*diag(1, dim(Z_prime[[j]])[1])
  V_prime_inv[[j]] <- solve(V_prime[[j]])
 }
 V_prime <- V_prime[!sapply(V_prime,is.null)] #exclude null list elements
 V_prime_inv <- V_prime_inv[!sapply(V_prime_inv,is.null)]
 V_est <- bdiag2(V_prime)
 V_est_inv <- bdiag2(V_prime_inv)

# The following beta calculation methods give same results
# beta = solve(t(X)%*%V_est_inv%*%X)%*%t(X)%*%V_est_inv%*%y

res = y - X%*%beta
 P = V_est - X %*% solve(t(X) %*% V_est_inv %*% X) %*% t(X)
list("res"=res, "V_est"=V_est, "V_est_inv"=V_est_inv, "n"=n, "X"=X, "yid"=yid, "yid_order"=yid_order, "dummy"=dummy, "P"=P)
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

