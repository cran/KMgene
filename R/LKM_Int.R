#' KM for Traits by Time Interaction in Longitudinal GWAS Data
#'
#' This function (LKM_Int) is used to perform Kernel Machine analysis for quantitative traits by time interaction in GWAS longitudinal data. \cr
#' # It considers random intercept and random time
#'
#' @name LKM_Int
#' @aliases LKM_Int
#' @param phenotype A vector of quantitative trait in the analysis (class: vector). The order should match the vector yid. No missing.
#' @param yid A vector of id (class: vector). Although it doesn't have to be sorted, observations from the same subject have to be connected with each other. The repeated id numbers indicate multiple time points for one subject. Make sure it is not a factor. No missing.
#' @param time A vector of time points (class: vector). The order should match the vector yid. No missing.
#' @param genotypes 1st column: gene name; 2nd column: snp name; 3rd-end columns: A matrix of genotypes for each subject (class: data.frame). The order of 3rd-end columns should match unique(yid). Coded as 0, 1, 2 and no missing. This genotype file can be a big file containing all genes or it can be files containing one single gene.
#' @param gid A vector of id mapping to samples in genotype file (class: vector). So the order of samples in gid must be the same as the order in genotypes. Make sure it is not a factor. Although gid doesn't have to be in the same order as yid, it is suggested to make them sorted in the same order in order to make all files easily to be tracked. No missing.
#' @param covariates A matrix of covariates (class: data.frame). The order of rows should match the vector yid. Default NULL. No missing.
#' @param acc Accuracy of numerical integration used in Davies' method. Default 1e-4.
#' @param append.write The name of pvalue output file. Write out p-values in real time. Don't need to wait until all genes are processed to get the final output.
#' @import MASS
#' @import nlme
#' @import CompQuadForm
#' @importFrom utils write.table
#' @importFrom stats as.formula
#' @return output: longitudinal time by marker interaction (LKM_Int) p-value
#' @examples
#' ######################################################################################
#' ### Examples for Marker by Time Interaction in Longitudinal Continuous Traits in GWAS
#' ######################################################################################
#' ###################
#' # Data using KM ###
#' ###################
#' data("LKM_numID_int")
#' pvalue1 <- LKM_Int(phenotype=lkm_int_n_y$y, genotypes=lkm_int_n_gene, time=lkm_int_n_y$time,
#' yid=lkm_int_n_y$id, gid=lkm_int_n_gid$gid, covariates=NULL, append.write="./pvalues.out")
#' @export
LKM_Int <- function(phenotype, time, yid, genotypes, gid, covariates=NULL, acc=1e-4, append.write=NULL){

# Regular checks
n<-length(unique(yid))
m<-length(phenotype)
if(!is.null(covariates)) {
  if(nrow(covariates)!=m) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}
#if(any(unique(yid)!=gid)) warning("id in phenotypes (yid) and id in genotypes (gid) are not in the same order. This would not affect the results")
if(length(gid)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
if(!is.data.frame(genotypes)) stop("genotypes should be a data.frame!")
if(ncol(genotypes)!=n+2) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")

Z.poly<-function(X){
  n<-nrow(X)
  g<-ncol(X)
  maf<-apply(X,2,sum)/(2*n)
  X.new<-matrix(NA,ncol=g,nrow=n)
  for(i in 1:g){
    X.new[,i]<-(X[,i]-2*maf[i])/sqrt(2*maf[i]*(1-maf[i]))
    }
  return(X.new/sqrt(g))
  }

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

genotypes2 <- split(genotypes, genotypes[,1])
genotypes3 <- genotypes4 <- list()
for (k in 1:length(genotypes2)){
 genotypes3[[k]] <- cbind(gid, t(genotypes2[[k]][,-c(1:2)]))
 colnames(genotypes3[[k]])[1] <- "yid"
 genotypes4[[k]] <- merge(yid_order, genotypes3[[k]], by="yid")
 genotypes4[[k]] <- genotypes4[[k]][order(genotypes4[[k]]$order),]
}

# Weight according to beta density function
#if(is.null(weights)){
#W <- list()
#for (k in 1:length(genotypes2)){
# w <- dbeta(colMeans(t(genotypes2[[k]][,-c(1:2)]))/2, 1, 25)
# W[[k]] <- diag(w^2)}}
#else if(!is.null(weights)){
#weights2 <- split(weights, weights[,1])
#W <- list()
#for (k in 1:length(weights2)){
# W[[k]] <- diag(weights2[[k]][,-c(1:2)])
#}}

if(is.null(covariates)){
   X <- cbind(intercept, time)
   X <- as.matrix(X)
   exprs<-paste("y ~ time")}
else if(!is.null(covariates)){
   X <- cbind(intercept, time, as.matrix(covariates))
   X <- as.matrix(X)
   for (i in 1:dim(covariates)[2]) assign(names(covariates)[i], covariates[,i])
   exprs<-paste("y ~ time +", paste(names(covariates),collapse=" + "))}

Q<-eig<-evals<-tmpout<-list()
pvalue<-vector()
for (k in 1:length(genotypes2)){
#class(data$yid) <- "character"
 G = as.matrix(genotypes4[[k]][,-c(1:2)])
 class(G)<-"numeric"
 Z.1 = Z.poly(G)
 Zt.1 = Z.1*time

 mix_model <- lme(as.formula(exprs), random=list(intercept=pdIdent(~Z.1-1), yid=pdSymm(~time)), method="REML", control=lmeControl(opt = "optim"))
 tau <- as.numeric(VarCorr(mix_model)[2,1])
 p <- dim(Z.1)[2]
 error <- as.numeric(VarCorr(mix_model)[p+2+3,1])
 var_int <- as.numeric(VarCorr(mix_model)[p+2+1,1])
 var_time <- as.numeric(VarCorr(mix_model)[p+2+2,1])
 sd_int <- as.numeric(VarCorr(mix_model)[p+2+1,2])
 sd_time <- as.numeric(VarCorr(mix_model)[p+2+2,2])
 corr <- as.numeric(VarCorr(mix_model)[p+2+2,3])
 covariance <- corr*sd_int*sd_time

 K_est = rbind(c(var_int, covariance), c(covariance, var_time))
 Z_prime <- V_prime <- V_prime_inv <- list()
 for (j in 1:length(Z_pre)){
  Z_prime[[j]] <- as.matrix(Z_pre[[j]][,2:3])
  class(Z_prime[[j]]) <- "numeric"
  V_prime[[j]] <- Z_prime[[j]]%*%K_est%*%t(Z_prime[[j]]) + error*diag(1, dim(Z_prime[[j]])[1])
 }
 V_prime <- V_prime[!sapply(V_prime,is.null)] #exclude null list elements
 V_est <- bdiag(V_prime)+tau*Z.1%*%t(Z.1)
 V_est_inv <- solve(V_est)

# The following beta calculation methods give same results
# beta = solve(t(X)%*%V_est_inv%*%X)%*%t(X)%*%V_est_inv%*%y

 beta = summary(mix_model)$tTable[,1]
 res = y - X%*%beta

  P = V_est - X %*% solve(t(X) %*% V_est_inv %*% X) %*% t(X)
  Q[[k]] = t(res) %*% V_est_inv %*% Zt.1 %*% t(Zt.1) %*% V_est_inv %*% res
  eig[[k]] = eigen(t(Zt.1) %*% V_est_inv %*% P %*% V_est_inv %*% Zt.1, symmetric=T, only.values=T)
  evals[[k]] = eig[[k]]$values[eig[[k]]$values>1e-6*eig[[k]]$values[1]]

  tmpout[[k]]<-davies(Q[[k]], evals[[k]], acc=acc)
  pvalue[k] <- tmpout[[k]]$Qq
  if(!is.null(append.write)){
  write.table(t(c(names(genotypes2)[k], signif(pvalue[k], digits=6))), file=append.write, row.names=F, col.names=F, append=T, quote=F) }
}
cbind(names(genotypes2), signif(pvalue, digits=6))
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

