#' Optimal KM for Quantitative Traits in Multivariate Family GWAS Data (fit null model)
#'
#' This function (MFKMO) is used to perform optimal KM analysis for quantitative traits in GWAS multivariate family data. \cr
#' # It takes familial correlation as a kinship matrix
#'
#' @name MFKMO_Null_Model
#' @aliases MFKMO_Null_Model
#' @param phenotype A vector of quantitative trait in the analysis (class: vector). The order should match the vector yid. No missing.
#' @param yid A vector of id (class: vector). Although it doesn't have to be sorted, observations from the same subject have to be connected with each other. The repeated id numbers indicate multiple time points for one subject. Make sure it is not a factor. No missing.
#' @param gid A vector of id mapping to samples in genotype file (class: vector). So the order of samples in gid must be the same as the order in genotypes. Make sure it is not a factor. Although gid doesn't have to be in the same order as yid, it is suggested to make them sorted in the same order in order to make all files easily to be tracked. No missing.
#' @param fa A vector of father id (class: vector). The father id indicates the father of each subject. If this subject has no father in this data, the value is set to "NA". Make sure it is not factor.
#' @param mo A vector of mother id (class: vector). The mother id indicates the mother of each subject. If this subject has no mother in this data, the value is set to "NA". Make sure it is not factor.
#' @param trait A vector of multivariate traits (class: vector). The order should match the vector yid. No missing.
#' @param covariates A matrix of covariates (class: data.frame). The order of rows should match the vector yid. Default NULL. No missing.
#' @param Ninitial The number of times to try initial values. The default is 10 times. If Ninitial=1, the initial value "cor" is always equal to correlation(trait1|covariates, trait2|covariates). One should try multiple initial values in order to find max log-likelihood. This could be time consuming, depends on the sample size. The good thing is that null model only needs to be fitted once for the whole genome, so it's worth trying many initial values.
#' @param LL.output Output all tried initial values and corresponding log-likelihoods. The initial value with max log-likelihood is used in the algorithm and it can be used for replication. The output file can be renamed.
#' @param cor Initial value. By default, it's not given, the program tries to find the best initial value. Once it's given, the program uses it as the only initial value. This is useful when one already knows the initial value corresponding to max log-likelihood.
#' @param eq.cov.effect Whether assume equal covariates effects on different traits (Default=False).
#' @param method The optimization method used in null model. The default method is an implementation of that of Nelder and Mead (1965), that uses only function values and is robust but relatively slow. Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each variable can be given a lower and/or upper bound.
#' @import kinship2
#' @importFrom utils write.table
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @importFrom stats resid
#' @importFrom stats lm
#' @importFrom stats runif
#' @importFrom stats optim
#' @return output: object as input for MFKMO
#' @export
MFKMO_Null_Model <- function(phenotype, trait, yid, gid, fa, mo, covariates=NULL, Ninitial=10, method="Nelder-Mead", LL.output='./LogLikelihood.txt', cor=NULL, eq.cov.effect=F){

if (method=="Nelder-Mead"){
reml.lik<-function(theta,y,X,K,I){
 #theta[1]: variance for family pheno1
 #theta[2]: variance for family pheno2
 #theta[3]: covariance for family
 #theta[4]: variance for family pheno1
 #theta[5]: variance for family pheno2
 #theta[6]: covariance for family
 K1<-exp(theta[1]) + abs(theta[3])
 K2<-exp(theta[2]) + abs(theta[3])
 K12<-theta[3]
 e1<-exp(theta[4]) + abs(theta[6])
 e2<-exp(theta[5]) + abs(theta[6])
 e12<-theta[6]
 K_est = matrix(c(K1,K12,K12,K2),nrow=2)
 e_est = matrix(c(e1,e12,e12,e2),nrow=2)
 V_K <- kronecker(K, K_est)
 V_e <- kronecker(I, e_est)
 V_est <- V_K + V_e
 V_est_inv <- solve(V_est)
 determinant1 = determinant(V_est, logarithm=T)[1:2]
 log_det_V_est <- determinant1$modulus[1] * determinant1$sign
 XV_estX=t(X)%*%V_est_inv%*%X
 determinant2 = determinant(XV_estX, logarithm=T)[1:2]
 log_det_XV_estX <- determinant2$modulus[1] * determinant2$sign
 beta=solve(XV_estX)%*%t(X)%*%V_est_inv%*%y
 mu_hat=X%*%beta
 logl<- log_det_V_est + log_det_XV_estX + t(y-mu_hat)%*%V_est_inv%*%(y-mu_hat)
 return(as.numeric(logl))
}

# Regular checks
n<-length(unique(yid))
m<-length(phenotype)
#if(any(unique(yid)!=gid)) warning("id in phenotypes (yid) and id in genotypes (gid) are not in the same order. This would not affect the results")
if(length(gid)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
if(!is.null(covariates)) {
  if(nrow(covariates)!=m) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

order <- seq(1, length(as.vector(yid)))
yid_order <- cbind.data.frame(yid, order)

type3 <- cbind(gid, fa, mo)
colnames(type3)[1] <- "yid"
type4 <- merge(yid_order, type3, by="yid")
type4 <- type4[order(type4$order),]

relation <- unique(cbind.data.frame(type4$yid, type4$fa, type4$mo))
relation <- as.matrix(relation)
K <- 2*kinship(relation[,1], relation[,2], relation[,3])
I <- diag(dim(relation)[1])

intercept <- rep(1,length(yid))
dummy <- model.matrix( ~ trait - 1)
y <- as.matrix(phenotype)

if(is.null(covariates)){
   X <- intercept
   X <- as.matrix(X)
   data <- data.frame(y=y, dummy)
   exprs<-paste("y ~ 1")}
else if(!is.null(covariates)){
   X <- cbind(intercept, as.matrix(covariates))
   X <- as.matrix(X)
   data <- data.frame(y=y, dummy, covariates)
   exprs<-paste("y ~ ", paste(names(covariates),collapse=" + "), sep="")}

 if(!eq.cov.effect){
 X1=X*rep(c(1,0), length(yid)/2)
 X2=X*rep(c(0,1), length(yid)/2)
 X=cbind(X1, X2)}

 data1<-data[data[,2]==1,]
 data2<-data[data[,3]==1,]

 if (Ninitial>1 & is.null(cor)){
  cor<-LL<-vector()
  all<-theta<-list()
  guess<-cor(resid(lm(as.formula(exprs), data=data1)), resid(lm(as.formula(exprs), data=data2)))
  if (sign(guess)==1)  cor<-runif(Ninitial, min=0, max=1)
  if (sign(guess)==-1)  cor<-runif(Ninitial, min=-1, max=0)
  cor[1]<-guess
  all[[1]] = optim(c(log(1-abs(cor[1])), log(1-abs(cor[1])), cor[1], log(1-abs(cor[1])), log(1-abs(cor[1])), cor[1]), reml.lik, y=y, X=X, K=K, I=I)
  theta[[1]] = all[[1]]$par
  LL[1] = -0.5*all[[1]]$value
  write.table(t(c(cor[1], LL[1])), file=LL.output, append=T, row.names=F, col.names=F, quote=F)
  for (i in 2:Ninitial){
   all[[i]] = optim(c(log(1-abs(cor[i])), log(1-abs(cor[i])), cor[i], log(1-abs(cor[i])), log(1-abs(cor[i])), cor[i]), reml.lik, y=y, X=X, K=K, I=I)
   theta[[i]] = all[[i]]$par
   LL[i] = -0.5*all[[i]]$value
   write.table(t(c(cor[i], LL[i])), file=LL.output, append=T, row.names=F, col.names=F, quote=F)
  }
  theta = theta[[which(LL==max(LL))]]
  LL = max(LL)
 }

 else if (Ninitial==1 & is.null(cor)){
  cor<-cor(resid(lm(as.formula(exprs), data=data1)), resid(lm(as.formula(exprs), data=data2)))
  all = optim(c(log(1-abs(cor)), log(1-abs(cor)), cor, log(1-abs(cor)), log(1-abs(cor)), cor), reml.lik, y=y, X=X, K=K, I=I)
  theta = all$par
  LL = -0.5*all$value
  write.table(t(c(cor, LL)), file=LL.output, row.names=F, col.names=F, quote=F)
 }
 else if (!is.null(cor)){
  all = optim(c(log(1-abs(cor)), log(1-abs(cor)), cor, log(1-abs(cor)), log(1-abs(cor)), cor), reml.lik, y=y, X=X, K=K, I=I)
  theta = all$par
  LL = -0.5*all$value
  write.table(t(c(cor, LL)), file=LL.output, row.names=F, col.names=F, quote=F)
 }

 K1<-exp(theta[1]) + abs(theta[3])
 K2<-exp(theta[2]) + abs(theta[3])
 K12<-theta[3]
 e1<-exp(theta[4]) + abs(theta[6])
 e2<-exp(theta[5]) + abs(theta[6])
 e12<-theta[6]
 K_est = matrix(c(K1,K12,K12,K2),nrow=2)
 e_est = matrix(c(e1,e12,e12,e2),nrow=2)
 V_K <- kronecker(K, K_est)
 V_e <- kronecker(I, e_est)
 V_est <- V_K + V_e
 V_est_inv <- solve(V_est)
 XV_estX=t(X)%*%V_est_inv%*%X
 beta=solve(XV_estX)%*%t(X)%*%V_est_inv%*%y
 res = y - X%*%beta
 P = V_est_inv - V_est_inv %*% X %*% solve(t(X) %*% V_est_inv %*% X) %*% t(X) %*% V_est_inv
 svdob<-eigen(P, symmetric=TRUE)
 Phalf<-svdob$vectors%*%diag(sqrt(round(svdob$values, 6)))%*%t(svdob$vectors)
 list("res"=res, "V_est"=V_est, "V_est_inv"=V_est_inv, "n"=n, "X"=X, "yid_order"=yid_order, "gid"=gid, "fa"=fa, "mo"=mo, "P"=P, "Phalf"=Phalf)
}

else if (method=="L-BFGS-B"){
reml.lik<-function(theta,y,X,K,I){
 #theta[1]: variance for family pheno1
 #theta[2]: variance for family pheno2
 #theta[3]: correlation for family
 #theta[4]: variance for family pheno1
 #theta[5]: variance for family pheno2
 #theta[6]: correlation for family
 K1<-theta[1]
 K2<-theta[2]
 K12<-theta[3]*sqrt(K1*K2)
 e1<-theta[4]
 e2<-theta[5]
 e12<-theta[6]*sqrt(e1*e2)
 K_est = matrix(c(K1,K12,K12,K2),nrow=2)
 e_est = matrix(c(e1,e12,e12,e2),nrow=2)
 V_K <- kronecker(K, K_est)
 V_e <- kronecker(I, e_est)
 V_est <- V_K + V_e
 V_est_inv <- solve(V_est)
 determinant1 = determinant(V_est, logarithm=T)[1:2]
 log_det_V_est <- determinant1$modulus[1] * determinant1$sign
 XV_estX=t(X)%*%V_est_inv%*%X
 determinant2 = determinant(XV_estX, logarithm=T)[1:2]
 log_det_XV_estX <- determinant2$modulus[1] * determinant2$sign
 beta=solve(XV_estX)%*%t(X)%*%V_est_inv%*%y
 mu_hat=X%*%beta
 logl<- log_det_V_est + log_det_XV_estX + t(y-mu_hat)%*%V_est_inv%*%(y-mu_hat)
 return(as.numeric(logl))
}

# Regular checks
n<-length(unique(yid))
m<-length(phenotype)
#if(any(unique(yid)!=gid)) warning("id in phenotypes (yid) and id in genotypes (gid) are not in the same order. This would not affect the results")
if(length(gid)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
if(!is.null(covariates)) {
  if(nrow(covariates)!=m) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

order <- seq(1, length(as.vector(yid)))
yid_order <- cbind.data.frame(yid, order)

type3 <- cbind(gid, fa, mo)
colnames(type3)[1] <- "yid"
type4 <- merge(yid_order, type3, by="yid")
type4 <- type4[order(type4$order),]

relation <- unique(cbind.data.frame(type4$yid, type4$fa, type4$mo))
relation <- as.matrix(relation)
K <- 2*kinship(relation[,1], relation[,2], relation[,3])
I <- diag(dim(relation)[1])

intercept <- rep(1,length(yid))
dummy <- model.matrix( ~ trait - 1)
y <- as.matrix(phenotype)

if(is.null(covariates)){
   X <- intercept
   X <- as.matrix(X)
   data <- data.frame(y=y, dummy)
   exprs<-paste("y ~ 1")}
else if(!is.null(covariates)){
   X <- cbind(intercept, as.matrix(covariates))
   X <- as.matrix(X)
   data <- data.frame(y=y, dummy, covariates)
   exprs<-paste("y ~ ", paste(names(covariates),collapse=" + "), sep="")}

 if(!eq.cov.effect){
 X1=X*rep(c(1,0), length(yid)/2)
 X2=X*rep(c(0,1), length(yid)/2)
 X=cbind(X1, X2)}

 data1<-data[data[,2]==1,]
 data2<-data[data[,3]==1,]

 if (Ninitial>1 & is.null(cor)){
  cor<-LL<-vector()
  all<-theta<-list()
  guess<-cor(resid(lm(as.formula(exprs), data=data1)), resid(lm(as.formula(exprs), data=data2)))
  if (sign(guess)==1)  cor<-runif(Ninitial, min=0, max=1)
  if (sign(guess)==-1)  cor<-runif(Ninitial, min=-1, max=0)
  cor[1]<-guess
  all[[1]] = optim(c(1, 1, cor[1], 1, 1, cor[1]), reml.lik, method="L-BFGS-B", lower=c(0,0,-1,0,0), upper=c(Inf,Inf,1,Inf,Inf),  y=y, X=X, K=K, I=I)
  theta[[1]] = all[[1]]$par
  LL[1] = -0.5*all[[1]]$value
  write.table(t(c(cor[1], LL[1])), file=LL.output, append=T, row.names=F, col.names=F, quote=F)
  for (i in 2:Ninitial){
   all[[i]] = optim(c(1, 1, cor[i], 1, 1, cor[i]), reml.lik, method="L-BFGS-B", lower=c(0,0,-1,0,0), upper=c(Inf,Inf,1,Inf,Inf), y=y, X=X, K=K, I=I)
   theta[[i]] = all[[i]]$par
   LL[i] = -0.5*all[[i]]$value
   write.table(t(c(cor[i], LL[i])), file=LL.output, append=T, row.names=F, col.names=F, quote=F)
  }
  theta = theta[[which(LL==max(LL))]]
  LL = max(LL)
 }

 else if (Ninitial==1 & is.null(cor)){
  cor<-cor(resid(lm(as.formula(exprs), data=data1)), resid(lm(as.formula(exprs), data=data2)))
  all = optim(c(1, 1, cor, 1, 1, cor), reml.lik, method="L-BFGS-B", lower=c(0,0,-1,0,0), upper=c(Inf,Inf,1,Inf,Inf), y=y, X=X, K=K, I=I)
  theta = all$par
  LL = -0.5*all$value
  write.table(t(c(cor, LL)), file=LL.output, row.names=F, col.names=F, quote=F)
 }
 else if (!is.null(cor)){
  all = optim(c(1, 1, cor, 1, 1, cor), reml.lik, method="L-BFGS-B", lower=c(0,0,-1,0,0), upper=c(Inf,Inf,1,Inf,Inf), y=y, X=X, K=K, I=I)
  theta = all$par
  LL = -0.5*all$value
  write.table(t(c(cor, LL)), file=LL.output, row.names=F, col.names=F, quote=F)
 }

 K1<-theta[1]
 K2<-theta[2]
 K12<-theta[3]*sqrt(K1*K2)
 e1<-theta[4]
 e2<-theta[5]
 e12<-theta[6]*sqrt(e1*e2)
 K_est = matrix(c(K1,K12,K12,K2),nrow=2)
 e_est = matrix(c(e1,e12,e12,e2),nrow=2)
 V_K <- kronecker(K, K_est)
 V_e <- kronecker(I, e_est)
 V_est <- V_K + V_e
 V_est_inv <- solve(V_est)
 XV_estX=t(X)%*%V_est_inv%*%X
 beta=solve(XV_estX)%*%t(X)%*%V_est_inv%*%y
 res = y - X%*%beta
 P = V_est_inv - V_est_inv %*% X %*% solve(t(X) %*% V_est_inv %*% X) %*% t(X) %*% V_est_inv
 svdob<-eigen(P, symmetric=TRUE)
 Phalf<-svdob$vectors%*%diag(sqrt(round(svdob$values, 6)))%*%t(svdob$vectors)
 list("res"=res, "V_est"=V_est, "V_est_inv"=V_est_inv, "n"=n, "X"=X, "yid_order"=yid_order, "gid"=gid, "fa"=fa, "mo"=mo, "P"=P, "Phalf"=Phalf)
}
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

