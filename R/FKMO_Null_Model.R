#' Optimal KM for continuous Traits in Familial GWAS Data (fit null model)
#'
#' This function (FKMO) is used to perform optimal famSKAT analysis for continuous traits in familial GWAS data \cr
#' # In the final correlation matrix, the covariance between parent and offspring is 0.5, the covariance between siblings is also 0.5
#'
#' @name FKMO_Null_Model
#' @aliases FKMO_Null_Model
#' @param phenotype A vector of quantitative trait in the analysis (class: vector). The order should match the vector id. Subjects with missing phenotypes are only used for kinship calculation.
#' @param id A vector of id (class: vector). It can be either numeric or character. The id indicates each subject. Make sure it is not factor. No missing.
#' @param fa A vector of father id (class: vector). It can be either numeric or character. The father id indicates the father of each subject. If this subject has no father in this data, the value is set to "NA". Make sure it is not factor.
#' @param mo A vector of mother id (class: vector). It can be either numeric or character. The mother id indicates the mother of each subject. If this subject has no mother in this data, the value is set to "NA". Make sure it is not factor.
#' @param covariates A matrix of covariates (class: data.frame). The order of rows should match the vector id. Default NULL. Subjects with missing covariates are only used for kinship calculation.
#' @import kinship2
#' @importFrom coxme lmekin
#' @importFrom stats as.formula
#' @return output: object as input for FKMO
#' @export
FKMO_Null_Model <- function(phenotype, id, fa, mo, covariates=NULL){

# Regular checks
n<-length(phenotype)
if(length(id)!=n) stop("Number of individuals inconsistent between phenotype and id. Check your data...")
if(!is.null(covariates)) {
  if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

order <- seq(1, length(as.vector(id)))
id_order <- cbind.data.frame(id, order)

y <- as.matrix(phenotype)
id <- as.character(id)
fa <- as.character(fa)
mo <- as.character(mo)
K <- kinship(id, fa, mo)

intercept <- rep(1,length(id))
if(is.null(covariates)){
   X <- as.matrix(intercept)
   exprs<-paste("y ~ 1 + (1|id)")}
else if(!is.null(covariates)){
   X <- cbind(intercept, as.matrix(covariates))
   X <- as.matrix(X)
   exprs<-paste("y ~ ", paste(names(covariates), collapse=" + "), "+ (1|id)")}

 data_pre <- data.frame(id=id,y=y)

 if(is.null(covariates)) data <- data_pre
 else if(!is.null(covariates)) data <- cbind(data_pre, covariates)

 # subjects with missing phenotypes are only used for kinship calculation, then they are deleted
 index = unique(which(is.na(data), arr.ind=T)[,1])
 if (!identical(index, integer(0))) {
   data=data[-index,]
   K=K[-index,-index]
   X=X[-index,]
   y=y[-index,]
   id_order=id_order[-index,]
 }

 mix_model <- lmekin(as.formula(exprs), data=data, varlist=list(K), method="REML")
 beta = fixef(mix_model)
 V <- mix_model$vcoef$id*K+mix_model$sigma^2*diag(dim(data)[1])
 V_inv <- solve(V)

 # The following beta calculation methods give same results
 #beta = solve(t(X)%*%V_est_inv%*%X)%*%t(X)%*%V_est_inv%*%y
 res = y - X%*%beta
 P = V_inv - V_inv %*% X %*% solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
 svdob<-eigen(P, symmetric=TRUE)
 Phalf<-svdob$vectors%*%diag(sqrt(round(svdob$values, 6)))%*%t(svdob$vectors)
 list("res"=res, "V"=V, "V_inv"=V_inv, "id_order"=id_order, "X"=X, "n"=n, "P"=P, "Phalf"=Phalf)
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.
# Chen H, Meigs JB, Dupuis J. 2013. Sequence kernel association test for quantitative traits in family samples. Genet Epidemiol. 37(2):196-204

