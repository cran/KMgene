#' KM for Binary Traits in Familial GWAS Data (fit null model)
#'
#' This function (FbKM) is used to perform KM analysis (Yan et al., 2015) for binary traits in familial GWAS data \cr
#' # In the final correlation matrix, the covariance between parent and offspring is 0.5, the covariance between siblings is also 0.5
#'
#' @name FbKM_Null_Model
#' @aliases FbKM_Null_Model
#' @param phenotype A vector of quantitative trait in the analysis (class: vector). The order should match the vector id. Subjects with missing phenotypes are only used for kinship calculation.
#' @param id A vector of id (class: vector). It can be either numeric or character. The id indicates each subject. Make sure it is not factor. No missing.
#' @param fa A vector of father id (class: vector). It can be either numeric or character. The father id indicates the father of each subject. If this subject has no father in this data, the value is set to "NA". Make sure it is not factor.
#' @param mo A vector of mother id (class: vector). It can be either numeric or character. The mother id indicates the mother of each subject. If this subject has no mother in this data, the value is set to "NA". Make sure it is not factor.
#' @param family Type of phenotype. (Default="binomial")
#' @param covariates A matrix of covariates (class: data.frame). The order of rows should match the vector id. Default NULL. Subjects with missing covariates are only used for kinship calculation.
#' @import MASS
#' @importFrom mgcv extract.lme.cov
#' @import kinship2
#' @importFrom stats as.formula
#' @include glmmPQL.R
#' @return output: object as input for FbKM
#' @export
FbKM_Null_Model <- function(phenotype, id, fa, mo, family="binomial", covariates=NULL){

# Regular checks
n<-length(phenotype)
if(length(id)!=n) stop("Number of individuals inconsistent between phenotype and id. Check your data...")
if(!is.null(covariates)) {
  if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

order <- seq(1, length(as.vector(id)))
id_order <- cbind.data.frame(id, order)

y <- as.matrix(phenotype)
K <- kinship(id, fa, mo)

  id <- order

intercept <- rep(1,length(id))
if(is.null(covariates)){
   X <- as.matrix(intercept)
   exprs<-paste("y ~ 1")}
else if(!is.null(covariates)){
   X <- cbind(intercept, as.matrix(covariates))
   X <- as.matrix(X)
   exprs<-paste("y ~ ", paste(names(covariates), collapse=" + "))}

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
   id=id[-index]
   id_order=id_order[-index,]
 }

 cs.K <- corSymm(2*K[lower.tri(K)],fixed=T)
 # id <- as.factor(id)
 id <- as.matrix(id)
 cs.K <- Initialize(cs.K,data=id)

 model <- glmmPQL2(as.formula(exprs), random=~1|id, correlation=cs.K, data=data, family=family, control=lmeControl(opt="optim"))

 beta <- model$fit$coefficients$fixed
 V <- extract.lme.cov(model$fit, data)
 V_inv <- solve(V)

 # test beta value
 #beta2 <- solve(t(X)%*%V_inv%*%X)%*%t(X)%*%V_inv%*%as.matrix(model$y_star)

 res = as.matrix(model$y_star - X%*%beta)
 P = V - X %*% solve(t(X) %*% V_inv %*% X) %*% t(X)
 list("res"=res, "V"=V, "V_inv"=V_inv, "id_order"=id_order, "X"=X, "n"=n, "P"=P)
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

