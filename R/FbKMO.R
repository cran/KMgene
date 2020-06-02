#' Optimal KM for Non-continuous Traits in Familial GWAS Data (calculate p-value)
#'
#' This function (FbKMO) is used to perform the optimal KM analysis for binary traits in familial GWAS data
#'
#' @name FbKMO
#' @aliases FbKMO
#' @param obj results saved from FbKMO_Null_Model.
#' @param genotypes 1st column: gene name; 2nd column: snp name; 3rd-end columns: A matrix of genotypes for each subject (class: data.frame). The order of 3rd-end columns should match id. Coded as 0, 1, 2 and no missing. This genotype file can be a big file containing all genes or it can be files containing one single gene.
#' @param gid A vector of id mapping to samples in genotype file (class: vector). So the order of samples in gid must be the same as the order in genotypes. Make sure it is not a factor. Although gid doesn't have to be in the same order as id, it is suggested to make them sorted in the same order in order to make all files easily to be tracked. No missing.
#' @param weights 1st column: gene name; 2nd column: snp name; 3rd column: A vector with the length equal to the number of variants in the test (class: data.frame). Default is Null indicating equal weight for all markers
#' @param r.all A list of predefined proportion of linear kernel and burden test. When r.all=0, regular kernel machine test (FbKM); when r.all=1, burden test.
#' @param acc Accuracy of numerical integration used in Davies' method for individual r.all p-values. Default 1e-4.
#' @param acc2 Accuracy of numerical integration used in Davies' method for the final p-value. Default 1e-4.
#' @param append.write The name of pvalue output file. Write out p-values in real time. Don't need to wait until all genes are processed to get the final output.
#' @import CompQuadForm
#' @importFrom utils write.table
#' @importFrom stats qchisq
#' @importFrom stats integrate
#' @include Get_Liu_Params_Mod.R
#' @include SKAT_Optimal_Integrate_Func_Davies.R
#' @return output: binary trait optimal family KM (Fb-KMO) p-value
#' @examples
#' #########################################################################
#' ### Examples for Binary Traits in Familial GWAS Data using optimal KM ###
#' #########################################################################
#' ### Subject IDs are numeric ###
#' data("FbKM_numID")
#' obj1 <- FbKMO_Null_Model(phenotype=fbkm_n_y$y, id=fbkm_n_y$id, fa=fbkm_n_y$fa,
#' mo=fbkm_n_y$mo, family="binomial", covariates=NULL)
#' pvalue1 <- FbKMO(obj=obj1, genotypes=fbkm_n_gene, gid=fbkm_n_geneid$gid, weights=NULL) 
#' # Read in a list of genes files instead of a big file containing all genes
#' obj <- FbKMO_Null_Model(phenotype=fbkm_n_y$y, id=fbkm_n_y$id, fa=fbkm_n_y$fa,
#' mo=fbkm_n_y$mo, family="binomial", covariates=NULL)
#' gene <- split(fbkm_n_gene, fbkm_n_gene[,1])
#' for (k in 1:2) {
#'   gene[[k]]$gene <- as.character(gene[[k]]$gene)
#'   pvalue1 <- FbKMO(obj=obj, genotypes=gene[[k]], gid=fbkm_n_geneid$gid, weights=NULL) 
#' }
#' ### Subject IDs are character ###
#' data("FbKM_charID")
#' obj1 <- FbKMO_Null_Model(phenotype=fbkm_c_y$y, id=as.character(fbkm_c_y$id),
#' fa=as.character(fbkm_c_y$fa), mo=as.character(fbkm_c_y$mo), family="binomial",
#' covariates=NULL)
#' pvalue1 <- FbKMO(obj=obj1, genotypes=fbkm_c_gene, gid=as.character(fbkm_c_geneid$gid),
#' weights=NULL)
#' @export
FbKMO <- function(obj, genotypes, gid, weights=NULL, acc=1e-4, acc2=1e-4, r.all=c(0, 0.25, 0.5 ,0.75, 1), append.write=NULL){
# Regular checks
n = obj$n
id_order = obj$id_order
if(!is.data.frame(genotypes)) stop("genotypes should be a data.frame!")
if(ncol(genotypes)!=n+2) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")

genotypes2 <- split(genotypes, genotypes[,1])
genotypes3 <- genotypes4 <- list()
for (k in 1:length(genotypes2)){
 genotypes3[[k]] <- cbind(gid, t(genotypes2[[k]][,-c(1:2)]))
 colnames(genotypes3[[k]])[1] <- "id"
 genotypes4[[k]] <- merge(id_order, genotypes3[[k]], by="id")
 genotypes4[[k]] <- genotypes4[[k]][order(genotypes4[[k]]$order),]
}

# Weight according to beta density function
if(is.null(weights)){
W <- list()
for (k in 1:length(genotypes2)){
# w <- dbeta(colMeans(t(genotypes2[[k]][,-c(1:2)]))/2, 1, 25)
# if(length(w)==1) W[[k]] <- w^2
# else W[[k]] <- diag(w^2)
 if (dim(as.matrix(t(genotypes2[[k]][,-c(1:2)])))[2] == 1) {W[[k]] <-1} else
 {W[[k]] <- diag(1, dim(as.matrix(t(genotypes2[[k]][,-c(1:2)])))[2])}
}}
else if(!is.null(weights)){
weights2 <- split(weights, weights[,1])
W <- list()
for (k in 1:length(weights2)){
 W[[k]] <- diag(weights2[[k]][,-c(1:2)])
}}

V = obj$V
V_inv = obj$V_inv
res = obj$res
X = obj$X

Q<-eig<-evals<-pskat<-list()
pvalue<-vector()
 P = obj$P
 Phalf = obj$Phalf

 IDX <- which(r.all >= 0.999)
 if (length(IDX) > 0) r.all[IDX] <- 0.999

for (k in 1:length(genotypes2)){
 G = as.matrix(genotypes4[[k]][,-c(1:2)])
 class(G)<-"numeric"
 one = matrix(1, nrow=dim(G)[2], ncol=dim(G)[2])
 for (i in 1:length(r.all)){
   R = (1-r.all[i])*W[[k]] + r.all[i]*sqrt(W[[k]])%*%one%*%sqrt(W[[k]])
   Q[[i]] = t(res) %*% V_inv %*% G %*% R %*% t(G) %*% V_inv %*% res
   svdob<-eigen(R, symmetric=TRUE)
   Rhalf<-svdob$vectors%*%diag(sqrt(round(svdob$values, 6)))%*%t(svdob$vectors)
   eig[[i]] = eigen(Rhalf %*% t(G) %*%  P %*% G %*% Rhalf, symmetric=T, only.values=T)
   evals[[i]] = eig[[i]]$values[eig[[i]]$values>1e-6*eig[[i]]$values[1]]

   tmpout<-davies(Q[[i]], evals[[i]], acc=acc)
   pskat[[i]]<-tmpout$Qq
 }

 pmin = min(unlist(pskat))
 if(pmin==2) pmin=1
 if(pmin<0){
  pmin=0
  warning("negative pvalue detected. change/decrease acc")
 }
 # bKM_Optimal_Param
    Z <- G%*%sqrt(W[[k]])
    Z1 <- Phalf%*%Z
    z_mean <- rowMeans(Z1)
    p.m <- dim(Z1)[2]
    Z_mean <- matrix(rep(z_mean, p.m), ncol = p.m, byrow = FALSE)
    cof1 <- (t(z_mean) %*% Z1)[1, ]/sum(z_mean^2)
    Z.item1 <- Z_mean %*% diag(cof1)
    Z.item2 <- Z1 - Z.item1
    W3.2.t <- t(Z.item2) %*% Z.item2
    eig_W3.2.t = eigen(W3.2.t, symmetric=T, only.values=T)
    lambda <- eig_W3.2.t$values[eig_W3.2.t$values>1e-6*eig_W3.2.t$values[1]]
    W3.3.item <- sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*%  Z.item2)) * 4
    MuQ <- sum(lambda)
    VarQ <- sum(lambda^2) * 2 + W3.3.item
    KerQ <- sum(lambda^4)/(sum(lambda^2))^2 * 12
    Df <- 12/KerQ
    tau <- rep(0, length(r.all))
    for (i in 1:length(r.all)) {
        r.corr <- r.all[i]
        term1 <- p.m^2 * r.corr + sum(cof1^2) * (1 - r.corr)
        tau[i] <- sum(term1) * sum(z_mean^2)
    }
    param.m <- list(MuQ = MuQ, VarQ = VarQ, KerQ = KerQ, lambda = lambda, VarRemain = W3.3.item, Df = Df, tau = tau)

 # Calculate the percentile (pmin.q)
 # bKM_Optimal_Each_Q
    c1 <- pmin.q <- vector()
    param.mat <- NULL
    for (i in 1:length(r.all)) {
        lambda.temp <- evals[[i]]
        c1[1] <- sum(lambda.temp)
        c1[2] <- sum(lambda.temp^2)
        c1[3] <- sum(lambda.temp^3)
        c1[4] <- sum(lambda.temp^4)
        param.temp <- Get_Liu_Params_Mod(c1)
        muQ <- param.temp$muQ
        varQ <- param.temp$sigmaQ^2
        df <- param.temp$l
        param.mat <- rbind(param.mat, c(muQ, varQ, df))
    }
    for (i in 1:length(r.all)) {
        muQ <- param.mat[i, 1]
        varQ <- param.mat[i, 2]
        df <- param.mat[i, 3]
        q.org <- qchisq(1 - pmin, df = df)
        q.q <- (q.org - df)/sqrt(2 * df) * sqrt(varQ) + muQ
        pmin.q[i] <- q.q
    }

 # bKM_Optimal_PValue_Davies
    re <- try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=40, subdivisions=1000, pmin.q=pmin.q, param.m=param.m, r.all=r.all, abs.tol=10^-25, acc=acc2), silent=TRUE)
    pvalue[k] <- 1 - re[[1]]
    if (!is.null(pmin)) {
        if (pmin * length(r.all) < pvalue[k]) {
            pvalue[k] = pmin * length(r.all)
        }
    }

 if(!is.null(append.write)){
  write.table(t(c(names(genotypes2)[k], signif(pvalue[k], digits=6))), file=append.write, row.names=F, col.names=F, append=T, quote=F) }
}
cbind(names(genotypes2), signif(pvalue, digits=6))
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

