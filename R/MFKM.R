#' KM for Quantitative Traits in Multivariate Family GWAS Data (calculate p-value)
#'
#' This function (MFKM) is used to perform KM analysis (Yan et al., 2015) for quantitative traits in GWAS multivariate family data.
#'
#' @name MFKM
#' @aliases MFKM
#' @param obj results saved from MFKM_Null_Model.
#' @param genotypes 1st column: gene name; 2nd column: snp name; 3rd-end columns: A matrix of genotypes for each subject (class: data.frame). The order of 3rd-end columns should match unique(yid). Coded as 0, 1, 2 and no missing. This genotype file can be a big file containing all genes or it can be files containing one single gene.
#' @param weights 1st column: gene name; 2nd column: snp name; 3rd column: A vector with the length equal to the number of variants in the test (class: data.frame). Default is Null indicating equal weight for all markers.
#' @param acc Accuracy of numerical integration used in Davies' method. Default 1e-4.
#' @param append.write The name of pvalue output file. Write out p-values in real time. Don't need to wait until all genes are processed to get the final output.
#' @param eq.gen.effect Whether assume equal genetic effects on different traits (default = False).
#' @import CompQuadForm
#' @importFrom utils write.table
#' @return output: multivariate family KM (MF-KM) p-value
#' @examples
#' ########################################################################################
#' ### Examples for Multivariate (two) Continuous Traits in Familial GWAS Data using KM ###
#' ########################################################################################
#' ### Subject IDs are numeric ###
#' data("MFKM_numID")
#' #If Ninitial=1, the initial value "cor" is always equal to
#' #correlation(trait1|covariates, trait2|covariates).
#' obj1 <- MFKM_Null_Model(phenotype=mfkm_n_y$y, trait=mfkm_n_y$trait, yid=mfkm_n_y$id,
#' gid=mfkm_n_geneid$gid, fa=mfkm_n_geneid$fa, mo=mfkm_n_geneid$mo, covariates=NULL,
#' Ninitial=1)
#' #One should try multiple initial values in order to find max log-likelihood. The default
#' #is 10 times.
#' #This could be time consuming, depends on the sample size. The good thing is that null
#' #model only needs to be fitted once for the whole genome, so it's worth trying many
#' #initial values. The initial value with max logl can be saved in ./LogLikelihood.txt
#' #for reuse.
#' pvalue1 <- MFKM(obj=obj1, genotypes=mfkm_n_gene, weights=NULL)
#' #If one wants to replicate the results, finds the initial value "cor" with the max
#' #loglikelihood in ./LogLikelihood.txt that is output by default. The 1st column is
#' #"cor"; 2nd column is "log-likelihood".
#' obj1 <- MFKM_Null_Model(phenotype=mfkm_n_y$y, trait=mfkm_n_y$trait, yid=mfkm_n_y$id,
#' gid=mfkm_n_geneid$gid, fa=mfkm_n_geneid$fa, mo=mfkm_n_geneid$mo, covariates=NULL,
#' cor=0.687439771651474)
#' #This "cor" is calcualted based on 3 initial values, "Ninitial=3".
#' pvalue1 <- MFKM(obj=obj1, genotypes=mfkm_n_gene, weights=NULL)
#' #Introduce missing into covariates and outcome.
#' #When samples have missing values in outcome or covariates, those samples are used only for
#' #kinship calculation.
#' mfkm_n_covariates[1,] = NA
#' mfkm_n_y$y[4] = NA
#' obj1 <- MFKM_Null_Model(phenotype=mfkm_n_y$y, trait=mfkm_n_y$trait, yid=mfkm_n_y$id,
#' gid=mfkm_n_geneid$gid, fa=mfkm_n_geneid$fa, mo=mfkm_n_geneid$mo, covariates=mfkm_n_covariates,
#' Ninitial=1)
#' pvalue1 <- MFKM(obj=obj1, genotypes=mfkm_n_gene, weights=NULL)
#' #Read in a list of genes files instead of a big file containing all genes
#' obj <- MFKM_Null_Model(phenotype=mfkm_n_y$y, trait=mfkm_n_y$trait, yid=mfkm_n_y$id,
#' gid=mfkm_n_geneid$gid, fa=mfkm_n_geneid$fa, mo=mfkm_n_geneid$mo, covariates=NULL,
#' Ninitial=1)
#' gene <- split(mfkm_n_gene, mfkm_n_gene[,1])
#' for (k in 1:2) {
#'   gene[[k]]$gene <- as.character(gene[[k]]$gene)
#'   pvalue1 <- MFKM(obj=obj, genotypes=gene[[k]], weights=NULL, append.write=
#'   "./pvalues.out")
#' }
#' ### Subject IDs are character ###
#' data("MFKM_charID")
#' obj1 <- MFKM_Null_Model(phenotype=mfkm_c_y$y, trait=mfkm_c_y$trait,
#' yid=as.character(mfkm_c_y$id),
#' gid=as.character(mfkm_c_geneid$gid), fa=as.character(mfkm_c_geneid$fa),
#' mo=as.character(mfkm_c_geneid$mo), covariates=NULL, Ninitial=1)
#' pvalue1 <- MFKM(obj=obj1, genotypes=mfkm_c_gene, weights=NULL)
#' @export
MFKM <- function(obj, genotypes, weights=NULL, acc=1e-4, append.write=NULL, eq.gen.effect=F){

n = obj$n
yid_order = obj$yid_order
gid = obj$gid
fa = obj$fa
mo = obj$mo

# Regular checks
if(!is.data.frame(genotypes)) stop("genotypes should be a data.frame!")
if(ncol(genotypes)!=n+2) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")

genotypes2 <- split(genotypes, genotypes[,1])
genotypes3 <- genotypes4 <- list()
for (k in 1:length(genotypes2)){
 genotypes3[[k]] <- cbind(gid, fa, mo, t(genotypes2[[k]][,-c(1:2)]))
 colnames(genotypes3[[k]])[1] <- "yid"
 genotypes4[[k]] <- merge(yid_order, genotypes3[[k]], by="yid")
 genotypes4[[k]] <- genotypes4[[k]][order(genotypes4[[k]]$order),]
}

# Weight according to beta density function
if(is.null(weights)){
W <- list()
for (k in 1:length(genotypes2)){
# w <- dbeta(colMeans(t(genotypes2[[k]][,-c(1:2)]))/2, 1, 25)
# if(eq.gen.effect) W[[k]] <- diag(w^2)
# else if(!eq.gen.effect) W[[k]] <- diag(rep(w^2, 2))}}
 w <- rep(1, dim(as.matrix(t(genotypes2[[k]][,-c(1:2)])))[2])
 if(eq.gen.effect) W[[k]] <- diag(w)
 else if(!eq.gen.effect) W[[k]] <- diag(rep(w, 2))}}
else if(!is.null(weights)){
weights2 <- split(weights, weights[,1])
W <- list()
for (k in 1:length(weights2)){
 if(eq.gen.effect) W[[k]] <- diag(weights2[[k]][,-c(1:2)])
 else if(!eq.gen.effect) W[[k]] <- diag(rep(weights2[[k]][,-c(1:2)], 2))
}}

res = obj$res
V_est = obj$V_est
V_est_inv = obj$V_est_inv
X = obj$X

Q<-eig<-evals<-tmpout<-list()
p<-vector()
 P = obj$P
for (k in 1:length(genotypes2)){
 G = as.matrix(genotypes4[[k]][,-c(1:4)])
 class(G)<-"numeric"
 if(!eq.gen.effect){
 G1=G*rep(c(1,0), length(yid_order)/2)
 G2=G*rep(c(0,1), length(yid_order)/2)
 G=cbind(G1, G2)}

 class(G)<-"numeric"
 Q[[k]] = t(res) %*% V_est_inv %*% G %*% W[[k]] %*% t(G) %*% V_est_inv %*% res
 eig[[k]] = eigen(sqrt(W[[k]]) %*% t(G) %*% V_est_inv %*% P %*% V_est_inv %*% G %*% sqrt(W[[k]]), symmetric=T, only.values=T)
 evals[[k]] = eig[[k]]$values[eig[[k]]$values>1e-6*eig[[k]]$values[1]]

 tmpout[[k]]<-davies(Q[[k]], evals[[k]], acc=acc)
 p[k] <- tmpout[[k]]$Qq
 if(!is.null(append.write)){
  write.table(t(c(names(genotypes2)[k], signif(p[k], digits=6))), file=append.write, row.names=F, col.names=F, append=T, quote=F) }
}
cbind(names(genotypes2), signif(p, digits=6))
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

