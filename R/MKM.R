#' KM for Quantitative Traits in Multivariate GWAS Data (calculate p-value)
#'
#' This function (MKM) is used to perform KM analysis (Tzeng et al. 2012) for quantitative traits in GWAS multivariate data.
#'
#' @name MKM
#' @aliases MKM
#' @param obj results saved from MKM_Null_Model.
#' @param genotypes 1st column: gene name; 2nd column: snp name; 3rd-end columns: A matrix of genotypes for each subject (class: data.frame). The order of 3rd-end columns should match unique(yid). Coded as 0, 1, 2 and no missing. This genotype file can be a big file containing all genes or it can be files containing one single gene.
#' @param weights 1st column: gene name; 2nd column: snp name; 3rd column: A vector with the length equal to the number of variants in the test (class: data.frame). Default is Null indicating equal weight for all markers
#' @param gid A vector of id mapping to samples in genotype file (class: vector). So the order of samples in gid must be the same as the order in genotypes. Make sure it is not a factor. Although gid doesn't have to be in the same order as yid, it is suggested to make them sorted in the same order in order to make all files easily to be tracked. No missing.
#' @param acc Accuracy of numerical integration used in Davies' method. Default 1e-4.
#' @param append.write The name of pvalue output file. Write out p-values in real time. Don't need to wait until all genes are processed to get the final output.
#' @param eq.gen.effect Whether assume equal genetic effects on different traits (default = False).
#' @import CompQuadForm
#' @importFrom utils write.table
#' @return output: multivariate KM (M-KM) p-value
#' @examples
#' ###############################################################################
#' ### Examples for Multivariate (two) Continuous Traits in GWAS Data using KM ###
#' ###############################################################################
#' ### Subject IDs are numeric ###
#' data("MKM_numID")
#' obj1 <- MKM_Null_Model(phenotype=mkm_n_ph$y, trait=mkm_n_ph$trait, yid=mkm_n_ph$id,
#' covariates=NULL)
#' pvalue1 <- MKM(obj=obj1, genotypes=mkm_n_gene, gid=mkm_n_geneid$gid, weights=NULL) 
#' # Read in a list of genes files instead of a big file containing all genes
#' obj <- MKM_Null_Model(phenotype=mkm_n_ph$y, trait=mkm_n_ph$trait, yid=mkm_n_ph$id,
#' covariates=NULL)
#' gene <- split(mkm_n_gene, mkm_n_gene[,1])
#' for (k in 1:2) {
#'   gene[[k]]$gene <- as.character(gene[[k]]$gene)
#'   pvalue1 <- MKM(obj=obj, genotypes=gene[[k]], gid=mkm_n_geneid$gid, weights=NULL)
#' }
#' ### Subject IDs are character ###
#' data("MKM_charID")
#' obj1 <- MKM_Null_Model(phenotype=mkm_c_ph$y, trait=mkm_c_ph$trait,
#' yid=as.character(mkm_c_ph$id), covariates=NULL)
#' pvalue1 <- MKM(obj=obj1, genotypes=mkm_c_gene, gid=as.character(mkm_c_geneid$gid),
#' weights=NULL)
#' @export
MKM <- function(obj, genotypes, gid, weights=NULL, acc=1e-4, append.write=NULL, eq.gen.effect=F){

yid = obj$yid
n = obj$n
yid_order = obj$yid_order
dummy = obj$dummy

# Regular checks
#if(any(unique(yid)!=gid)) warning("id in phenotypes (yid) and id in genotypes (gid) are not in the same order. This would not affect the results")
if(length(gid)!=n) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
if(!is.data.frame(genotypes)) stop("genotypes should be a data.frame!")
if(ncol(genotypes)!=n+2) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")

genotypes2 <- split(genotypes, genotypes[,1])
genotypes3 <- genotypes4 <- list()
for (k in 1:length(genotypes2)){
 genotypes3[[k]] <- cbind(gid, t(genotypes2[[k]][,-c(1:2)]))
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
dummy = obj$dummy

Q<-eig<-evals<-tmpout<-list()
p<-vector()
 P = obj$P
for (k in 1:length(genotypes2)){
 G = as.matrix(genotypes4[[k]][,-c(1:2)])
 class(G)<-"numeric"
 if(!eq.gen.effect){
   G1=G*dummy[,1]
   G2=G*dummy[,2]
   G=cbind(G1,G2)}

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

