#' KM for continuous Traits in Familial GWAS Data (calculate p-value)
#'
#' This function (FKM) is used to perform famSKAT analysis (Chen et al., 2013) for continuous traits in familial GWAS data
#'
#' @name FKM
#' @aliases FKM
#' @param obj results saved from FKM_Null_Model.R.
#' @param genotypes 1st column: gene name; 2nd column: snp name; 3rd-end columns: A matrix of genotypes for each subject (class: data.frame). The order of 3rd-end columns should match id. Coded as 0, 1, 2 and no missing.
#' This genotype file can be a big file containing all genes or it can be files containing one single gene.
#' @param gid A vector of id mapping to samples in genotype file (class: vector). So the order of samples in gid must be the same as the order in genotypes. Make sure it is not a factor. Although gid doesn't have to be in the same order as id, it is suggested to make them sorted in the same order in order to make all files easily to be tracked. No missing.
#' @param weights 1st column: gene name; 2nd column: snp name; 3rd column: A vector with the length equal to the number of variants in the test (class: data.frame). Default is Null indicating equal weight for all markers
#' @param acc Accuracy of numerical integration used in Davies' method. Default 1e-4.
#' @param append.write The name of pvalue output file. Write out p-values in real time. Don't need to wait until all genes are processed to get the final output.
#' @import CompQuadForm
#' @importFrom utils write.table
#' @return output: continuous trait family KM (F-KM) p-value
#' @examples
#' #####################################################################
#' ### Examples for Continuous Traits in Familial GWAS Data using KM ###
#' #####################################################################
#' ### Subject IDs are numeric ###
#' data("FKM_numID")
#' obj1 <- FKM_Null_Model(phenotype=fkm_n_y$y, id=fkm_n_y$id, fa=fkm_n_y$fa,
#' mo=fkm_n_y$mo, covariates=NULL)
#' pvalue1 <- FKM(obj=obj1, genotypes=fkm_n_gene, gid=fkm_n_geneid$gid, weights=NULL, 
#' append.write="./pvalues.out")
#' # Read in a list of genes files instead of a big file containing all genes
#' obj <- FKM_Null_Model(phenotype=fkm_n_y$y, id=fkm_n_y$id, fa=fkm_n_y$fa, mo=fkm_n_y$mo,
#' covariates=NULL)
#' gene <- split(fkm_n_gene, fkm_n_gene[,1])
#' for (k in 1:2) {
#'   gene[[k]]$gene <- as.character(gene[[k]]$gene)
#'   pvalue1 <- FKM(obj=obj, genotypes=gene[[k]], gid=fkm_n_geneid$gid, weights=NULL,
#'   append.write="./pvalues.out")
#' }
#' ### Subject IDs are character ###
#' data("FKM_charID")
#' obj1 <- FKM_Null_Model(phenotype=fkm_c_y$y, id=as.character(fkm_c_y$id),
#' fa=as.character(fkm_c_y$fa), mo=as.character(fkm_c_y$mo),  covariates=NULL)
#' pvalue1 <- FKM(obj=obj1, genotypes=fkm_c_gene, gid=as.character(fkm_c_geneid$gid),
#' weights=NULL)
#' @export
FKM <- function(obj, genotypes, gid, weights=NULL, acc=1e-4, append.write=NULL){
# Regular checks
n = obj$n
id_order = obj$id_order
if(!is.data.frame(genotypes)) stop("genotypes should be a data.frame!")
if(ncol(genotypes)!=n+2) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")

gid <- as.character(gid)
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

Q<-eig<-evals<-tmpout<-list()
p<-vector()
 P = obj$P
for (k in 1:length(genotypes2)){
 G = as.matrix(genotypes4[[k]][,-c(1:2)])
 class(G)<-"numeric"
 Q[[k]] = t(res) %*% V_inv %*% G %*% W[[k]] %*% t(G) %*% V_inv %*% res
 eig[[k]] = eigen(sqrt(W[[k]]) %*% t(G) %*% V_inv %*% P %*% V_inv %*% G %*% sqrt(W[[k]]), symmetric=T, only.values=T)
 evals[[k]] = eig[[k]]$values[eig[[k]]$values>1e-6*eig[[k]]$values[1]]

 tmpout[[k]]<-davies(Q[[k]], evals[[k]], acc=acc)
 p[k]<-tmpout[[k]]$Qq
 if(!is.null(append.write)){
 write.table(t(c(names(genotypes2)[k], signif(p[k], digits=6))), file=append.write, row.names=F, col.names=F, append=T, quote=F) }
}
cbind(names(genotypes2), signif(p, digits=6))
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.
# Chen H, Meigs JB, Dupuis J. 2013. Sequence kernel association test for quantitative traits in family samples. Genet Epidemiol. 37(2):196-204

