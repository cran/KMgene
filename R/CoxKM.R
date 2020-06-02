#' KM for Survival Traits in GWAS Data (calculate p-value)
#'
#' This function (CoxKM) is used to perform the KM analysis (Chen et al., 2014) for survival traits in GWAS data
#'
#' @name CoxKM
#' @aliases CoxKM
#' @param object Object saved from prepCoxKM.
#' @param SNPInfo The SNP Info file. This should contain the fields listed in snpNames and aggregateBy. Only SNPs in this table will be meta analyzed, so this may be used to restrict the analysis.
#' @param wts Either a function to calculate testing weights, or a character specifying a vector of weights in the SNPInfo file. For skatMeta the default are the ‘beta’ weights.
#' @param method The p-value calculation method. Default is "saddlepoint", "integration" is the Davies method used in the SKAT package.
#' @param snpNames The field of SNPInfo where the SNP identifiers are found. Default is "Name"
#' @param aggregateBy The field of SNPInfo on which the skat results were aggregated. Default is "gene".
#' @param mafRange Range of MAF’s to include in the analysis (endpoints included). Default is all SNPs (0 <= MAF <= 0.5).
#' @param verbose Logical. Whether progress bars should be printed.
#' @importFrom seqMeta skatMeta 
#' @importFrom stats dbeta
#' @return output: survival trait KM (CoxKM) p-value
#' @examples
#' ##########################################################
#' ### Examples for Survival Traits in GWAS Data using KM ###
#' ##########################################################
#' data("CoxKM_data")
#' library(seqMeta)
#' cov <- as.matrix(pheno1[,2:3])
#' p<-vector()
#' for (i in 1:2) {
#'   Gene <- unique(SNPInfo$gene)[i]
#'   choose <- SNPInfo$gene == Gene
#'   geno <- Z1[,choose]
#'   cohort1 <- prepCoxKM(Z=Z1, Surv(time, status)~strata(sex)+bmi, 
#'   SNPInfo=SNPInfo[choose,], data=pheno1)
#'   p[i] <- CoxKM(cohort1, SNPInfo=SNPInfo[choose,])
#' }
#' @export

CoxKM <- function(object, SNPInfo=NULL, wts=function(maf){dbeta(maf, 1, 25)}, 
         method="saddlepoint", snpNames="Name", aggregateBy="gene", mafRange=c(0, 0.5), verbose=FALSE) 
{
 skatMeta(object, SNPInfo=SNPInfo, wts=wts, method=method, snpNames=snpNames, aggregateBy=aggregateBy, mafRange=mafRange, verbose=verbose)
}
