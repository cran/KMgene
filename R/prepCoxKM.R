#' KM for Survival Traits in GWAS Data (fit null model)
#'
#' This function (CoxKM) is used to perform KM analysis (Chen et al., 2014) for survival traits in GWAS data \cr
#'
#' @name prepCoxKM
#' @aliases prepCoxKM
#' @param Z A matrix of genotypes (class: matrix).  rows correspond to individuals and columns correspond to SNPs. Use ’NA’ for missing values. The column names of this matrix should correspond to SNP names in the SNP information file.
#' @param formula Base formula under null hypothesis (class: formula). For Cox models, the formula follows that of the coxph() function.
#' @param SNPInfo SNP Info file (class: data.frame). Must contain fields given in ’snpName’ and ’aggregateBy’.
#' @param snpNames The field of SNPInfo where the SNP identifiers are found. (Default="Name")
#' @param aggregateBy The field of SNPInfo on which the skat results were aggregated. (Default="gene")
#' @param data The data file in which to find variables in the formula (class: data.frame).
#' @param verbose Whether or not to print the progress bar (class: logical).
#' @importFrom seqMeta prepCox 
#' @return output: object as input for FbKM
#' @export
prepCoxKM <- function(Z, formula, SNPInfo=NULL, snpNames="Name", aggregateBy="gene", data=parent.frame(), verbose=FALSE) 
{
 prepCox(Z, formula, SNPInfo=SNPInfo, snpNames=snpNames, aggregateBy=aggregateBy, data=data, verbose=verbose) 
}

