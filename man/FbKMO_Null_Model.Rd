% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FbKMO_Null_Model.R
\name{FbKMO_Null_Model}
\alias{FbKMO_Null_Model}
\title{Optimal KM for Binary Traits in Familial GWAS Data (fit null model)}
\usage{
FbKMO_Null_Model(phenotype, id, fa, mo, family = "binomial", covariates = NULL)
}
\arguments{
\item{phenotype}{A vector of quantitative trait in the analysis (class: vector). The order should match the vector id. Subjects with missing phenotypes are only used for kinship calculation.}

\item{id}{A vector of id (class: vector). It can be either numeric or character. The id indicates each subject. Make sure it is not factor. No missing.}

\item{fa}{A vector of father id (class: vector). It can be either numeric or character. The father id indicates the father of each subject. If this subject has no father in this data, the value is set to "NA". Make sure it is not factor.}

\item{mo}{A vector of mother id (class: vector). It can be either numeric or character. The mother id indicates the mother of each subject. If this subject has no mother in this data, the value is set to "NA". Make sure it is not factor.}

\item{family}{Type of phenotype. (Default="binomial")}

\item{covariates}{A matrix of covariates (class: data.frame). The order of rows should match the vector id. Default NULL. Subjects with missing covariates are only used for kinship calculation.}
}
\value{
output: object as input for FbKMO
}
\description{
This function (FbKMO) is used to perform optimal KM analysis for binary traits in familial GWAS data \cr
# In the final correlation matrix, the covariance between parent and offspring is 0.5, the covariance between siblings is also 0.5
}
