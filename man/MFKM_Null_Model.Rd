% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MFKM_Null_Model.R
\name{MFKM_Null_Model}
\alias{MFKM_Null_Model}
\title{KM for Quantitative Traits in Multivariate Family GWAS Data (fit null model)}
\usage{
MFKM_Null_Model(
  phenotype,
  trait,
  yid,
  gid,
  fa,
  mo,
  covariates = NULL,
  Ninitial = 10,
  method = "Nelder-Mead",
  LL.output = NULL,
  cor = NULL,
  eq.cov.effect = F
)
}
\arguments{
\item{phenotype}{A vector of quantitative trait in the analysis (class: vector). The order should match the vector yid. Subjects with missing phenotypes are only used for kinship calculation.}

\item{trait}{A vector of multivariate traits (class: vector). The order should match the vector yid. No missing.}

\item{yid}{A vector of id (class: vector). Although it doesn't have to be sorted, observations from the same subject have to be connected with each other. The repeated id numbers indicate mutiple time points for one subject. Make sure it is not a factor. No missing.}

\item{gid}{A vector of id mapping to samples in genotype file (class: vector). So the order of samples in gid must be the same as the order in genotypes. Make sure it is not a factor. Although gid doesn't have to be in the same order as yid, it is suggested to make them sorted in the same order in order to make all files easily to be tracked. No missing.}

\item{fa}{A vector of father id (class: vector). The father id indicates the father of each subject. If this subject has no father in this data, the value is set to "NA". Make sure it is not factor.}

\item{mo}{A vector of mother id (class: vector). The mother id indicates the mother of each subject. If this subject has no mother in this data, the value is set to "NA". Make sure it is not factor.}

\item{covariates}{A matrix of covariates (class: data.frame). The order of rows should match the vector yid. Default NULL. Subjects with missing covariates are only used for kinship calculation.}

\item{Ninitial}{The number of times to try initial values. The default is 10 times. If Ninitial=1, the initial value "cor" is always equal to correlation(trait1|covariates, trait2|covariates). One should try multiple initial values in order to find max log-likelihood. This could be time consuming, depends on the sample size. The good thing is that null model only needs to be fitted once for the whole genome, so it's worth trying many initial values.}

\item{method}{The optimization method used in null model. The default method is an implementation of that of Nelder and Mead (1965), that uses only function values and is robust but relatively slow. Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each variable can be given a lower and/or upper bound.}

\item{LL.output}{Output all tried initial values and corresponding log-likelihoods. The initial value with max log-likelihood is used in the algorithm and it can be used for replication. The output file can be renamed.}

\item{cor}{Initial value. By default, it's not given, the program tries to find the best initial value. Once it's given, the program uses it as the only initial value. This is useful when one already knows the initial value corresponding to max log-likelihood.}

\item{eq.cov.effect}{Whether assume equal covariates effects on different traits (Default=False).}
}
\value{
output: object as input for MFKM
}
\description{
This function (MFKM) is used to perform KM analysis (Yan et al., 2015) for quantitative traits in GWAS multivariate family data. \cr
# It takes familial correlation as a kinship matrix
}
