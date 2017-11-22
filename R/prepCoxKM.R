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
#' @importFrom survival coxph
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame
#' @importFrom stats coef
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom survival coxph.control
#' @importFrom stats na.omit
#' @importFrom stats var
#' @importFrom stats cov
#' @importFrom Matrix Matrix
#' @useDynLib KMgene
#' @include coxlr.fit.R
#' @include ginv_s.R
#' @return output: object as input for FbKM
#' @exportMethod c.KMgene c.KMgene
#' @export
prepCoxKM <- function (Z, formula, SNPInfo = NULL, snpNames = "Name", aggregateBy = "gene", 
    data = parent.frame(), verbose = FALSE) 
{
    env <- environment()
    if (is.null(SNPInfo)) {
        warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
    }
    nullmodel <- coxph(formula = formula, data = data) ### use survival package
    nullmodel$strata <- eval(parse(text = rownames(attr(nullmodel$terms, 
        "factors"))[attr(nullmodel$terms, "specials")$strata]), 
        envir = data)
    X <- model.matrix(nullmodel, data)
    rn <- row.names(model.frame(nullmodel, data = data))
    nullcoef <- coef(nullmodel)
    mysnps <- colnames(Z)
    SNPInfo[, aggregateBy] <- as.character(SNPInfo[, aggregateBy])
    which.snps.Z <- colnames(Z) %in% SNPInfo[, snpNames]
    ZtoSI <- match(SNPInfo[, snpNames], mysnps[which.snps.Z])
    nsnps <- sum(!is.na(ZtoSI))
    if (nsnps == 0) {
        stop("no column names in Z match SNP names in the SNP Info file!")
    }
    n = nrow(Z)
    if (verbose) {
        cat("\n Calculating signed LRTs... Progress:\n")
        pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
        pb.i <- 0
    }
    maf0 <- colMeans(Z, na.rm = TRUE)[which.snps.Z]/2
    maf0[is.nan(maf0)] <- -1
    maf <- maf0[ZtoSI]
    names(maf) <- SNPInfo[, snpNames]
    zlrt <- apply(Z[, which.snps.Z, drop = FALSE], 2, function(z) {
        if (any(is.na(z))) {
            if (all(is.na(z))) 
                z <- rep(0, length(z))
            mz <- mean(z, na.rm = TRUE)
            z[is.na(z)] <- mz
        }
        if (verbose) {
            assign("pb.i", get("pb.i", env) + 1, env)
            if (get("pb.i", env)%%ceiling(nsnps/100) == 0) 
                setTxtProgressBar(get("pb", env), get("pb.i", 
                  env))
        }
        model <- coxlr.fit(cbind(z, X), nullmodel$y, nullmodel$strata, 
            NULL, init = c(0, nullcoef), coxph.control(iter.max = 100), 
            NULL, "efron", rn) 
        return(sign(coef(model)[1]) * sqrt(2 * diff(model$loglik)))
    })[ZtoSI]
    names(zlrt) <- SNPInfo[, snpNames]
    zlrt[is.na(zlrt)] <- 0
    if (verbose) 
        close(pb)
    zlrt[maf == 0] <- 0
    maf[!(SNPInfo[, snpNames] %in% colnames(Z))] <- -1
    zlrt <- split(zlrt, SNPInfo[, aggregateBy])
    maf <- split(maf, SNPInfo[, aggregateBy])
    ngenes <- length(unique(SNPInfo[, aggregateBy]))
    if (verbose) {
        cat("\n Calculating covariance... Progress:\n")
        pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
        pb.i <- 0
    }
    re <- as.list(by(SNPInfo[, snpNames], SNPInfo[, aggregateBy], 
        function(snp.names) {
            inds <- match(snp.names, colnames(Z))
            mcov <- matrix(0, length(snp.names), length(snp.names))
            if (length(na.omit(inds)) > 0) {
                Z0 <- as.matrix(Z[, na.omit(inds), drop = FALSE])
                if (any(is.na(Z0))) 
                  Z0 <- apply(Z0, 2, function(z) {
                    if (all(is.na(z))) 
                      z <- rep(0, length(z))
                    mz <- mean(z, na.rm = TRUE)
                    z[is.na(z)] <- mz
                    z
                  })
                zvar <- apply(Z0, 2, var)
                mod1 <- coxlr.fit(cbind(Z0[, zvar != 0], X), 
                  nullmodel$y, nullmodel$strata, NULL, init = c(rep(0, 
                    ncol(Z0[, zvar != 0, drop = FALSE])), nullcoef), 
                  coxph.control(iter.max = 0), NULL, "efron", 
                  rn)
                mcov[which(!is.na(inds))[zvar != 0], which(!is.na(inds))[zvar != 
                  0]] <- tryCatch(ginv_s(mod1$var[1:sum(zvar != 
                  0), 1:sum(zvar != 0), drop = FALSE]), error = function(e) {
                  cov(Z0[nullmodel$y[, "status"] == 1, zvar != 
                    0, drop = FALSE]) * (sum(nullmodel$y[, "status"] == 
                    1) - 1)
                })
            }
            rownames(mcov) <- colnames(mcov) <- snp.names
            if (verbose) {
                assign("pb.i", get("pb.i", env) + 1, env)
                if (get("pb.i", env)%%ceiling(ngenes/100) == 
                  0) 
                  setTxtProgressBar(get("pb", env), get("pb.i", 
                    env))
            }
            return(forceSymmetric(Matrix(mcov, sparse = TRUE)))
        }), simplify = FALSE)
    for (k in 1:length(re)) {
        re[[k]] <- list(scores = zlrt[[k]] * sqrt(diag(as.matrix(re[[k]]))), 
            cov = as.matrix(re[[k]]), n = n, maf = maf[[k]], sey = 1)
    }
    if (verbose) 
        close(pb)
    return(re)
}
