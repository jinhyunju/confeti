#' Function that uses lrgpr's glmApply to test associations between IC coefficients and genotypes.
#'
#' Testing associations between IC coefficients and genotypes.
#'
#' @param ica_object ICA object that is created with the `run_ica()` function from the `picaplot` package.
#' @param genotype_mx Input genotype matrix with N samples x G genotypes.
#' @param n_cores Number of cores to use for this process.
#' @param sig_threshold Bonferroni significance threshold to use in calling significant p-values.
#' @param return_pval TRUE or FALSE value indicating whether the calculated p-value matrix should be returned.
#' @return ICA object with confounding factor information added.
#' @keywords keywords
#'
#' @import lrgpr
#' @import formula.tools
#' @export
ica_genotype_test <- function(ica_object, genotype_mx, n_cores = 1, sig_threshold = 0.05, return_pval = FALSE){
    ic_coefficients <- t(ica_object$A)

    pval_mx <- glmApply(ic_coefficients ~ SNP,
                        features = genotype_mx,
                        nthreads = n_cores)$pValues

    colnames(pval_mx) <- rownames(ica_object$A)
    sig <- which(pval_mx < (sig_threshold/length(pval_mx) ), arr.ind = TRUE)

    genetic_factors <- colnames(pval_mx)[unique(sig[,"col"])]
    non_genetic <- colnames(pval_mx)[which(!(colnames(pval_mx) %in% genetic_factors))]

    if(return_pval){
        ica_object$geno_pval = pval_mx
        ica_object$genetic_factors = genetic_factors
        ica_object$confounding_factors = non_genetic
        return(ica_object)
    } else {
        ica_object$genetic_factors = genetic_factors
        ica_object$confounding_factors = non_genetic
        return(ica_object)
    }
}


#' Function to calculate the covariance structure (K) matrix.
#'
#' Creating the K matrix based on non-genetic independent components.
#'
#' @param ica_object ICA object that contains the confounding factor labels generated with `ica_genotype_test`.
#' @return N x N covariance structure matrix.
#' @keywords keywords
#'
#' @export
get_similarity_mx <- function(ica_object){

    cf_idx = ica_object$confounding_factors

	message(length(cf_idx), " out of ", nrow(ica_object$A), " ICs used as confounding factors \n")

	if(length(cf_idx) == 0){
		stop("No ICs estimated as hidden factors")
	}

	non_genetic_factors <- match(ica_object$confounding_factors, colnames(ica_object$S))
	y_star <- ica_object$S[,non_genetic_factors] %*% ica_object$A[non_genetic_factors,]
	K_mx <- cov(y_star)
    return(K_mx)
}
