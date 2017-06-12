#' Main function to run CONFETI.
#'
#' Creating a sample covariance matrix based on non-genetic independent components.
#'
#' @param expr_data Expression matrix with dimensions g x n.
#' @param snp_data Genotype matrix with dimensions s x n.
#' @param return_all If TRUE, the full ICA results are returned. If set FALSE
#' only the sample covariance matrix and lower dimensional phenotype matrix are returned.
#' @param k Number of components to be estimated or method to estimate it.
#' @param var_cutoff Percent variance threshold to use when <k_est> is not supplied.
#' @param ica_runs Number of times to run ICA. If this is set to a number larger than 1,
#' only ICs that replicate between runs are going to be returned
#' @param cores Number of cores to use for genotype association testing, and ICA multi runs (if ica_runs \> 1).
#' @param scale_pheno If set to TRUE the pre-processing step of the data will include a scaling step.
#' This will divide all phenotypes by their standard deviation.
#' @param h_clust_cutoff is the cutoff value used in hierarchical clustering. Default is set to 0.3.
#' @param similarity_measure How to measure the similarity between ICs.
#' @param threshold Bonferroni significance threshold for genotype association testing
#' @param return_pval If set to TRUE, p-values for genotype ICA association testing are returned
#'
#' @return N x N covariance structure matrix.
#' @keywords keywords
#'
#' @import lrgpr
#' @import picaplot
#' @import bigmemory
#' @export
confeti <- function(expr_data, snp_data, return_all = TRUE,
                    cores = 1,
                    k = NULL, var_cutoff = 95,
                    ica_runs = 1, scale_pheno = FALSE,
                    h_clust_cutoff = 0.3,
                    similarity_measure = "peaks",
                    threshold = 0.05,
                    return_pval = FALSE){

    ica_object = picaplot::runICA(pheno_mx = expr_data,
                         k_est = k,
                         var_cutoff = var_cutoff,
                         n_runs = ica_runs,
                         n_cores = cores,
                         scale_pheno = scale_pheno,
                         h_clust_cutoff = h_clust_cutoff,
                         similarity_measure = similarity_measure)

    ica_object = ica_genotype_test(ica_object = ica_object,
                                   genotype_mx = t(snp_data),
                                   n_cores = cores,
                                   sig_threshold = threshold,
                                   return_pval = return_pval)
    confeti_results = get_similarity_mx(ica_object)

    if(return_all){
        ica_object$confeti = confeti_results
        return(ica_object)
    } else {
        return(confeti_results)
    }

}

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
#' @import picaplot
#' @import bigmemory
#' @export
ica_genotype_test <- function(ica_object,
                              genotype_mx,
                              n_cores = 1,
                              sig_threshold = 0.05,
                              return_pval = FALSE){
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
#' @param normalize_ystar If set to TRUE phenotype values of the reconstructed y will be mean centered and divided by their standard deviation.
#' @return N x N covariance structure matrix.
#' @keywords keywords
#'
#' @export
get_similarity_mx <- function(ica_object,
                              normalize_ystar = FALSE){

    cf_idx = ica_object$confounding_factors

	message(length(cf_idx), " out of ", nrow(ica_object$A), " ICs used as confounding factors \n")

	if(length(cf_idx) == 0){
		stop("No ICs estimated as hidden factors")
	}

	non_genetic_factors <- match(ica_object$confounding_factors, colnames(ica_object$S))
	y_star <- ica_object$S[,non_genetic_factors] %*% ica_object$A[non_genetic_factors,]

	if(normalize_ystar){
	    y_norm = t(apply(y_star, 1, function(x) (x-mean(x)) / sd(x)))
	    K_mx = cov(y_norm)
	} else {
	    K_mx <- cov(y_star)
	}

    return(list("Kmx" = K_mx, "y_star" = y_star))
}



