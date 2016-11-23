## CONfounding Factor Estimation Through Independent Component Analysis (CONFETI) R package

## Installation Instructions

### Dependencies

`confeti` requires the installation of the following two dependencies. 
Please follow the installation instructions to install both packages prior to installing `confeti`.

- R package `lrgpr` (http://lrgpr.r-forge.r-project.org)

- R package `picaplot` (https://github.com/jinhyunju/picaplot)


### Installation Instructions

After installing the above packages, run the following code to install `confeti`.

```r

install.packages("devtools")

library("devtools")

devtools::install_github("jinhyunju/confeti")

```

### Note 

A dependency issue is currently being worked out that causes an error in calling the `ica_genotype_test()` function. 
The temporary work around is to load the package `lrgpr` into the main workspace with `confeti`. I do understand that this is not an optimal solution and apologize for the inconvenience. I am currently trying to fix the issue, but in the mean time `lrgpr` will be listed in the `Depends` field in the DESCRIPTION file. 


## Usage

#### 1) Creating a sample covariance matrix using CONFETI

A confeti sample covariance matrix can be constructed by the function `confeti()`.

You will need two objects in your R environment.

- `expr_data` which is an `g x n` expression matrix with `g` gene measurements and `n` samples.

- `snp_data` which is a `s x n` genotype matrix with `s` genotypes coded as 0,1,2 and `n` samples.


```r

library(confeti)

confeti_results = confeti(expr_data, snp_data, return_all = FALSE)


```

`confeti_results` will be a list with 2 levels.

- `K_mx` : `n x n` sample covariance matrix calculated by confeti.

- `y_star` : Lower dimensional phenotype matrix with candidate genetic effects removed.


#### 2) Fitting a linear mixed model using `lrgpr`

After generating a sample covariance matrix using the `confeti()` function you can use it to fit a linear mixed model. 
For example, to fit an lmm for a single phenotype and single genotype you can use the `lrgpr()` function from the `lrgpr` package.

```r

single_pheno = expr_data[1,]

single_geno = snp_data[1,]

lmm_fit = lrgpr(single_pheno ~ single_geno, decomp = svd(confeti_results$Kmx))

```

