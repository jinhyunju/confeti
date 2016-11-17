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

## Usage

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



