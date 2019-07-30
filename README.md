
<!-- README.md is generated from README.Rmd. Please edit that file -->
modi
====

The package **modi** provides several functions for multivariate outlier detection and imputation. They may be used when analysing multivariate quantitative survey data. The distribution of such data is often not multivariate normal. Futhermore, the data is often skewed and exhibits features of a semi-continuous distribution. Finally, missing values and non-response is common. The functions provided in **modi** address those problems.

Overview
--------

The following outlier detection and imputation functions are provided in **modi**:

-   `BEM()` is an implementation of the BACON-EEM algorithm to detect outliers under the assumption of a multivariate normal distribution.
-   `TRC()` is an implementation of the transformed rank correlation (TRC) algorithm to detect outliers.
-   `EAdet()` is an outlier detection method based on the epidemic algorithm (EA).
-   `EAimp()` is an imputation method based on the epidemic algorithm (EA).
-   `GIMCD()` is an outlier detection method based on non-robust Gaussian imputation (GI) and the highly robust minimum covariance determinant (MCD) algorithm.
-   `POEM()` is a nearest neighbor imputation method for outliers and missing values.
-   `winsimp()` is an imputation method for outliers and missing values based on winsorization and Gaussian imputation.
-   `ER()` is a robust multivariate outlier detection algorithm that can cope with missing values.

Installation
------------

You can install modi from github by runing the following line of code in R:

``` r
# install.packages("devtools")
devtools::install_github("martinSter/modi")
```

Example
-------

The following simple example shows how the BACON-EEM algorithm can be applied to detect outliers in the Bushfire dataset:

``` r
# Bushfire data set with 20% of observations missing completely at random (MCAR)

library(modi)

# load the data containing missing values and the corresponding weights
data(bushfirem, bushfire.weights)

# run BEM algorithm to detect multivariate outliers
bem.res <- BEM(bushfirem, bushfire.weights, alpha = (1 - 0.01 / nrow(bushfirem)))
#> alpha should be less than 0.5:
#>         alpha set to 1 - alpha
#>  0.9997368
#>  BEM has detected 15 outlier(s) in 0.07 seconds. 
#> 

# show the outliers as detected by BEM
bushfirem[bem.res$outind, ]
#>     X1  X2  X3  X4  X5
#> 7   92 110  46 165  NA
#> 8   94  95  29 113 190
#> 9   94  94  29 110 188
#> 10 100 104  21 133 208
#> 11 108 115  17 144  NA
#> 12 134 156  10 163 233
#> 15  81 137 426  NA 306
#> 31 136 155  NA 246 301
#> 32 103  97 552  NA 364
#> 33  NA  66 576 340 377
#> 34  79  66 572 340 376
#> 35  79  66 577 341 379
#> 36  78  NA 574 342 377
#> 37  78  66 571 343 379
#> 38  NA  66  NA 344 380

# show mean per column as computed in BEM
print(bem.res$output$center)
#> [1] 109.2987 149.6075 262.0240 215.3746 276.8635
```
