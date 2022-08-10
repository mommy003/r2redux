# r2redux

The ‘r2redux’ package can be used to derive test statistics for R^2 values from polygenic risk score (PGS) models (variance and covariance of R^2 values, p-value and 95% confidence intervals (CI)). For example, it can test if two sets of R^2 values from two different PGS models are significantly different to each other whether the two sets of PGS are independent or dependent. Because R^2 value is often regarded as the predictive ability of PGS, r2redux package can be useful to assess the performances of PGS methods or multiple sets of PGS based on different information sources. Furthermore, the package can derive the information matrix of β ̂_1^2and β ̂_2^2 from a multiple regression (see olkin_beta1_2 or olkin_beta_info function in the manual), which is a basis of a novel PGS-based genomic partitioning method (see r2_enrich or r2_enrich_beta function in the manual). It is recommended that the target sample size in the PGS study should be more than 2,000 for quantitative traits  and more than 5,000 for binary responses or case-control studies. The P-value generated from the r2redux package provides two types of p-values (for one- and two-tailed test) unless the comparison is for nested models (e.g. y=PGS_1+PGS_2+e vs. y=PGS_2+e) where the R^2 of the full model is expected to be always higher than the reduced model.  When there are multiple covariates (e.g. age, sex and other demographic variables), the phenotypes can be adjusted for the covariates, and pre-adjusted phenotypes (residuals) should be used in the r2redux.  

# INSTALLATION
To use r2redux:
```
install.packages("r2redux") 
library(r2redux)
```
 or
```
install.packages("devtools")
library(devtools)
devtools::install_github("mommy003/r2redux")
library(r2redux)
```


# QUICK START
We illustrate the usage of r2redux using multiple sets of PGS estimated based on GWAS summary statistics from UK Biobank or Biobank Japan (reference datasets). In a target dataset, the phenotypes of target samples (y) can be predicted with PGS (a PGS model, e.g. y=PRS+e, where y and PGS are column-standardised 1. Note that the target individuals should be independent from reference individuals. We can test the significant differences of the predictive ability (R^2) between a pair of PGS (see r2_diff function and example in the manual).

# DATA PREPARATION
**a.	Statistical testing of significant difference between R2 values for p-value thresholds:** 
r2redux requires only phenotype and estimated PGS (from PLINK or any other software). Note that any missing value in the phenotypes and PGS tested in the model should be removed. If we want to test the significant difference of R^2 values for p-value thresholds, r2_diff function can be used with an input file that includes the following fields (also see test_ukbb_thresholds_scaled in the example directory form github  (https://github.com/mommy003/r2redux)  or read dat1 file embedded within the package and r2_diff function in the manual).


- Phenotype (y)
- PGS for p value 1 (x1)
- PGS for p value 0.5 (x2)
- PGS for p value 0.4 (x3)
- PGS for p value 0.3 (x4)
- PGS for p value 0.2 (x5)
- PGS for p value 0.1 (x6)
- PGS for p value 0.05 (x7)
- PGS for p value 0.01 (x8)
- PGS for p value 0.001 (x9)
- PGS for p value 0.0001 (x10)

To get the test statistics for the difference between R2(y=x[,v1]) and R2(yx[,v2]). (here we define R_1^2= R^2(y=x[,v1])) and R_2^2=R^2(y=x[,v2])))
```
dat=read.table("test_ukbb_thresholds_scaled") #(see example files) or
dat=dat1 #(this example embedded within the package)
nv=length(dat$V1)
v1=c(1)
v2=c(2)
output=r2_diff(dat,v1,v2,nv)
```

- r2redux output
- output$var1 (variance of R_1^2)
- 0.0001437583
- output$var2 (variance of R_2^2)
- 0.0001452828
- output$var_diff (variance of difference between R_1^2and R_2^2)
- 5.678517e-07
- output$r2_based_p (p-value for significant difference between R_1^2  and R_2^2)
- 0.5514562
- output$mean_diff (differences between R_1^2 and R_2^2)
- -0.0004488044
- output$upper_diff (upper limit of 95% CI for the difference)
- 0.001028172
- output$lower_diff (lower limit of 95% CI for the difference)
- -0.001925781

 
**b. PGS-based genomic enrichment analysis:**
If we want to perform some enrichment analysis (e.g., regulatory vs non_regulatory) in the PGS context to test significantly different from the expectation (p_exp= # SNPs in the regulatory / total # SNPs = 4%). We simultaneously fit two sets of PGS from regulatory and non-regulatory to get β ̂_regu^2 and β ̂_(non-regu)^2, using a multiple regression, and assess if the ratio,  (β ̂_1^2)/(r_(y,〖(x〗_1,x_2))^2 ) are significantly different from the expectation, p_exp. To test this, we need to prepare input file for r2redux that includes the following fields (e.g. test_ukbb_enrichment_choles in example directory or read dat2 file embedded within the package and r2_enrich_beta function in the manual).

- Phenotype (y)
- PGS for regulatory region (x1)
- PGS for non-regulatory region (x2)      

To get the test statistic for the ratio which is significantly different from the expectation. var(β ̂_1^2/r_(y,(x_1,x_2))^2), where β ̂_1^2  is the squared regression coefficient of x_1 from a multiple regression model, i.e. y=x_1 β_1+ x_2 β_2+e, and r_(y,(x_1,x_2))^2 is the coefficient of determination of the model. It is noted that y, x_1 and x_2 are column standardised (mean 0 and variance 1).

- dat=read.table("test_ukbb_enrichment_choles") (see example file) or 
- dat=dat2 (this example embedded within the package)
- nv=length(dat$V1)
- v1=c(1)
- v2=c(2)
- dat=dat2
- nv=length(dat$V1)
- v1=c(1)
- v2=c(2)
- output=r2_beta_var(dat,v1,v2,nv)
- r2redux output
- output$beta1_sq (beta1^2)
- 0.01118301
- output$beta2_sq (beta2^2)
- 0.004980285
- output$var1 (variance of beta1^2)
- 7.072931e-05
- output$var2 (variance of beta2^2)
- 3.161929e-05
- output$var1_2 (variance of difference between beta1^2 and beta2^2)
- 0.000162113
- output$cov (covariance between beta1^2 and beta2^2)
- -2.988221e-05
- output$upper_beta1_sq (upper limit of 95% CI for beta1^2)
- 0.03037793
- output$lower_beta1_sq (lower limit of 95% CI for beta1^2)
- -0.00123582
- output$upper_beta2_sq (upper limit of 95% CI for beta2^2)
- 0.02490076
- output$lower_beta2_sq (lower limit of 95% CI for beta2^2)
- -0.005127546

- dat=dat2 (this example embedded within the package)
- nv=length(dat$V1)
- v1=c(1)
- v2=c(2)
- expected_ratio=0.04
- output=r2_enrich_beta(dat,v1,v2,nv,expected_ratio)
- r2redux output
- output$beta1_sq (beta1^2)
- 0.01118301
- output$beta2_sq (beta2^2)
- 0.004980285
- output$ratio1 (beta1^2/R^2)
- 0.4392572
- output$ratio2 (beta2^2/R^2)
- 0.1956205
- output$ratio_var1 (variance of ratio 1)
- 0.08042288
- output$ratio_var2 (variance of ratio 2)
- 0.0431134
- output$upper_ratio1 (upper limit of 95% CI for ratio 1)
- 0.9950922
- output$lower_ratio1 (lower limit of 95% CI for ratio 1)
- -0.1165778
- output$upper_ratio2 upper limit of 95% CI for ratio 2)
- 0.6025904
- output$lower_ratio2 (lower limit of 95% CI for ratio 2)
- -0.2113493
- output$enrich_p1 (two tailed P-value for beta1^2/R^2 is significantly different from exp1)
- 0.1591692
- output$enrich_p1_one_tail (one tailed P-value for beta1^2/R^2 is significantly different from exp1)
- 0.07958459
- output$enrich_p2 (P-value for beta2^2/R2 is significantly different from (1-exp1))
- 0.000232035
- output$enrich_p2_one_tail (one tailed P-value for beta2^2/R2 is significantly different from (1-exp1))
- 0.0001160175


# References
1. Olkin, I. and  Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
2. Momin, M.M., Lee, S., Wray, N.R. and Lee S.h. 2022. Significance tests for R2 of out-of-sample prediction using polygenic scores. bioRxiv.

# Contact information
Please contact Hong Lee (hong.lee@unisa.edu.au) or Moksedul Momin (cvasu.momin@gmail.com) if you have any queries.
