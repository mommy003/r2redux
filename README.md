# r2redux
The ‘r2redux’ package can be used to derive test statistics for R2 values from polygenic risk score (PRS) models (variance and covariance of R2 values, p-value and 95% confidence intervals (CI)). For example, it can test if two sets of R2 values from two different PRS models are significantly different to each other whether the two sets of PRS are independent or dependent. Because R2 value is often regarded as the predictive ability of PRS, r2redux package can be useful to assess the performances of PRS methods or multiple sets of PRS based on different information sources. Furthermore, the package can derive the information matrix of beta1^2 and beta2^2 from a multiple regression (see olkin_beta1_2 or olkin_beta_info function in the manual), which is a basis of a novel PRS-based genomic partitioning method (see r2_enrich or r2_enrich_beta function in the manual).  

# INSTALLATION
To use r2redux:
- install.packages("devtools")
- library(devtools)
- devtools::install_github("mommy003/r2redux") or
- install.packages("r2redux")  
- library(r2redux)

# QUICK START
We illustrate the usage of r2redux using multiple sets of PRS estimated based on GWAS summary statistics from UK Biobank or biobank Japan (reference datasets). In a target dataset, the phenotypes of target samples (y) can be predicted with PRS (a PRS model, e.g. y = PRS + e where y and PRS are column-standardised (Olkin and Finn 1995)). Note that the target individuals should be independent from reference individuals. We can test the significant differences of the predictive ability (R2) between a pair of PRS (see r2_diff function and example in the manual). 


# DATA PREPARATION
**a.	Statistical testing of significant difference between R2 values for p-value thresholds:** 
r2redux requires only phenotype and estimated PRS (from PLINK or any other software). Note that any missing value in the phenotypes should be removed. Phenotype and PRSs should be column-standardised before using r2redux (Olkin and Finn 1995). If we want to test the significant difference of R2 values for p-value thresholds, r2_diff function can be used with an input file that includes the following fields (also see test_ukbb_thresholds_scaled in the example directory and r2_diff function in the manual). 
- Phenotype (y)
- PRS for p value 1 (x1)
- PRS for p value 0.5 (x2)
- PRS for p value 0.4 (x3)
- PRS for p value 0.3 (x4)
- PRS for p value 0.2 (x5)
- PRS for p value 0.1 (x6)
- PRS for p value 0.05 (x7)
- PRS for p value 0.01 (x8)
- PRS for p value 0.001 (x9)
- PRS for p value 0.0001 (x10)
 
**b. PRS-based genomic enrichment analysis:**
If we want to perform some enrichment analysis (e.g., regulatory vs non_regulatory) in the PRS context to test significantly different from the expectation (4% = # SNPs in the regulatory / total # SNPs). We simultaneously fit two sets of PRS from regulatory and non-regulatory to get Î²_regu^2 and Î²_non_regu^2, using a multiple regression, and assess if the ratio, (Î²_regu^2)/(Î²_regu^2 + Î²_(non_regu)^2 ) , is significantly different from the expectation. To test this, we need to prepare input file for r2redux that includes the following fields (e.g. test_ukbb_enrichment_choles in example directory and r2_enrich or r2_enrich_beta function in the manual).
- Phenotype (y)
- PRS for regulatory region (x1)
- PRS for non-regulatory region (x2)      

# Refrences
Olkin, I. and J.D. Finn, Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.

# Contact information
Please contact Hong Lee (hong.lee@unisa.edu.au) or Moksedul Momin (momin@cvasu.ac.bd) if you have any queries.
