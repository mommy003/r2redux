# r2redux
The r2redux package will be used to test the significant difference between two PRS. We use predictive ability (R^2) measures and their variance covariance matrix to assess if the predictive ability of PRS based on different sources is significantly different to each other. 

# INSTALLATION
To use r2redux:
- install.packages("devtools")
- install.packages("r2redux")
- library(devtools)
- library(r2redux)
 
# QUICK START
We illustrate the usage of r2redux2 using the GWAS summary statistics from UK Biobank (reference) to predict white British as target. Note that the target individuals were independent from reference individuals. We can test the significance differences of predictive ability among different P-value threshold of SNP effect. We can also test the differences of predictive ability of UKBB reference and BBJ reference in target individuals or we can test the differences of joint model and single model [e.g., R_((UKBB+BBJ))^2 vs  R_((UKBB))^2  and/or R_((UKBB+BBJ))^2 and R_((BBJ))^2]. 


# DATA PREPARATION
a.	To estimate R2 for each p-value threshold 
r2redux requires only phenotype and estimated PRS from PLINK or any other software of interest. Please note that, any missing values in the phenotypes should be removed. Phenotype and PRSs should be scaled before using r2redux2. If we want to test significant difference between/among thresholds, need to prepare input file for r2redux2 that includes following fields (e.g. test_ukbb_thresholds_scaled in example file). 
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

b. Genomic enrichment analysis
If we want to perform some enrichment analysis (regulatory vs non_regulatory) to in the PRS context to test significantly different from the expectation (e.g., contribution of regulatory SNP is 4%). We simultaneously fit two sets of PRS from regulatory and non-regulatory to get β ̂_regu^2 and β ̂_(non_regu)^2, using a multiple regression, and assess if the ratio, (β ̂_regu^2)/(β ̂_regu^2+β ̂_(non_regu)^2 ) , is significantly different from the expectation. To test this, we need to prepare input file for r2redux2 that includes following fields (e.g. test_ukbb_enrichment_choles in example file).
- Phenotype (y)
- PRS for regulatory region (x1)
- PRS for non-regulatory region (x2)      

