% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_enrich2.r
\name{r2_enrich}
\alias{r2_enrich}
\title{r2_enrich}
\usage{
r2_enrich(dat, v1, v2, nv, exp1)
}
\arguments{
\item{dat}{N by (M+1) matrix having variables in the order of cbind(y,x)}

\item{nv}{sample size}

\item{exp1}{The expectation of the ratio (e.g. # SNPs for the genomic region of interest / total # SNPs in genomic partitioning)}

\item{v1/v2}{These can be set as v1=1 and v2=2, v1=2 and v2=1, v1=3 and v2=2, or any combination as long as the value is between 1 - M}
}
\description{
This function estimates var(t1/(t1+t2))
where t1 = R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])) and
t2 = R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2]))
where R2 is the R squared value of the model,
y is N by 1 matrix having the dependent variable, and
x is N by M matrix having M explanatory variables.
v1 or v2 indicates the ith column in the x matrix
(v1 or v2 should be a single interger between 1 - M, see Arguments below)
}
\examples{
To get test statistics for the ratio, i.e. t1/(t1+t2). 
t1 = R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])) and
t2 = R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2]))
(here we define R2_1=R2(y~x[,v1]), R2_2=R2(y~x[,v2]), R2_12=R2(y~x[,v1]+x[,v2])

dat=read.table("test_ukbb_enrichment_choles") (see example file)
nv=length(dat$V1)
v1=c(1)
v2=c(2)
expected_ratio=0.04 (# SNPs for the regulatory/total # SNPs)
output=r2_enrich(dat,v1,v2,nv,expected_ratio)
output

r2redux output

output$var1 (variance of R2_1)
8.758455e-05

output$var2 (variance of R2_2)
7.36385e-05

output$var12 (variance of R2_12)
0.000102236

output$var_diff1_2 (var of difference of R2(y~x[,v1]) - R2(y~x[,v2])))
6.074567e-05

output$var_diff12_1 (var of difference of R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])))
1.184853e-05

output$var_diff12_2 (var of difference of R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2])))
2.650564e-05

output$mean_diff12_1 (difference of R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])))
0.003048595

output$mean_diff12_2 (difference of  R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2])))
0.006845484

output$ratio (ratio =  t1/(t1+t2))
0.6918768

output$ratio_var (variance of ratio, var(t1/(t1+t2))
0.1324076

output$enrich_p (p-value for testing the ratio significantly different 
from the expectation (exp1))
0.07321821

output$upper_ratio (upper limit of 95\% CI for the ratio)
1.405079

output$lower_ratio (lower limit of 95\% CI for the ratio)
-0.02132515
}
\keyword{R2}
\keyword{between}
\keyword{of}
\keyword{ratio}
\keyword{variance}
