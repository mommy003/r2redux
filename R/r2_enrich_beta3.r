#' r2_enrich_beta
#'
#' This function estimates var(beta1^2/R^2), 
#' beta1 and R^2 are regression coefficient and the coefficient
#' of determination from a multiple regression model,
#' i.e. y = x1 * beta1 + x2 * beta2 +e, where y, x1 and x2 are column-standardised
#' (see Olkin and Finn 1995).
#' y is N by 1 matrix having the dependent variable, and
#' x1 is N by 1 matrix having the ith explanatory variables.
#' x2 is N by 1 matrix having the jth explanatory variables.
#' v1 and v2 indicates the ith and jth column in the data
#' (v1 or v2 should be a single interger between 1 - M, see Arguments below).
#' @references
#' Olkin, I. and Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
#' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
#' @param v1 These can be set as v1=1, v1=2, v1=3  or any value between 1 - M based on combination
#' @param v2 These can be set as v2=1, v2=2, v2=3, or any value between 1 - M based on combination 
#' @param nv Sample size
#' @param exp1 The expectation of the ratio (e.g. ratio of # SNPs in genomic partitioning)
#' @keywords variance of ratio between beta^2 from a multiple regression
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @return  This function will estimate var(beta1^2/R^2), beta1 and R^2 are regression coefficient and the coefficient of determination from a multiple regression model, i.e. y = x1 * beta1 + x2 * beta2 +e, where y, x1 and x2 are column-standardised. The outputs are listed as follows.
#' \item{beta1_sq}{beta1_sq}
#' \item{beta2_sq}{beta2_sq}
#' \item{ratio1}{beta1_sq/R^2}
#' \item{ratio2}{beta2_sq/R^2}
#' \item{ratio_var1}{variance of ratio 1}
#' \item{ratio_var2}{variance of ratio 2}
#' \item{upper_ratio1}{upper limit of 95% CI for ratio 1}
#' \item{lower_ratio1}{lower limit of 95% CI for ratio 1}
#' \item{upper_ratio2}{upper limit of 95% CI for ratio 2}
#' \item{lower_ratio2}{lower limit of 95% CI for ratio 2}
#' \item{enrich_p1}{two tailed P-value for beta1_sq/R^2 is significantly different from exp1}
#' \item{enrich_p1_one_tail}{one tailed P-value for beta1_sq/R^2 is significantly different from exp1}
#' \item{enrich_p2}{P-value for beta2_sq/R2 is significantly different from (1-exp1)}
#' \item{enrich_p2_one_tail}{one tailed P-value for beta2_sq/R2 is significantly different from (1-exp1)}
#' @examples
#' #To get the test statistic for the ratio which is significantly
#' #different from the expectation, this function estiamtes 
#' #var (beta1^2/R^2), where 
#' #beta1^2 and R^2 are regression coefficients and the 
#' #coefficient of dterminationfrom a multiple regression model,
#' #i.e. y = x1 * beta1 + x2 * beta2 +e, where y, x1 and x2 are 
#' #column-standardised.
#'
#' dat=dat2
#' nv=length(dat$V1)
#' v1=c(1)
#' v2=c(2)
#' expected_ratio=0.04
#' output=r2_enrich_beta(dat,v1,v2,nv,expected_ratio)
#' output
#' 
#' #r2redux output
#'
#' #output$beta1_sq (beta1_sq)
#' #0.01118301
#' 
#' #output$beta2_sq (beta2_sq)
#' #0.004980285
#' 
#' #output$ratio1 (beta1_sq/R^2)
#' #0.4392572
#' 
#' #output$ratio2 (beta2_sq/R^2)
#' #0.1956205
#' 
#' #output$ratio_var1 (variance of ratio 1)
#' #0.08042288
#' 
#' #output$ratio_var2 (variance of ratio 2)
#' #0.0431134
#' 
#' #output$upper_ratio1 (upper limit of 95% CI for ratio 1)
#' #0.9950922
#' 
#' #output$lower_ratio1 (lower limit of 95% CI for ratio 1)
#' #-0.1165778
#' 
#' #output$upper_ratio2 upper limit of 95% CI for ratio 2)
#' #0.6025904
#' 
#' #output$lower_ratio2 (lower limit of 95% CI for ratio 2)
#' #-0.2113493
#' 
#' #output$enrich_p1 (two tailed P-value for beta1_sq/R^2 is 
#' #significantly different from exp1)
#' #0.1591692
#' 
#' #output$enrich_p1_one_tail (one tailed P-value for beta1_sq/R^2 
#' #is significantly different from exp1)
#' #0.07958459
#' 
#' #output$enrich_p2 (two tailed P-value for beta2_sq/R2 is 
#' #significantly different from (1-exp1))
#' #0.000232035
#' 
#' #output$enrich_p2_one_tail (one tailed P-value for beta2_sq/R2  
#' #is significantly different from (1-exp1))
#' #0.0001160175


r2_enrich_beta = function (dat,v1,v2,nv,exp1) {
  
  dat=scale(dat);omat=cor(dat)
  m3=lm(dat[,1]~dat[,1+v1]+dat[,1+v2])
  s3=summary(m3)
  dvr1=s3$coefficients[2,1]^2
  dvr2=s3$coefficients[3,1]^2
  
  #enrichment p-value
  dvrt=s3$r.squared
  ratio1=dvr1/dvrt 
  ratio2=dvr2/dvrt
  ord=c(1,(1+v1),(1+v2))
  
  aoa1=olkin_beta_ratio(omat[ord,ord],nv)
  ratio1_var=aoa1$ratio_var
  
  chisq1=((ratio1-exp1)^2)/ratio1_var
  p3=pchisq(chisq1,1,lower.tail=F)
  
  #95% CI
  uci1=(ratio1)+1.96*ratio1_var^.5
  lci1=(ratio1)-1.96*ratio1_var^.5
  
  ord=c(1,(1+v2),(1+v1))
  
  aoa2=olkin_beta_ratio(omat[ord,ord],nv)
  ratio2_var=aoa2$ratio_var
  
  chisq2=((ratio2-(1-exp1))^2)/ratio2_var
  p4=pchisq(chisq2,1,lower.tail=F)
  
  #95% CI
  uci2=(ratio2)+1.96*ratio2_var^.5
  lci2=(ratio2)-1.96*ratio2_var^.5
  z=list(beta1_sq=dvr1,beta2_sq=dvr2,ratio1=ratio1,ratio2=ratio2, ratio_var1=ratio1_var,ratio_var2=ratio2_var,upper_ratio1=uci1,lower_ratio1=lci1, upper_ratio2=uci2,lower_ratio2=lci2,enrich_p1=p3,enrich_p1_one_tail=p3/2, enrich_p2=p4,enrich_p2_one_tail=p4/2)
  
  #NOTE: enrich_p1=p3, enrichment_p2=p4 due to normal distribution
  return(z)
  

} 