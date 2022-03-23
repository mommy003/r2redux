  #' r2_enrich_beta
  #'
  #' This function estimates var(t1/(t1+t2))
  #' where t1 = beta1^2 and t2 = beta2^2, and
  #' beta1 and 2 are regression coefficients from a multiple regression model,
  #' i.e. y = x1•beta1 + x2•beta2 + e, where y, x1 and x2 are column-standardised
  #' (see Olkin and Finn 1995).
  #' y is N by 1 matrix having the dependent variable, and
  #' x1 is N by 1 matrix having the ith explanatory variables.
  #' x2 is N by 1 matrix having the jth explanatory variables.
  #' v1 and v2 indicates the ith and jth column in the data
  #' (v1 or v2 should be a single interger between 1 - M, see Arguments below).
  #' Note that r2_enrich (above) and r2_enrich_beta is equivalent (identical p-value derived).
  #' @references
  #' Olkin, I. and J.D. Finn, Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
  #' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
  #' @param v1/v2 These can be set as v1=1 and v2=2, v1=2 and v2=1, v1=3 and v2=2, or any combination as long as the value is between 1 - M
  #' @param nv sample size
  #' @param exp1 The expectation of the ratio (e.g. ratio of # SNPs in genomic partitioning)
  #' @keywords variance of ratio between beta^2 from a multiple regression
  #' @export
  #' @examples
  #' To get the test statistic for the ratio which is significantly
  #' different from the expectation.
  #' var(t1/(t1+t2)), where t1 = beta1^2 and t2 = beta2^2.
  #' beta1 and beta2 are regression coefficients from a multiple regression model,
  #' i.e. y = x1•beta1 + x2•beta2 + e, where y, x1 and x2 are column-standardised
  #'
  #' dat=read.table("test_ukbb_enrichment_choles") (see example file)
  #' nv=length(dat$V1)
  #' v1=c(1)
  #' v2=c(2)
  #' expected_ratio=0.04
  #' output=r2_enrich_beta(dat,v1,v2,nv,expected_ratio)
  #' output
  #'
  #' r2redux output
  #'
  #' output$var1 (variance of t1)
  #' 7.072931e-05
  #' 
  #' output$var2 (variance of t2)
  #' 3.161929e-05
  #' 
  #' output$var1_2 (difference between t1 and t2)
  #' 0.000162113
  #' 
  #' output$beta1_sq (t1)
  #' 0.01118301
  #' 
  #' output$beta2_sq (t2)
  #' 0.004980285
  #' 
  #' output$cov (covariance between t1 and t2)
  #' -2.988221e-05
  #' 
  #' output$ratio (ratio = t1/(t1+t2_2))
  #' 0.6918768
  #' 
  #' output$ratio_var (variance of ratio, var(t1/(t1+t2))
  #' 0.1324076
  #' 
  #' output$enrich_p (p-value for testing the ratio significantly different from 
  #' the expectation (exp1))
  #' 0.07321821
  #' 
  #' output$upper_ratio (upper limit of 95% CI for the ratio)
  #' 1.405079
  #' 
  #' output$lower_ratio (lower limit of 95% CI for the ratio)
  #' -0.02132515


  r2_enrich_beta = function (dat,v1,v2,nv,exp1) {
    #source("aoa12_1.r")
    #source("aoa12_13.r")
    #source("aoa1_2.r")
    #source("aoa12_3.r")
    #source("aoa12_34.r")
    source("aoa_beta1_2.r")

    dat=scale(dat);omat=cor(dat)

      m1=lm(dat[,1]~dat[,(1+v1)])
      s1=summary(m1)
      m2=lm(dat[,1]~dat[,(1+v2)])
      s2=summary(m2)
      m3=lm(dat[,1]~dat[,1+v1]+dat[,1+v2])
      s3=summary(m3)

      ord=c(1,(1+v1),(1+v2))
      aoa=olkin_beta1_2(omat[ord,ord],nv)
      var1=aoa$var1
      var2=aoa$var2
      var1_2=aoa$var1_2

      cov=-0.5*(var1_2-var1-var2)

      dvr1=s3$coefficients[2,1]^2
      dvr2=s3$coefficients[3,1]^2
      #dvr1=s3$coefficients[2,1]
      #dvr2=s3$coefficients[3,1]

      #enrichment p-value
      dvrt=dvr1+dvr2; ratio1=dvr1/dvrt; ratio2=dvr2/dvrt

      vart=var1+var2+2*cov
      covt1=var1+cov
      covt2=var2+cov

      ratio1_var=(ratio1^2)*(var1/dvr1^2+vart/dvrt^2-2*covt1/(dvr1*dvrt))
      ratio2_var=(ratio2^2)*(var2/dvr2^2+vart/dvrt^2-2*covt2/(dvr2*dvrt))

      chisq1=((ratio1-exp1)^2)/ratio1_var
      p3=pchisq(chisq1,1,lower.tail=F)

      #95% CI
      uci=ratio1+1.96*ratio1_var^.5
      lci=ratio1-1.96*ratio1_var^.5

      #z=list(var1=var1,var2=var2,var1_2=var1_2,beta1_sq=dvr1,beta2_sq=dvr2,cov=cov,ratio1=ratio1,ratio2=ratio2,ratio1_var=ratio1_var,ratio2_var=ratio2_var,enrich_p=p3,upper_ratio1=uci,lower_ratio1=lci)
      z=list(var1=var1,var2=var2,var1_2=var1_2,beta1_sq=dvr1,beta2_sq=dvr2,cov=cov,ratio=ratio1,ratio_var=ratio1_var,enrich_p=p3,upper_ratio=uci,lower_ratio=lci)
      return(z)


  }



