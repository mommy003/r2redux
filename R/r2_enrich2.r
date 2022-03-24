  #' r2_enrich
  #'
  #' This function estimates var(t1/(t1+t2)) 
  #' where t1 = R2(y~x\[,v1]+x\[,v2]) - R2(y~x\[,v1])) and
  #' t2 = R2(y~x\[,v1]+x\[,v2]) - R2(y~x\[,v2]))
  #' where R2 is the R squared value of the model,
  #' y is N by 1 matrix having the dependent variable, and
  #' x is N by M matrix having M explanatory variables.
  #' v1 or v2 indicates the ith column in the x matrix
  #' (v1 or v2 should be a single interger between 1 - M, see Arguments below)
  #' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
  #' @param v1/v2 These can be set as v1=1 and v2=2, v1=2 and v2=1, v1=3 and v2=2, or any combination as long as the value is between 1 - M
  #' @param nv sample size
  #' @param exp1 The expectation of the ratio (e.g. # SNPs for the genomic region of interest / total # SNPs in genomic partitioning) 
  #' @keywords variance of ratio between R2
  #' @export
  #' @examples 
  #' To get test statistics for the ratio, i.e. t1/(t1+t2). 
  #' t1 = R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])) and
  #' t2 = R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2]))
  #' (here we define R2_1=R2(y~x[,v1]), R2_2=R2(y~x[,v2]), R2_12=R2(y~x[,v1]+x[,v2])
  #'
  #' dat=read.table("test_ukbb_enrichment_choles") (see example file)
  #' nv=length(dat$V1)
  #' v1=c(1)
  #' v2=c(2)
  #' expected_ratio=0.04 (# SNPs for the regulatory/total # SNPs)
  #' output=r2_enrich(dat,v1,v2,nv,expected_ratio)
  #' output
  #'
  #' r2redux output
  #'
  #' output$var1 (variance of R2_1)
  #' 8.758455e-05
  #' 
  #' output$var2 (variance of R2_2)
  #' 7.36385e-05
  #' 
  #' output$var12 (variance of R2_12)
  #' 0.000102236
  #' 
  #' output$var_diff1_2 (var of difference of R2(y~x[,v1]) - R2(y~x[,v2])))
  #' 6.074567e-05
  #' 
  #' output$var_diff12_1 (var of difference of R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])))
  #' 1.184853e-05
  #' 
  #' output$var_diff12_2 (var of difference of R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2])))
  #' 2.650564e-05
  #' 
  #' output$mean_diff12_1 (difference of R2(y~x[,v1]+x[,v2]) - R2(y~x[,v1])))
  #' 0.003048595
  #' 
  #' output$mean_diff12_2 (difference of  R2(y~x[,v1]+x[,v2]) - R2(y~x[,v2])))
  #' 0.006845484
  #' 
  #' output$ratio (ratio =  t1/(t1+t2))
  #' 0.6918768
  #' 
  #' output$ratio_var (variance of ratio, var(t1/(t1+t2))
  #' 0.1324076
  #' 
  #' output$enrich_p (p-value for testing the ratio significantly different 
  #' from the expectation (exp1))
  #' 0.07321821
  #' 
  #' output$upper_ratio (upper limit of 95% CI for the ratio)
  #' 1.405079
  #' 
  #' output$lower_ratio (lower limit of 95% CI for the ratio)
  #' -0.02132515




  r2_enrich = function (dat,v1,v2,nv,exp1) {
    source("aoa12_1.r")
    source("aoa12_13.r")
    source("aoa1_2.r")
    source("aoa12_3.r")
    source("aoa12_34.r")

    dat=scale(dat);omat=cor(dat)

      m1=lm(dat[,1]~dat[,(1+v1)])
      s1=summary(m1)
      m2=lm(dat[,1]~dat[,(1+v2)])
      s2=summary(m2)
      m3=lm(dat[,1]~dat[,1+v1]+dat[,1+v2])
      s3=summary(m3)

      R2=s1$r.squared;mv2=1    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100  #*(1-R2)

      R2=s2$r.squared;mv2=1    #expected variance for s1r2
      t200=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t200           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t200=t200^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var2=t200  #*(1-R2)

      R2=s3$r.squared;mv2=3    #expected variance for s1r2
      t300=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t300           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t300=t300^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var12=t300  #*(1-R2)

      dvr1=s3$r.squared-s1$r.squared
      chi_dum=dvr1/(1/(nv-1)*(1-dvr1)^2) #NCP
      p1=pchisq(chi_dum,1,lower.tail=F)

      dvr2=s3$r.squared-s2$r.squared
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      ord=c(1,(1+v1),(1+v2))
      #omat=omat[ord,ord]
      aoa1=olkin12_1(omat[ord,ord],nv)
      aoa0=olkin1_2(omat[ord,ord],nv)

      ord=c(1,(1+v2),(1+v1))
      #omat=omat[ord,ord]
      aoa2=olkin12_1(omat[ord,ord],nv)

      # cov (s3r2-s1r2,s3r2-s2r2) ***********************************************
      # for mean(res[,3]*res[,3])
      v300=s3$r.squared^2+t300
      # for mean(res[,2]*res[,1])=mean()*mean()+cov(1,2)
      v300=v300+s2$r.squared*s1$r.squared+(t100+t200-aoa0)/2
      #for mean(res[,3]*res[,2])
      v300=v300-(s3$r.squared*s2$r.squared+(t300+t200-aoa2)/2)
      #for mean(res[,3]*res[,1])
      v300=v300-(s3$r.squared*s1$r.squared+(t300+t100-aoa1)/2)

      v102=dvr2*dvr1
      cov=v300-v102

      #enrichment p-value
      dvrt=dvr1+dvr2; ratio1=dvr1/dvrt; ratio2=dvr2/dvrt

      vart=aoa1+aoa2+2*cov
      covt1=aoa1+cov
      covt2=aoa2+cov

      ratio1_var=(ratio1^2)*(aoa1/dvr1^2+vart/dvrt^2-2*covt1/(dvr1*dvrt))
      ratio2_var=(ratio2^2)*(aoa2/dvr2^2+vart/dvrt^2-2*covt2/(dvr2*dvrt))

      chisq2=((ratio2-exp1)^2)/ratio2_var
      p3=pchisq(chisq2,1,lower.tail=F)

      #95% CI
      uci=ratio2+1.96*ratio2_var^.5
      lci=ratio2-1.96*ratio2_var^.5

      #z=list(var1=var1,var2=var2,var12=var12,mean_diff12_1=dvr1,mean_diff12_2=dvr2,var_diff1_2=aoa0,var_diff12_1=aoa1,var_diff12_2=aoa2,cov=cov,ratio1=ratio1,ratio2=ratio2,ratio1_var=ratio1_var,ratio2_var=ratio2_var,enrich_p=p3,upper_ratio1=uci,lower_ratio1=lci)
      z=list(var1=var1,var2=var2,var12=var12,mean_diff12_1=dvr1,mean_diff12_2=dvr2,var_diff1_2=aoa0,var_diff12_1=aoa1,var_diff12_2=aoa2,cov=cov,ratio=ratio2,ratio_var=ratio2_var,enrich_p=p3,upper_ratio=uci,lower_ratio=lci)
      return(z)


  }



