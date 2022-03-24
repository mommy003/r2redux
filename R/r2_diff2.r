  #' r2_diff function
  #'
  #' This function estimates var(R2(y~x\[,v1]) - R2(y~x\[,v2]))
  #' where R2 is the R squared value of the model,
  #' y is N by 1 matrix having the dependent variable, and
  #' x is N by M matrix having M explanatory variables.
  #' v1 or v2 indicates the ith column in the x matrix 
  #' (v1 or v2 can be multiple values between 1 - M, see Arguments below)
  #' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
  #' @param v1 This can be set as v1=c(1) or v1=c(1,2)
  #' @param v2 This can be set as v2=c(2), v2=c(3), v2=c(1,3) or v2=c(3,4)
  #' @param nv sample size
  #' @keywords R2 variance information matrix
  #' @export
  #' @examples
  #' To get the test statistics for the difference between R2(y~x[,v1]) and 
  #' R2(y~x[,v2]). (here we define R2_1=R2(y~x[,v1])) and R2_2=R2(y~x[,v2])))
  #'
  #' dat=read.table("test_ukbb_thresholds_scaled") (see example files)
  #' nv=length(dat$V1)
  #' v1=c(1)
  #' v2=c(2)
  #' output=r2_diff(dat,v1,v2,nv)
  #' output
  #'
  #' r2redux output
  #' 
  #' output$var1 (variance of R2_1)
  #' 0.0001437583
  #' 
  #' output$var2 (variance of R2_2)
  #' 0.0001452828
  #' 
  #' output$var_diff (variance of difference between R2_1 and R2_2)
  #' 5.678517e-07
  #' 
  #' output$r2_based_p (p-value for significant difference between R2_1 and R2_2)
  #' 0.5514562
  #' 
  #' output$mean_diff (differences between R2_1 and R2_2)
  #' -0.0004488044
  #' 
  #' output$upper_diff (upper limit of 95% CI for the difference)
  #' 0.001028172
  #' 
  #' output$lower_diff (lower limit of 95% CI for the difference)
  #' -0.001925781
  #'
  #'
  #'
  #' To get the test statistics for the difference between R2(y~x[,v1]+x[,v2]) and 
  #' R2(y~x[,v2]). (here R2_1=R2(y~x[,v1]+x[,v2]) and R2_2=R2(y~x[,v1]))
  #' 
  #' dat=read.table("test_ukbb_thresholds_scaled") (see example files)
  #' nv=length(dat$V1)
  #' v1=c(1,2)
  #' v2=c(1)
  #' output=r2_diff(dat,v1,v2,nv)
  #' output
  #'
  #' r2redux output
  #' 
  #' output$var1 (variance of R2_1)
  #' 0.0001475195
  #' 
  #' output$var2 (variance of R2_2)
  #' 0.0001437583
  #' 
  #' output$var_diff (variance of difference between R2_1 and R2_2)
  #' 2.321425e-06
  #' 
  #' output$r2_based_p (p-value for significant difference between R2_1 and R2_2)
  #' 0.4369177
  #' 
  #' output$mean_diff (differences between R2_1 and R2_2)
  #' 0.0006042383
  #' 
  #' output$upper_diff (upper limit of 95% CI for the difference)
  #' 0.004887989
  #' 
  #' output$lower_diff (lower limit of 95% CI for the difference)
  #' -0.0005574975
  #'
  #' 


  r2_diff = function (dat,v1,v2,nv) {
    source("aoa12_1.r")
    source("aoa12_13.r")
    source("aoa1_2.r")
    source("aoa12_3.r")
    source("aoa12_34.r")

    dat=scale(dat);omat=cor(dat)

    if (length(v1)==1 & length(v2)==1) {
      ord=c(1,(1+v1),(1+v2))
      #omat=omat[ord,ord]

      m1=lm(dat[,1]~dat[,(1+v1)])
      s1=summary(m1)

      m2=lm(dat[,1]~dat[,(1+v2)])
      s2=summary(m2)

      R2=s1$r.squared;mv2=1    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100 #*(1-R2)

      R2=s2$r.squared;mv2=1    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var2=t100   #*(1-R2)

      LR=-2*(logLik(m2)-logLik(m1))
      p1=pchisq(LR,1,lower.tail=F)

      dvr2=s1$r.squared-s2$r.squared
      dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      aoa=olkin1_2(omat[ord,ord],nv)
      chi_dum=dvr2^2/aoa
      p3=pchisq(chi_dum,1,lower.tail=F)

      uci=dvr2+1.96*aoa^.5
      lci=dvr2-1.96*aoa^.5

      z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=aoa,LRT_p=p1,p2=p2,r2_based_p=p3,mean_diff=dvr2,upper_diff=uci,lower_diff=lci)
      #NOTE: r2_based_p=p3 due to normal distribution
      return(z)
    }

    if (length(v1)==2 & length(v2)==1 & length(unique(c(v1,v2)))==2) {
      if (v1[1]==v2[1]) {ord=c(1,(1+v1[1]),(1+v1[2]))}
      if (v1[2]==v2[1]) {ord=c(1,(1+v1[2]),(1+v1[1]))}
      #omat=omat[ord,ord]

      m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
      s1=summary(m1)

      m2=lm(dat[,1]~dat[,(1+v2)])
      s2=summary(m2)

      R2=s1$r.squared;mv2=2    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100  #*(1-R2)

      R2=s2$r.squared;mv2=1    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var2=t100   #*(1-R2)

      LR=-2*(logLik(m2)-logLik(m1))
      p1=pchisq(LR,1,lower.tail=F)

      dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
      dvr2=s1$r.squared-s2$r.squared
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      aoa=olkin12_1(omat[ord,ord],nv)
      chi_dum2=dvr2^2/aoa
      p3=pchisq(chi_dum2,1,lower.tail=F)

      #95% CI
      mv=1;lamda=chi_dum
      uci=qchisq(0.975,1,ncp=lamda)
      uci=(uci-lamda-mv)/(2*(mv+2*lamda))^.5
      uci=uci*aoa^.5+dvr2
      lci=qchisq(0.025,1,ncp=lamda)
      lci=(lci-lamda-mv)/(2*(mv+2*lamda))^.5
      lci=lci*aoa^.5+dvr2

      z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=aoa,LRT_p=p1,p3=p3,r2_based_p=p2,mean_diff=dvr2,upper_diff=uci,lower_diff=lci)
      #NOTE: r2_based_p=p2 due to chi^2 distribution
      return(z)
    }

    if (length(v1)==2 & length(v2)==1 & length(unique(c(v1,v2)))==3) {
      ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]))
      #omat=omat[ord,ord]

      m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
      s1=summary(m1)

      m2=lm(dat[,1]~dat[,(1+v2)])
      s2=summary(m2)

      R2=s1$r.squared;mv2=2    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100  #*(1-R2)

      R2=s2$r.squared;mv2=1    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var2=t100  #*(1-R2)

      LR=-2*(logLik(m2)-logLik(m1))
      p1=pchisq(LR,1,lower.tail=F)

      dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
      dvr2=s1$r.squared-s2$r.squared
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      aoa=olkin12_3(omat[ord,ord],nv)    #check
      chi_dum2=dvr2^2/aoa
      p3=pchisq(chi_dum2,1,lower.tail=F)

      #95% CI
      uci=dvr2+1.96*aoa^.5
      lci=dvr2-1.96*aoa^.5

      z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=aoa,LRT_p=p1,p2=p2,r2_based_p=p3,mean_diff=dvr2,upper_diff=uci,lower_diff=lci)
      #NOTE: r2_based_p=p3 due to normal distribution
      return(z)
    }


    if (length(v1)==2 & length(v2)==2 & length(unique(c(v1,v2)))==3) {
      if (v1[1]==v2[1]) {ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[2]))}
      if (v1[1]==v2[2]) {ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]))}
      if (v1[2]==v2[1]) {ord=c(1,(1+v1[2]),(1+v1[1]),(1+v2[2]))}
      if (v1[2]==v2[2]) {ord=c(1,(1+v1[2]),(1+v1[1]),(1+v2[1]))}
      #omat=omat[ord,ord]

      m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
      s1=summary(m1)

      m2=lm(dat[,1]~dat[,1+v2[1]]+dat[,1+v2[2]])
      s2=summary(m2)

      R2=s1$r.squared;mv2=2    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100  #*(1-R2)

      R2=s2$r.squared;mv2=2    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var2=t100  #*(1-R2)

      LR=-2*(logLik(m2)-logLik(m1))
      p1=pchisq(LR,1,lower.tail=F)

      dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
      dvr2=s1$r.squared-s2$r.squared
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      aoa=olkin12_13(omat[ord,ord],nv)
      chi_dum2=dvr2^2/aoa
      p3=pchisq(chi_dum2,1,lower.tail=F)

      #95% CI
      uci=dvr2+1.96*aoa^.5
      lci=dvr2-1.96*aoa^.5

      z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=aoa,LRT_p=p1,p2=p2,r2_based_p=p3,mean_diff=dvr2,upper_diff=uci,lower_diff=lci)
      #NOTE: r2_based_p=p3 due to normal distribution
      return(z)
    }

    if (length(v1)==2 & length(v2)==2 & length(unique(c(v1,v2)))==4) {
      ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]),(1+v2[2]))
      #omat=omat[ord,ord]

      m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
      s1=summary(m1)

      m2=lm(dat[,1]~dat[,1+v2[1]]+dat[,1+v2[2]])
      s2=summary(m2)

      R2=s1$r.squared;mv2=2    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100  #*(1-R2)

      R2=s2$r.squared;mv2=2    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var2=t100  #*(1-R2)

      LR=-2*(logLik(m2)-logLik(m1))
      p1=pchisq(LR,1,lower.tail=F)

      dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
      dvr2=s1$r.squared-s2$r.squared
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      aoa=olkin12_34(omat[ord,ord],nv)
      chi_dum2=dvr2^2/aoa
      p3=pchisq(chi_dum2,1,lower.tail=F)

      #95% CI
      uci=dvr2+1.96*aoa^.5
      lci=dvr2-1.96*aoa^.5

      z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=aoa,LRT_p=p1,p2=p2,r2_based_p=p3,mean_diff=dvr2,upper_diff=uci,lower_diff=lci)
      #NOTE: r2_based_p=p3 due to normal distribution
      return(z)
    }



  }


