
#' r_diff function
#'
#' This function estimates var(R(y~x\[,v1]) - R(y~x\[,v2]))
#' where R is the correlation between y and x,
#' y is N by 1 matrix having the dependent variable, and
#' x is N by M matrix having M explanatory variables.
#' v1 or v2 indicates the ith column in the x matrix 
#' (v1 or v2 can be multiple values between 1 - M, see Arguments below)
#' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
#' @param v1 This can be set as v1=c(1) or v1=c(1,2)
#' @param v2 This can be set as v2=c(2), v2=c(3), v2=c(1,3) or v2=c(3,4)
#' @param nv Sample size
#' @keywords R2 variance information matrix
#' 
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @return  This function will estimate significant difference between two PGS (either dependent or independent and joint or single). To get the test statistics for the difference between R2(y~x\[,v1]) and R2(y~x\[,v2]). (here we define R2_1=R2(y~x\[,v1])) and R2_2=R2(y~x\[,v2]))). The outputs are listed as follows.
#' \item{r1}{R2_1}
#' \item{r2}{R2_2}
#' \item{var1}{Variance of R_1}
#' \item{var2}{variance of R_2}
#' \item{var_diff}{Variance of difference between R_1 and R_2}
#' \item{r2_based_p}{P-value for significant difference between R_1 and R_2 for two tailed test}
#' \item{r_based_p_one_tail}{P-value for significant difference between R_1 and R_2 for one tailed test}
#' \item{mean_diff}{Differences between R_1 and R_2}
#' \item{upper_diff}{Upper limit of 95% CI for the difference}
#' \item{lower_diff}{Lower limit of 95% CI for the difference}
#' @examples
#' #To get the test statistics for the difference between R(y~x[,v1]) and 
#' #R(y~x[,v2]). (here we define R_1=R(y~x[,v1])) and R_2=R(y~x[,v2])))
#' 
#' dat=dat1
#' nv=length(dat$V1)
#' v1=c(1)
#' v2=c(2)
#' output=r_diff(dat,v1,v2,nv)
#' output
#' 
#' #r2redux output
#' 
#' #output$r1 (R_1)
#' #0.1958636
#' 
#' #output$r2 (R-2)
#' #0.197006
#' 
#' #output$var1 (variance of R_1)
#' #0.0009247466
#' 
#' #output$var2 (variance of R_1)
#' #0.0001451358
#' 
#' #output$var_diff (variance of difference between R_1 and R_2)
#' #3.65286e-06
#' 
#' #output$r_based_p (two tailed p-value for significant difference between R2_1 and R2_2)
#' #0.5500319
#' 
#' #output$r_based_p_one_tail (one tailed p-value for significant difference between R2_1 and R2_2)
#' #0.2750159
#' 
#' #output$mean_diff
#' #-0.001142375 (differences between R2_1 and R2_2)
#' 
#' #output$upper_diff (upper limit of 95% CI for the difference)
#' #0.002603666
#' 
#' #output$lower_diff (lower limit of 95% CI for the difference)
#' #-0.004888417
#' 
#' 
#' #To get the test statistics for the difference between R(y~x[,v1]+[,v2]) and 
#' #R(y~x[,v2]). (here R_1=R(y~x[,v1]+x[,v2]) and R_2=R(y~x[,v1]))
#' 
#' nv=length(dat$V1)
#' v1=c(1,2)
#' v2=c(2)
#' output=r_diff(dat,v1,v2,nv)
#' output
#' 
#' #output$rsq1
#' #0.03896678
#'
#' #output$rsq2
#' #0.03881135
#'
#' #output$var1
#' #0.0001473686
#'
#' #output$var2
#' #0.0001451358
#'
#' #output$var_diff
#' #5.974193e-07
#'
#' #output$r2_based_p
#' #0.6933505
#'
#' #output$mean_diff
#' #0.0001554339
#'
#' #output$upper_diff
#' #0.002355053
#'
#' #output$lower_diff
#' #-0.0003955781



r_diff = function (dat,v1,v2,nv) {
  
  
  dat=scale(dat);omat=cor(dat)
  
if (length(v1)==1 & length(v2)==1) {
    ord=c(1,(1+v1),(1+v2))
    m1=lm(dat[,1]~dat[,(1+v1)])
    s1=summary(m1)
    m2=lm(dat[,1]~dat[,(1+v2)])
    s2=summary(m2)
        
    R2=s1$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv ) *(1-R2)^2) #Infor matrix
    var1=t100
    R2=s2$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv ) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var2=t100 
    dvr=(s1$r.squared^.5-s2$r.squared^.5)
    #aoa=olkin1_2(omat[ord,ord],nv)
    r_aoa=r_olkin1_2(omat[ord,ord],nv)
    #chi_dum=dvr2^2/aoa
    r_chi_dum=dvr^2/r_aoa
    #p3=pchisq(chi_dum,1,lower.tail=F)
    r_p3=pchisq(r_chi_dum,1,lower.tail=F )

    #uci=dvr2+1.96*aoa^.5
    #lci=dvr2-1.96*aoa^.5
    r_uci=dvr+1.96*r_aoa^.5
    r_lci=dvr-1.96*r_aoa^.5
    z=list(r1=s1$r.squared^.5,r2=s2$r.squared^.5,var1=var1,var2=var2,var_diff=r_aoa,r_based_p=r_p3,r_based_p_one_tail=r_p3/2 ,mean_diff=dvr,upper_diff=r_uci,lower_diff=r_lci)
    #NOTE: r_based_p=r_p3 due to normal distribution
    return(z)
    
  }
  
  
  if (length(v1)==2 & length(v2)==1 & length(unique(c(v1,v2)))==2) {
    if (v1[1]==v2[1]) {ord=c(1,(1+v1[1]),(1+v1[2]))}
    if (v1[2]==v2[1]) {ord=c(1,(1+v1[2]),(1+v1[1]))}
    
    
    m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
    s1=summary(m1)
    m2=lm(dat[,1]~dat[,(1+v2)])
    s2=summary(m2)
    
    
    
    R2=s1$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var1=t100
    R2=s2$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv ) *(1-R2)^2) #Infor matrix
    var2=t100 
    LR=-2*(logLik(m2)-logLik(m1))
    p1=pchisq(LR,1,lower.tail=F)
    dvr=(s1$r.squared^.5-s2$r.squared^.5)
    
    #chi_dum=dvr2/(1/(nv )*(1-dvr2)^2) #NCP
    #p2=pchisq(chi_dum,1,lower.tail=F)
    r_aoa=r_olkin12_1(omat[ord,ord],nv)
    r_chi_dum=dvr^2/r_aoa
    r_p3=pchisq(r_chi_dum,1,lower.tail=F)
    
    
    #95% CI
    r_uci=dvr+1.96*r_aoa^.5
    r_lci=dvr-1.96*r_aoa^.5
    
    z=list(r1=s1$r.squared^.5,r2=s2$r.squared^.5,var1=var1,var2=var2,var_diff=r_aoa,LRT_p=p1,r_based_p=r_p3,mean_diff=dvr,upper_diff=r_uci,lower_diff=r_lci)
    #NOTE: r_based_p=r_p3 due to normal  
    return(z)
    
  }
  
  
  
  if (length(v1)==2 & length(v2)==1 & length(unique(c(v1,v2)))==3) {
    ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]))

    m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
    s1=summary(m1)
    m2=lm(dat[,1]~dat[,(1+v2)])
    s2=summary(m2)
    
    
    
    R2=s1$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var1=t100 
    R2=s2$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var2=t100 
    
    
    dvr=(s1$r.squared^.5-s2$r.squared^.5)
    r_aoa=r_olkin12_3(omat[ord,ord],nv) #check
    r_chi_dum=dvr^2/r_aoa
    r_p3=pchisq(r_chi_dum,1,lower.tail=F)
    
    
    
    #95% CI
    r_uci=dvr+1.96*r_aoa^.5
    r_lci=dvr-1.96*r_aoa^.5
    
    z=list(r1=s1$r.squared^.5,r2=s2$r.squared^.5,var1=var1,var2=var2,var_diff=r_aoa,r_based_p=r_p3,r_based_p_one_tail=r_p3/2,mean_diff=dvr,upper_diff=r_uci,lower_diff=r_lci)
    #NOTE: r_based_p=r_p3 due to normal distribution
    return(z)
    
  }
  
  
  if (length(v1)==2 & length(v2)==2 & length(unique(c(v1,v2)))==3) {
    if (v1[1]==v2[1]) {ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[2]))}
    if (v1[1]==v2[2]) {ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]))}
    if (v1[2]==v2[1]) {ord=c(1,(1+v1[2]),(1+v1[1]),(1+v2[2]))}
    if (v1[2]==v2[2]) {ord=c(1,(1+v1[2]),(1+v1[1]),(1+v2[1]))}
    
    m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
    s1=summary(m1)
    m2=lm(dat[,1]~dat[,1+v2[1]]+dat[,1+v2[2]])
    s2=summary(m2)
    
    R2=s1$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var1=t100 
    R2=s2$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var2=t100 
    dvr=(s1$r.squared^.5-s2$r.squared^.5)
    r_aoa=r_olkin12_13(omat[ord,ord],nv)
    r_chi_dum=dvr^2/r_aoa
    r_p3=pchisq(r_chi_dum,1,lower.tail=F)
    
    
    #95% CI
    r_uci=dvr+1.96*r_aoa^.5
    r_lci=dvr-1.96*r_aoa^.5
    
    z=list(r1=s1$r.squared^.5,r2=s2$r.squared^.5,var1=var1,var2=var2,var_diff=r_aoa,r_based_p=r_p3, r_based_p_one_tail=r_p3/2,mean_diff=dvr,upper_diff=r_uci,lower_diff=r_lci)
    #NOTE: r_based_p=r_p3 due to normal distribution
    return(z)
    
  }
  
  
  if (length(v1)==2 & length(v2)==2 & length(unique(c(v1,v2)))==4) {
      ord=c(1,(1+v1[1]),(1+v1[2]),(1+v2[1]),(1+v2[2]))
    
    
    m1=lm(dat[,1]~dat[,1+v1[1]]+dat[,1+v1[2]])
    s1=summary(m1)
    m2=lm(dat[,1]~dat[,1+v2[1]]+dat[,1+v2[2]])
    s2=summary(m2)
    
    R2=s1$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var1=t100 
    R2=s2$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    var2=t100 
    dvr=(s1$r.squared^.5-s2$r.squared^.5)
    r_aoa=r_olkin12_34(omat[ord,ord],nv)
    r_chi_dum=dvr^2/r_aoa
    r_p3=pchisq(r_chi_dum,1,lower.tail=F)
    
    
    
    #95% CI
    r_uci=dvr+1.96*r_aoa^.5
    r_lci=dvr-1.96*r_aoa^.5
    
    z=list(r1=s1$r.squared^.5,r2=s2$r.squared^.5,var1=var1,var2=var2,var_diff=r_aoa,r_based_p=r_p3, r_based_p_one_tail=r_p3/2,mean_diff=dvr,upper_diff=r_uci,lower_diff=r_lci)
    #NOTE: r_based_p=p3 due to normal distribution
    return(z)
    
  }
  
}

