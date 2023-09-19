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
#' @param nv Sample size
#' @keywords R2 variance information matrix
#' 
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @return  This function will estimate significant difference between two PGS (either dependent or independent and joint or single). To get the test statistics for the difference between R2(y~x\[,v1]) and R2(y~x\[,v2]). (here we define R2_1=R2(y~x\[,v1])) and R2_2=R2(y~x\[,v2]))). The outputs are listed as follows.
#' \item{rsq1}{R2_1}
#' \item{rsq2}{R2_2}
#' \item{var1}{Variance of R2_1}
#' \item{var2}{variance of R2_2}
#' \item{var_diff}{Variance of difference between R2_1 and R2_2}
#' \item{r2_based_p}{two tailed P-value for significant difference between R2_1 and R2_2}
#' \item{r2_based_p_one_tail}{one tailed P-value for significant difference}
#' \item{mean_diff}{Differences between R2_1 and R2_2}
#' \item{upper_diff}{Upper limit of 95% CI for the difference}
#' \item{lower_diff}{Lower limit of 95% CI for the difference}
#' @examples
#' #To get the test statistics for the difference between R2(y~x[,1]) and 
#' #R2(y~x[,2]). (here we define R2_1=R2(y~x[,1])) and R2_2=R2(y~x[,2])))
#' 
#' dat=dat1
#' nv=length(dat$V1)
#' v1=c(1)
#' v2=c(2)
#' output=r2_diff(dat,v1,v2,nv)
#' output
#'
#' #r2redux output
#'
#' #output$rsq1 (R2_1)
#' #0.03836254
#' 
#' #output$rsq2 (R2_2)
#' #0.03881135
#' 
#' #output$var1 (variance of R2_1)
#' #0.0001436128
#' 
#' #output$var2 (variance of R2_2)
#' #0.0001451358
#' 
#' #output$var_diff (variance of difference between R2_1 and R2_2)
#' #5.678517e-07
#' 
#' #output$r2_based_p (two tailed p-value for significant difference)
#' #0.5514562
#'
#' #output$r2_based_p_one_tail(one tailed p-value for significant difference)
#' #0.2757281
#' 
#' #output$mean_diff (differences between R2_1 and R2_2)
#' #-0.0004488044
#' 
#' #output$upper_diff (upper limit of 95% CI for the difference)
#' #0.001028172
#' 
#' #output$lower_diff (lower limit of 95% CI for the difference)
#' #-0.001925781
#' 
#' #output$$p$nested
#' #1
#'
#' #output$$p$nonnested
#' #0.5514562
#'
#' #output$$p$LRT
#' #1
#'
#' #To get the test statistics for the difference between R2(y~x[,1]+x[,2]) and 
#' #R2(y~x[,2]). (here R2_1=R2(y~x[,1]+x[,2]) and R2_2=R2(y~x[,1]))
#' 
#' dat=dat1
#' nv=length(dat$V1)
#' v1=c(1,2)
#' v2=c(1)
#' output=r2_diff(dat,v1,v2,nv)
#'
#' #r2redux output
#'
#' #output$rsq1 (R2_1)
#' #0.03896678
#' 
#' #output$rsq2 (R2_2)
#' #0.03836254
#'
#' #output$var1 (variance of R2_1)
#' #0.0001473686
#' 
#' #output$var2 (variance of R2_2)
#' #0.0001436128
#' 
#' #output$var_diff (variance of difference between R2_1 and R2_2)
#' #2.321425e-06
#' 
#' #output$r2_based_p (p-value for significant difference between R2_1 and R2_2)
#' #0.4366883
#' 
#' #output$mean_diff (differences between R2_1 and R2_2)
#' #0.0006042383
#' 
#' #output$upper_diff (upper limit of 95% CI for the difference)
#' #0.00488788
#' 
#' #output$lower_diff (lower limit of 95% CI for the difference)
#' #-0.0005576171
#'
#'
#' #When faced with multiple predictors common between two models, for example, 
#' #y = any_cov1 + any_cov2 + ... + any_covN  + e   vs.
#' #y = PRS + any_cov1 + any_cov2 +...+ any_covN + e
#'
#' #A more streamlined approach can be adopted by consolidating the various
#' #predictors into a single predictor (see R code below).
#'
#' #R
#' #dat=dat1
#' #here let's assume, we wanted to test one PRS (dat$V2) 
#' #with 5 covariates (dat$V7 to dat$V11)
#' #mod1 <- lm(dat$V1~dat$V2 + dat$V7+ dat$V8+ dat$V9+ dat$V10+ dat$V11)
#' #merged_predictor1 <- mod1$fitted.values
#' #mod2 <- lm(dat$V1~ dat$V7+ dat$V8+ dat$V9+ dat$V10+ dat$V11)
#' #merged_predictor2 <- mod2$fitted.values
#' #dat=data.frame(dat$V1,merged_predictor1,merged_predictor2)
#'
#' #the comparison can be equivalently expressed as:
#' #y = merged_predictor1 + e   vs. 
#' #y = merged_predictor2 + e
#'   
#' #This comparison can be simply achieved using the r2_diff function, e.g.
#' 
#' #To get the test statistics for the difference between R2(y~x[,1]) and
#' #R2(y~x[,2]). (here x[,1]= merged_predictor2 (from full model), 
#' #and x[,2]= merged_predictor1(from reduced model))
#' #v1=c(1)
#' #v2=c(2)
#' #output=r2_diff(dat,v1,v2,nv)
#' #note that the merged predictor from the full model (v1) should be the first.
#' 
#' #str(output)
#' #List of 11
#' #$ rsq1               : num 0.0428
#' #$ rsq2               : num 0.042
#' #$ var1               : num 0.0.000158
#' #$ var2               : num 0.0.000156
#' #$ var_diff           : num 2.87e-06
#' #$ r2_based_p         : num 0.658
#' #$ r2_based_p_one_tail: num 0.329
#' #$ mean_diff          : num 0.000751
#' #$ upper_diff         : num 0.00407
#' #$ lower_diff         : num -0.00257
#' #$ p                  :List of 3
#' #..$ nested   : num 0.386
#' #..$ nonnested: num 0.658
#' #..$ LRT      : num 0.376
#' 
#' 
#' #Importantly note that in this case, merged_predictor1 is nested within
#' #merged_predictor2 (see mod1 vs. mod2 above). Therefore, this is 
#' #nested model comparison. So, output$p$nested (0.386) should be used 
#' #instead of output$p$nonnested (0.658). 
#' #Note that r2_based_p is the same as output$p$nonnested (0.658) here. 
#'   
#'
#' ##For this scenario, alternatively, the outcome variable (y) can be preadjusted
#' #with covariate(s), following the procedure in R:
#'
#' #mod <- lm(y ~ any_cov1 + any_cov2 + ... + any_covN)
#' #y_adj=scale(mod$residuals)
#' #then, the comparative significance test can be approximated by using
#' #the following model y_adj = PRS (r2_var(dat, v1, nv)) 
#' 
#' #R
#' #dat=dat1
#' #mod <- lm(dat$V1~dat$V7+ dat$V8+ dat$V9+ dat$V10+ dat$V11)
#' #y_adj=scale(mod$residuals)
#' #dat=data.frame(y_adj,dat$V2)
#' #v1=c(1)
#' #output=r2_var(dat, v1, nv)
#'
#' #str(output)
#' #$ var       : num 2e-06
#' #$ LRT_p     :Class 'logLik' : 0.98 (df=2)
#' #$ r2_based_p: num 0.977
#' #$ rsq       : num 8.21e-07
#' #$ upper_r2  : num 0.00403
#' #$ lower_r2  : num -0.000999
#'
#'
#' #In another scenario where the same covariates, but different
#' #PRS1 and PRS2 are compared,   
#' #y = PRS1 + any_cov1 + any_cov2 + ... + any_covN  + e   vs.
#' #y = PRS2 + any_cov1 + any_cov2 + ... + any_covN + e
#'      
#' #following approach can be employed (see R code below).
#'    
#' #R
#' #dat=dat1
#' #here let's assume dat$V2 as PRS1, dat$V3 as PRS2 and dat$V7 to dat$V11 as covariates
#' #mod1 <- lm(dat$V1~dat$V2 + dat$V7+ dat$V8+ dat$V9+ dat$V10+ dat$V11)
#' #merged_predictor1 <- mod1$fitted.values
#' #mod2 <- lm(dat$V1~dat$V3 + dat$V7+ dat$V8+ dat$V9+ dat$V10+ dat$V11)
#' #merged_predictor2 <- mod2$fitted.values
#' #dat=data.frame(dat$V1,merged_predictor2,merged_predictor1)
#'   
#' #the comparison can be equivalently expressed as:
#' #y = merged_predictor1 + e   vs. 
#' #y = merged_predictor2 + e
#' 
#' #This comparison can be simply achieved using the r2_diff function, e.g.
#' 
#' #To get the test statistics for the difference between R2(y~x[,1]) and
#' #R2(y~x[,2]). (here x[,1]= merged_predictor2, and x[,2]= merged_predictor1)
#' #v1=c(1)
#' #v2=c(2)
#' #output=r2_diff(dat,v1,v2,nv)
#'
#' #str(output)
#' #List of 11
#' #$ rsq1               : num 0.043
#' #$ rsq2               : num 0.0428
#' #$ var1               : num 0.000159
#' #$ var2               : num 0.000158
#' #$ var_diff           : num 2.6e-07
#' #$ r2_based_p         : num 0.657
#' #$ r2_based_p_one_tail: num 0.328
#' #$ mean_diff          : num 0.000227
#' #$ upper_diff         : num 0.00123
#' #$ lower_diff         : num 0.000773
#' #$ p                  :List of 3
#' #..$ nested   : num 0.634
#' #..$ nonnested: num 0.657
#' #..$ LRT      : num 0.627
#' 
#' #Importantly note that in this case, merged_predictor1 and merged_predictor2 
#' #are not nested to each other (see mod1 vs. mod2 above). 
#' #Therefore, this is nonnested model comparison. 
#' #So, output$p$nonnested (0.657) should be used instead of 
#' #output$p$nested (0.634). Note that r2_based_p is the same 
#' #as output$p$nonnested (0.657) here. 
#'
#'
#' #For the above non-nested scenario, alternatively, the outcome variable (y) 
#' #can be preadjusted with covariate(s), following the procedure in R:
#' #mod <- lm(y ~ any_cov1 + any_cov2 + ... + any_covN)
#' #y_adj=scale(mod$residuals)
#'
#' #R
#' #dat=dat1
#' #mod <- lm(dat$V1~dat$V7+ dat$V8+ dat$V9+ dat$V10+ dat$V11)
#' #y_adj=scale(mod$residuals)
#' #dat=data.frame(y_adj,dat$V3,dat$V2)
#'
#' #the comparison can be equivalently expressed as:
#' #y_adj = PRS1 + e   vs. 
#' #y_adj = PRS2 + e
#' #then, the comparative significance test can be approximated by using r2_diff function
#
#' #To get the test statistics for the difference between R2(y~x[,1]) and
#' #R2(y~x[,2]). (here x[,1]= PRS1 and x[,2]= PRS2)
#' #v1=c(1)
#' #v2=c(2)
#' #output=r2_diff(dat,v1,v2,nv)
#'
#' #str(output)
#' #List of 11
#' #$ rsq1               : num 5.16e-05
#' #$ rsq2               : num 4.63e-05
#' #$ var1               : num 2.21e-06
#' #$ var2               : num 2.18e-06
#' #$ var_diff           : num 1.31e-09
#' #$ r2_based_p         : num 0.884
#' #$ r2_based_p_one_tail: num 0.442
#' #$ mean_diff          : num 5.28e-06
#' #$ upper_diff         : num 7.63e-05
#' #$ lower_diff         : num -6.57e-05
#' #$ p                  :List of 3
#' #..$ nested   : num 0.942
#' #..$ nonnested: num 0.884
#' #..$ LRT      : num 0.942
#' 

  r2_diff = function (dat,v1,v2,nv) {
  
  dat=scale(dat);omat=cor(dat)
  
  
  if (length(v1)==1 & length(v2)==1) {
    ord=c(1,(1+v1),(1+v2))
    m1=lm(dat[,1]~dat[,(1+v1)])
    s1=summary(m1)
    m2=lm(dat[,1]~dat[,(1+v2)])
    s2=summary(m2)
    
    R2=s1$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv ) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var1=t100 
    
    
    R2=s2$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv ) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var2=t100 
    
    
    
    LR=-2*(logLik(m2)-logLik(m1))
    p1=pchisq(LR,1,lower.tail=F)
    
    dvr2=s1$r.squared-s2$r.squared
    #dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
    chi_dum=dvr2/(1/(nv)*(1-dvr2)^2) #NCP
    p2=pchisq(chi_dum,1,lower.tail=F)
    
    
    aoa=olkin1_2(omat[ord,ord],nv)
    chi_dum=dvr2^2/aoa
    p3=pchisq(chi_dum,1,lower.tail=F)
    uci=dvr2+1.96*aoa^.5
    lci=dvr2-1.96*aoa^.5
    
    

z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=as.numeric(aoa),r2_based_p=as.numeric(p3),r2_based_p_one_tail=as.numeric(p3/2),mean_diff=dvr2,upper_diff=as.numeric(uci),lower_diff=as.numeric(lci),p=list(nested=p2,nonnested=as.numeric(p3),LRT=as.numeric(p1)))
    #NOTE: r2_based_p=p3 due to normal distribution
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
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var1=t100 
    
    R2=s2$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv ) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var2=t100 
    
    LR=-2*(logLik(m2)-logLik(m1))
    p1=pchisq(LR,1,lower.tail=F)
    
    dvr=(s1$r.squared^.5-s2$r.squared^.5)^2 
    dvr2=s1$r.squared-s2$r.squared
    chi_dum=dvr2/(1/(nv )*(1-dvr2)^2) #NCP
    p2=pchisq(chi_dum,1,lower.tail=F)
    
    
    aoa=olkin12_1(omat[ord,ord],nv)
    #chi_dum2=dvr2^2/aoa
    #p3=pchisq(chi_dum2,1,lower.tail=F)
    
    
    #95% CI
    mv=1;lamda=chi_dum
    uci=qchisq(0.975,1,ncp=lamda)
    uci=(uci-lamda-mv)/(2*(mv+2*lamda))^.5
    uci=uci*aoa^.5+dvr2
    lci=qchisq(0.025,1,ncp=lamda)
    lci=(lci-lamda-mv)/(2*(mv+2*lamda))^.5
    lci=lci*aoa^.5+dvr2
    
    z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var=var2,var_diff=as.numeric(aoa),LRT_p=as.numeric(p1),r2_based_p=as.numeric(p2),r2_based_p_one_tail=as.numeric(p2/2),mean_diff=dvr2,upper_diff=as.numeric(uci),lower_diff=as.numeric(lci))
    #NOTE: r2_based_p=p2 due to chi^2 distribution
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
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var1=t100 
    
    R2=s2$r.squared;mv2=1 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var2=t100 
    
    
    
    #LR=-2*(logLik(m2)-logLik(m1))
    #p1=pchisq(LR,1,lower.tail=F)
    #dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
    dvr2=s1$r.squared-s2$r.squared
    #chi_dum=dvr2/(1/(nv)*(1-dvr2)^2) #NCP
    #p2=pchisq(chi_dum,1,lower.tail=F)
    
    aoa=olkin12_3(omat[ord,ord],nv) #check
    chi_dum2=dvr2^2/aoa
    p3=pchisq(chi_dum2,1,lower.tail=F)
    
    
    #95% CI
    uci=dvr2+1.96*aoa^.5
    lci=dvr2-1.96*aoa^.5
    
    
    z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=as.numeric(aoa),r2_based_p=as.numeric(p3),r2_based_p_one_tail=as.numeric(p3/2),mean_diff=dvr2,upper_diff=as.numeric(uci),lower_diff=as.numeric(lci),p=list(nested=p2,nonnested=as.numeric(p3),LRT=as.numeric(p1)))
    #NOTE: r2_based_p=p3 due to normal distribution
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
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var1=t100 
    
    
    R2=s2$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var2=t100 
    
    
    #LR=-2*(logLik(m2)-logLik(m1))
    #p1=pchisq(LR,1,lower.tail=F)
    #dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
    dvr2=s1$r.squared-s2$r.squared
    #chi_dum=dvr2/(1/(nv)*(1-dvr2)^2) #NCP
    #p2=pchisq(chi_dum,1,lower.tail=F)
    
    aoa=olkin12_13(omat[ord,ord],nv)
    chi_dum2=dvr2^2/aoa
    p3=pchisq(chi_dum2,1,lower.tail=F)
    
    
    #95% CI
    uci=dvr2+1.96*aoa^.5
    lci=dvr2-1.96*aoa^.5
    z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=as.numeric(aoa),r2_based_p=as.numeric(p3),r2_based_p_one_tail=as.numeric(p3/2),mean_diff=dvr2,upper_diff=as.numeric(uci),lower_diff=as.numeric(lci),p=list(nested=p2,nonnested=as.numeric(p3),LRT=as.numeric(p1)))
    #NOTE: r2_based_p=p3 due to normal distribution
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
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var1=t100 
    
    R2=s2$r.squared;mv2=2 #expected variance for s1r2
    t100=(1/(nv) *(1-R2)^2) #Infor matrix
    lamda=R2/t100 #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
    t100=t100^2*2*(mv2+2*lamda) #var(beta)^2*var(non-cental chi)
    var2=t100 
    
    
    #LR=-2*(logLik(m2)-logLik(m1))
    #p1=pchisq(LR,1,lower.tail=F)
    #dvr=(s1$r.squared^.5-s2$r.squared^.5)^2
    dvr2=s1$r.squared-s2$r.squared
    #chi_dum=dvr2/(1/(nv)*(1-dvr2)^2) #NCP
    #p2=pchisq(chi_dum,1,lower.tail=F)
    
    aoa=olkin12_34(omat[ord,ord],nv)
    chi_dum2=dvr2^2/aoa
    p3=pchisq(chi_dum2,1,lower.tail=F)
    
    #95% CI
    uci=dvr2+1.96*aoa^.5
    lci=dvr2-1.96*aoa^.5
    z=list(rsq1=s1$r.squared,rsq2=s2$r.squared,var1=var1,var2=var2,var_diff=as.numeric(aoa),r2_based_p=as.numeric(p3),r2_based_p_one_tail=as.numeric(p3/2),mean_diff=dvr2,upper_diff=as.numeric(uci),lower_diff=as.numeric(lci),p=list(nested=p2,nonnested=as.numeric(p3),LRT=as.numeric(p1)))
    #NOTE: r2_based_p=p3 due to normal distribution
    return(z)
    
  }
  

}





