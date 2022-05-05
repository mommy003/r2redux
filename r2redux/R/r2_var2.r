  #' r2_var function
  #'
  #' This function estimates var(R2(y~x\[,v1]))
  #' where R2 is the R squared value of the model,
  #' where R2 is the R squared value of the model,
  #' y is N by 1 matrix having the dependent variable, and
  #' x is N by M matrix having M explanatory variables.
  #' v1 indicates the ith column in the x matrix
  #' (v1 can be multiple values between 1 - M, see Arguments below)
  #' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
  #' @param v1 This can be set as v1=c(1), v1=c(1,2) or possibly with more values
  #' @param nv sample size
  #' @keywords R2 variance information matrix
  #' @export
  #' @examples 
  #' To get the test statistics for R2(y~x[,v1])
  #' 
  #' dat=read.table("test_ukbb_thresholds_scaled") (see example file)
  #' nv=length(dat$V1)
  #' v1=c(1) 
  #' output=r2_var(dat,v1,nv)
  #'
  #' r2redux output
  #'
  #' output$var (variance of R2) 
  #' 0.0001437583
  #' 
  #' output$r2_based_p (P-value under the null hypothesis, i.e. R2=0)
  #' 1.213645e-10
  #' 
  #' output$mean_r2 (R2)
  #' 0.03836254
  #' 
  #' output$upper_r2 (upper limit of 95% CI for R2)
  #' 0.06435214
  #' 
  #' output$lower_r2 (lower limit of 95% CI for R2)
  #' 0.01763347
  #'
  #'
  #'
  #' To get the test statistic for R2(y~x[,v1]+x[,v2]+x[,v3])
  #' 
  #' dat=read.table("test_ukbb_thresholds_scaled") (see example file)
  #' nv=length(dat$V1)
  #' v1=c(1,2,3) 
  #' outout=r2_var(dat,v1,nv)
  #' output
  #'
  #' r2redux output
  #'
  #' output$var (variance of R2)
  #' 0.0001499374
  #' 
  #' output$r2_based_p (R2 based P-value)
  #' 7.461267e-11
  #' 
  #' output$mean_r2 (R2)
  #' 0.03917668
  #' 
  #' output$upper_r2 (upper limit of 95% CI for R2)
  #' 0.06538839
  #' 
  #' output$lower_r2 (lower limit of 95% CI for R2)
  #' 0.01821657

  r2_var = function (dat,v1,nv) {
    #source("aoa12_1.r")
    #source("aoa12_13.r")
    #source("aoa1_2.r")
    #source("aoa12_3.r")
    #source("aoa12_34.r")

    dat=scale(dat);omat=cor(dat)

      ord=c(1,1+v1)

      m0=lm(dat[,1]~1)
      s0=summary(m0)
      m1=lm(dat[,1]~as.matrix(dat[,(1+v1)]))
      s1=summary(m1)

      R2=s1$r.squared;mv2=length(v1)    #expected variance for s1r2
      t100=(1/(nv-1) *(1-R2)^2) #Infor matrix
      lamda=R2/t100           #non-central chi^2 (mu=k+lamda, var=2*(k+2*lamda))
      t100=t100^2*2*(mv2+2*lamda)  #var(beta)^2*var(non-cental chi)
      var1=t100  #*(1-R2)

      LR=-2*(logLik(m0)-logLik(m1))
      p1=pchisq(LR,mv2,lower.tail=F)

      dvr2=R2
      chi_dum=dvr2/(1/(nv-1)*(1-dvr2)^2) #NCP
      p2=pchisq(chi_dum,1,lower.tail=F)

      #95% CI
      mv=mv2;lamda=chi_dum
      uci=qchisq(0.975,1,ncp=lamda)
      uci=(uci-lamda-1)/(2*(mv+2*lamda))^.5
      uci=uci*var1^.5+dvr2
      lci=qchisq(0.025,1,ncp=lamda)
      lci=(lci-lamda-1)/(2*(mv+2*lamda))^.5
      lci=lci*var1^.5+dvr2

      z=list(var=var1,LRT_p=p1,r2_based_p=p2,mean_r2=dvr2,upper_r2=uci,lower_r2=lci)
      #NOTE: r2_based_p=p2 due to chi^2 distribution
      return(z)

  }
