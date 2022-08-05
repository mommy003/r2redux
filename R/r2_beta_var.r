  
#' r2_beta_var 
#'
#' This function estimates var(beta1^2) and (beta2^2), and
#' beta1 and 2 are regression coefficients from a multiple regression model,
#' i.e. y = x1 * beta1 + x2 * beta2 +e, y, x1 and x2 are column-standardised
#' (see Olkin and Finn 1995).
#' y is N by 1 matrix having the dependent variable, 
#' x1 is N by 1 matrix having the ith explanatory variable.
#' x2 is N by 1 matrix having the jth explanatory variable.
#' v1 and v2 indicates the ith and jth column in the data
#' (v1 or v2 should be a single interger between 1 - M, see Arguments below).
#' @references
#' Olkin, I. and Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
#' @param dat N by (M+1) matrix having variables in the order of cbind(y,x)
#' @param v1 This can be set as v1=1, v1=2, v1=3  or any value between 1 - M based on combination
#' @param v2 This can be set as v2=1, v2=2, v2=3, or any value between 1 - M based on combination 
#' @param nv Sample size
#' @keywords variance of beta^2 from a multiple regression
#' @export
#' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
#' @return  This function will estiamte the variance of beta1^2 and beta2^2, and the covariance between beta1^2 and beta2^2, i.e. the information matrix of squared regression coefficients. beta1 and beta2 are regression coefficients from a multiple regression model, i.e. y = x1 * beta1 + x2 * beta2 +e, where y, x1 and x2 are column-standardised. The outputs are listed as follows.
#' \item{beta1_sq}{beta1_sq}
#' \item{beta2_sq}{beta2_sq}
#' \item{var1}{Variance of beta1_sq}
#' \item{var2}{Variance of beta2_sq}
#' \item{var1_2}{Variance of difference between beta1_sq and beta2_sq}
#' \item{cov}{Covariance between beta1_sq and beta2_sq}
#' \item{upper_beta1_sq}{upper limit of 95% CI for beta1_sq}
#' \item{lower_beta1_sq}{lower limit of 95% CI for beta1_sq}
#' \item{upper_beta2_sq}{upper limit of 95% CI for beta2_sq}
#' \item{lower_beta2_sq}{lower limit of 95% CI for beta2_sq}
#' @examples
#' #To get the 95% CI of beta1_sq and beta2_sq
#' #beta1 and beta2 are regression coefficients from a multiple regression model,
#' #i.e. y = x1 * beta1 + x2 * beta2 +e, where y, x1 and x2 are column-standardised.
#'
#' dat=dat2
#' nv=length(dat$V1)
#' v1=c(1)
#' v2=c(2)
#' output=r2_beta_var(dat,v1,v2,nv)
#' output
#' #r2redux output
#' #output$beta1_sq (beta1_sq)
#' #0.01118301
#' 
#' #output$beta2_sq (beta2_sq)
#' #0.004980285
#' 
#' #output$var1 (variance of beta1_sq)
#' #7.072931e-05
#' 
#' #output$var2 (variance of beta2_sq)
#' #3.161929e-05
#' 
#' #output$var1_2 (variance of difference between beta1_sq and beta2_sq)
#' #0.000162113
#' 
#' #output$cov (covariance between beta1_sq and beta2_sq)
#' #-2.988221e-05
#' 
#' #output$upper_beta1_sq (upper limit of 95% CI for beta1_sq)
#' #0.03037793
#' 
#' #output$lower_beta1_sq (lower limit of 95% CI for beta1_sq)
#' #-0.00123582
#' 
#' #output$upper_beta2_sq (upper limit of 95% CI for beta2_sq)
#' #0.02490076
#' 
#' #output$lower_beta2_sq (lower limit of 95% CI for beta2_sq)
#' #-0.005127546


r2_beta_var = function (dat,v1,v2,nv) {


    dat=scale(dat);omat=cor(dat)
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


      #95% CI
      mv=2;lamda=dvr1/(1/(nv)*(1-dvr1)^2)
      uci1=qchisq(0.975,1,ncp=lamda)
      uci1=(uci1-lamda-1)/(2*(mv+2*lamda))^.5
      uci1=uci1*var1^.5+dvr1
      lci1=qchisq(0.025,1,ncp=lamda)
      lci1=(lci1-lamda-1)/(2*(mv+2*lamda))^.5
      lci1=lci1*var1^.5+dvr1

      mv=2;lamda=dvr2/(1/(nv)*(1-dvr2)^2)
      uci2=qchisq(0.975,1,ncp=lamda)
      uci2=(uci2-lamda-1)/(2*(mv+2*lamda))^.5
      uci2=uci2*var1^.5+dvr2
      lci2=qchisq(0.025,1,ncp=lamda)
      lci2=(lci2-lamda-1)/(2*(mv+2*lamda))^.5
      lci2=lci2*var1^.5+dvr2


      z=list(beta1_sq=dvr1,beta2_sq=dvr2,var1=var1,var2=var2,var1_2=var1_2,cov=cov,upper_beta1_sq=uci1,lower_beta1_sq=lci1,upper_beta2_sq=uci2,lower_beta2_sq=lci2)

      return(z)

  }