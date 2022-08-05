  #' olkin_beta_ratio function
  #'
  #' This function derives variance of beta1^2 / R^2 
  #' where beta1 and 2 are regression coefficients from a multiple regression model,
  #' i.e. y = x1 * beta1 + x2 * beta2 + e, where y, x1 and x2 are column-standardised
  #' (see Olkin and Finn 1995).
  #' @references
  #' Olkin, I. and Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
  #' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
  #' @param nv sampel size
  #' @keywords information matrix in the context of correlation
  #' @export
  #' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
  #' @return  This function will generate the variance of the proportion, i.e. beta1_2/R^2.The outputs are listed as follows.
  #' \item{ratio_var}{Variance of ratio}
  #' @examples
  #' #To get information (variance-covariance) matrix of beta1 and beta2 where 
  #' #beta1 and 2 are regression coefficients from a multiple regression model.
  #' dat=dat2
  #' omat=cor(dat)[1:3,1:3]
  #' #omat
  #' #1.0000000 0.1497007 0.136431
  #' #0.1497007 1.0000000 0.622790
  #' #0.1364310 0.6227900 1.000000
  #' 
  #' nv=length(dat$V1)
  #' output=olkin_beta_ratio(omat,nv)
  #' output 
  #'
  #' #r2redux output
  #'
  #' #output$ratio_var (Variance of ratio)
  #' #0.08042288

  olkin_beta_ratio = function (omat,nv) {

     #aova in p158 in Olkin and Finn


f1=expression(((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 - c22 * 
    c33)) * c31)^2 /  (c22 * ((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 -
    c22 * c33)) * c31)^2 + 2 * c32 * (((c33/(c22 * c33 - c32^2)) *
    c21 + (c32/(c32^2 - c22 * c33)) * c31) * ((c32/(c32^2 - c22 *
    c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)) + c33 * ((c32/(c32^2 -
    c22 * c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)^2)     )

  c11=omat[1,1]
  c21=omat[2,1]
  c22=omat[2,2]
  c31=omat[3,1]
  c32=omat[3,2]
  c33=omat[3,3]

  av1=array(0,3)
  av1[1]=eval(D(f1,'c21'))
  av1[2]=eval(D(f1,'c31'))
  av1[3]=eval(D(f1,'c32'))

  ov=matrix(0,3,3)
  ov[1,1]=(1-omat[2,1]^2)^2/nv
  ov[2,2]=(1-omat[3,1]^2)^2/nv
  ov[3,3]=(1-omat[3,2]^2)^2/nv
  ov[2,1]=(0.5*(2*omat[3,2]-omat[2,1]*omat[3,1])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,2]^3)/nv
  ov[1,2]=ov[2,1]
  ov[3,1]=(0.5*(2*omat[3,1]-omat[2,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,1]^3)/nv
  ov[1,3]=ov[3,1]
  ov[3,2]=(0.5*(2*omat[2,1]-omat[3,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[2,1]^3)/nv
  ov[2,3]=ov[3,2]

  #variance of the difference
  aova1=t(av1)%*%ov%*%(av1)

  z=list(ratio_var=aova1)
  return(z) 

  }

