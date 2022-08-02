  #' olkin_beta_inf function
  #'
  #' This function derives Information matrix for beta1 and beta2
  #' where beta1 and 2 are regression coefficients from a multiple regression model,
  #' i.e. y = x1 * beta1 + x2 * beta2 + e, where y, x1 and x2 are column-standardised   
  #' (see Olkin and Finn 1995).
  #' @references
  #' Olkin, I. and Finn, J.D. Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
  #' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
  #' @param nv Sample size
  #' @keywords information matrix in the context of correlation
  #' @export
  #' @importFrom stats D cor dnorm lm logLik pchisq qchisq qnorm
  #' @return  This function will generate information (variance-covariance) matrix of beta1 and beta2.The outputs are listed as follows.
  #' \item{info}{2x2 information (variance-covariance) matrix}
  #' \item{var1}{Variance of beta1}
  #' \item{var2}{Variance of beta2}
  #' \item{var1_2}{Variance of difference between beta1 and beta2}
  #' @examples 
  #' #To get information (variance-covariance) matrix of beta1 and beta2 where 
  #' #beta1 and 2 are regression coefficients from a multiple regression model.
  #' dat=dat1
  #' omat=cor(dat)[1:3,1:3]
  #' #omat
  #' #1.0000000 0.1958636 0.1970060
  #' #0.1958636 1.0000000 0.9981003
  #' #0.1970060 0.9981003 1.0000000
  #' 
  #' nv=length(dat$V1)
  #' output=olkin_beta_inf(omat,nv)
  #' output 
  #' 
  #' #output$info (2x2 information (variance-covariance) matrix)
  #' #0.2531406 -0.2526212
  #' #-0.2526212  0.2530269            
  #' 
  #' #output$var1 (variance of beta1)
  #' #0.2531406
  #'             
  #' #output$var2 (variance of beta2)
  #' #0.2530269
  #' 
  #' #output$var1_2 (variance of difference between beta1 and beta2)
  #' #1.01141  
  
       
  
 
   
 

  olkin_beta_inf = function (omat,nv) {

     #aova in p158 in Olkin and Finn

f=expression(((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 - c22 * 
    c33)) * c31) - ((c32/(c32^2 - c22 * c33)) * c21 + (c22/(c22 * 
    c33 - c32^2)) * c31))

f1=expression(((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 - c22 * 
    c33)) * c31))

f2=expression(((c32/(c32^2 - c22 * c33)) * c21 + (c22/(c22 * 
    c33 - c32^2)) * c31))

  c11=omat[1,1]
  c21=omat[2,1]
  c22=omat[2,2]
  c31=omat[3,1]
  c32=omat[3,2]
  c33=omat[3,3]

  av=array(0,3)
  av[1]=eval(D(f,'c21'))
  av[2]=eval(D(f,'c31'))
  av[3]=eval(D(f,'c32'))
  av1=array(0,3)
  av1[1]=eval(D(f1,'c21'))
  av1[2]=eval(D(f1,'c31'))
  av1[3]=eval(D(f1,'c32'))
  av2=array(0,3)
  av2[1]=eval(D(f2,'c21'))
  av2[2]=eval(D(f2,'c31'))
  av2[3]=eval(D(f2,'c32'))

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
  aova=t(av)%*%ov%*%(av)
  aova1=t(av1)%*%ov%*%(av1)
  aova2=t(av2)%*%ov%*%(av2)
  info=matrix(0,2,2)
  info[1,1]=aova1;info[2,2]=aova2
  info[2,1]=(aova-aova1-aova2)*-0.5;info[1,2]=info[2,1]
  cov=t(av1)%*%ov%*%(av2) #this is he same as info[2,1]
  z=list(var1_2=aova,var1=aova1,var2=aova2,info=info)
  return(z) 

  }

