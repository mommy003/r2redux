  #' olkin_beta_inf function
  #'
  #' This function derives Information matrix for beta1 and beta2
  #' where beta1 and 2 are regression coefficients from a multiple regression model,
  #' i.e. y = x1•beta1 + x2•beta2 + e, where y, x1 and x2 are column-standardised   
  #' (see Olkin and Finn 1995).
  #' @references
  #' Olkin, I. and J.D. Finn, Correlations redux. Psychological Bulletin, 1995. 118(1): p. 155.
  #' @param omat 3 by 3 matrix having the correlation coefficients between y, x1 and x2, i.e. omat=cor(dat) where dat is N by 3 matrix having variables in the order of cbind (y,x1,x2)
  #' @param nv sample size
  #' @keywords information matrix in the context of correlation
  #' @export
  #' @examples
  #' olkin_beta_inf (omat,nv) (see example file)
  #' dat=read.table("test_ukbb_thresholds_scaled") (see example files)
  #' omat=cor(dat)
  #' omat
  #' 1.000 0.195 0.197 0.197 0.197 0.194 0.191 0.183 0.168 0.130 0.103
  #' 0.195 1.000 0.998 0.995 0.992 0.984 0.964 0.936 0.860 0.757 0.669
  #' 0.197 0.998 1.000 0.998 0.996 0.988 0.969 0.942 0.867 0.764 0.674
  #' 0.197 0.995 0.998 1.000 0.998 0.991 0.972 0.946 0.872 0.769 0.679
  #' 0.197 0.992 0.996 0.998 1.000 0.995 0.978 0.953 0.879 0.777 0.687
  #' 0.194 0.984 0.988 0.991 0.995 1.000 0.988 0.966 0.897 0.796 0.706
  #' 0.191 0.964 0.969 0.972 0.978 0.988 1.000 0.987 0.928 0.831 0.742
  #' 0.183 0.936 0.942 0.946 0.953 0.966 0.987 1.000 0.959 0.869 0.780
  #' 0.168 0.860 0.867 0.872 0.879 0.897 0.928 0.959 1.000 0.942 0.859
  #' 0.130 0.757 0.764 0.769 0.777 0.796 0.831 0.869 0.942 1.000 0.951
  #' 0.103 0.669 0.674 0.679 0.687 0.706 0.742 0.780 0.859 0.951 1.000


  olkin_beta_inf = function (omat,nv) {

     #aøa in p158 in Olkin and Finn

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

  ø=matrix(0,3,3)
  ø[1,1]=(1-omat[2,1]^2)^2/nv
  ø[2,2]=(1-omat[3,1]^2)^2/nv
  ø[3,3]=(1-omat[3,2]^2)^2/nv
  ø[2,1]=(0.5*(2*omat[3,2]-omat[2,1]*omat[3,1])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,2]^3)/nv
  ø[1,2]=ø[2,1]
  ø[3,1]=(0.5*(2*omat[3,1]-omat[2,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,1]^3)/nv
  ø[1,3]=ø[3,1]
  ø[3,2]=(0.5*(2*omat[2,1]-omat[3,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[2,1]^3)/nv
  ø[2,3]=ø[3,2]

  #variance of the difference
  aøa=t(av)%*%ø%*%(av)
  aøa1=t(av1)%*%ø%*%(av1)
  aøa2=t(av2)%*%ø%*%(av2)
  info=matrix(0,2,2)
  info[1,1]=aøa1;info[2,2]=aøa2
  info[2,1]=(aøa-aøa1-aøa2)*-0.5;info[1,2]=info[2,1]
  cov=t(av1)%*%ø%*%(av2) #this is he same as info[2,1]
  z=list(var1_2=aøa,var1=aøa1,var2=aøa2,info=info)
  return(z) 

  }

