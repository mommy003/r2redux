
  olkin12_3 = function (omat,nv) {

  #aøa in p158 in Olkin and Finn, but this is for var(r2(0,1,2) - r2(0,3,4))

  f=expression(c22*( ( (c33/(c22*c33-c32^2))*c21 + (c32/(c32^2-c22*c33))*c31 )^2 +  2*c32*(( (c33/(c22*c33-c32^2))*c21 + (c32/(c32^2-c22*c33))*c31 ) *  ( (c32/(c32^2-c22*c33))*c21 + (c22/(c22*c33-c32^2))*c31 )) +  c33*( (c32/(c32^2-c22*c33))*c21 + (c22/(c22*c33-c32^2))*c31 )^2) - (c41^2))

  c11=omat[1,1]
  c21=omat[2,1]
  c22=omat[2,2]
  c31=omat[3,1]
  c32=omat[3,2]
  c33=omat[3,3]
  c41=omat[4,1]
  c42=omat[4,2]
  c43=omat[4,3]
  c44=omat[4,4]

  av=array(0,4)
  av[1]=eval(D(f,'c21'))
  av[2]=eval(D(f,'c31'))
  av[3]=eval(D(f,'c41'))
  av[4]=eval(D(f,'c32'))

  ø=matrix(0,4,4)
  ø[1,1]=(1-omat[2,1]^2)^2/nv
  ø[2,2]=(1-omat[3,1]^2)^2/nv
  ø[3,3]=(1-omat[4,1]^2)^2/nv
  ø[4,4]=(1-omat[3,2]^2)^2/nv

  ø[2,1]=(0.5*(2*omat[3,2]-omat[2,1]*omat[3,1])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,2]^3)/nv
  ø[1,2]=ø[2,1]

  ø[3,1]=(0.5*(2*omat[4,2]-omat[2,1]*omat[4,1])*(1-omat[4,2]^2-omat[2,1]^2-omat[4,1]^2)+omat[4,2]^3)/nv
  ø[1,3]=ø[3,1]
  ø[3,2]=(0.5*(2*omat[4,3]-omat[3,1]*omat[4,1])*(1-omat[4,3]^2-omat[3,1]^2-omat[4,1]^2)+omat[4,3]^3)/nv
  ø[2,3]=ø[3,2]

  ø[4,1]=(0.5*(2*omat[3,1]-omat[2,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[3,1]^3)/nv
  ø[1,4]=ø[4,1]
  ø[4,2]=(0.5*(2*omat[2,1]-omat[3,1]*omat[3,2])*(1-omat[3,2]^2-omat[2,1]^2-omat[3,1]^2)+omat[2,1]^3)/nv
  ø[2,4]=ø[4,2]

  ø[4,3]=omat[2,1]^2+omat[3,1]^2+omat[4,2]^2+omat[4,3]^2
  ø[4,3]=ø[4,3]*0.5*omat[4,1]*omat[3,2]
  ø[4,3]=ø[4,3]+omat[2,1]*omat[4,3]+omat[3,1]*omat[4,2]
  ø[4,3]=ø[4,3]-omat[4,1]*omat[2,1]*omat[3,1]
  ø[4,3]=ø[4,3]-omat[4,1]*omat[4,2]*omat[4,3]
  ø[4,3]=ø[4,3]-omat[2,1]*omat[4,2]*omat[3,2]
  ø[4,3]=ø[4,3]-omat[3,1]*omat[4,3]*omat[3,2]
  ø[4,3]=ø[4,3]/nv 
  ø[3,4]=ø[4,3] 

  #variance of the difference
  aøa=t(av)%*%ø%*%(av)
  #return(aøa) 

  }

