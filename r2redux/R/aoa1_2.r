

  olkin1_2 = function (omat,nv) {

     #aøa in p158 in Olkin and Finn
  av=array(0,3)
  av[1]=2*omat[2,1]
  av[2]=-2*omat[3,1]
  av[3]=0

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
  #return(aøa) 

  }

