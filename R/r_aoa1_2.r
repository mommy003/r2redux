

  r_olkin1_2 = function (omat,nv) {

  #for correlation, not r2
  #aova in p158 in Olkin and Finn
  av=array(0,3)
  av[1]=1  #2*omat[2,1]
  av[2]=-1  #-2*omat[3,1]
  av[3]=0

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
  #return(aova) 

  }

