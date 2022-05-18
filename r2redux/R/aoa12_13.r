
  olkin12_13 = function (omat,nv) {

  #aøa in p158 in Olkin and Finn (using my own code)
  av=array(0,5)
  f=expression((c22 * ((c33/(c22 * c33 - c32^2)) * c21 + (c32/(c32^2 - 
    c22 * c33)) * c31)^2 + 2 * c32 * (((c33/(c22 * c33 - c32^2)) * 
    c21 + (c32/(c32^2 - c22 * c33)) * c31) * ((c32/(c32^2 - c22 * 
    c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)) + c33 * ((c32/(c32^2 - 
    c22 * c33)) * c21 + (c22/(c22 * c33 - c32^2)) * c31)^2) - 
    (c22 * ((c44/(c22 * c44 - c42^2)) * c21 + (c42/(c42^2 - c22 * 
        c44)) * c41)^2 + 2 * c42 * (((c44/(c22 * c44 - c42^2)) * 
        c21 + (c42/(c42^2 - c22 * c44)) * c41) * ((c42/(c42^2 - 
        c22 * c44)) * c21 + (c22/(c22 * c44 - c42^2)) * c41)) + 
        c44 * ((c42/(c42^2 - c22 * c44)) * c21 + (c22/(c22 * 
            c44 - c42^2)) * c41)^2))
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

  av[1]=eval(D(f,'c21'))
  av[2]=eval(D(f,'c31'))
  av[3]=eval(D(f,'c41'))
  av[4]=eval(D(f,'c32'))
  av[5]=eval(D(f,'c42'))


  ø=matrix(0,5,5)
  ø[1,1]=(1-omat[2,1]^2)^2/nv
  ø[2,2]=(1-omat[3,1]^2)^2/nv
  ø[3,3]=(1-omat[4,1]^2)^2/nv
  ø[4,4]=(1-omat[3,2]^2)^2/nv
  ø[5,5]=(1-omat[4,2]^2)^2/nv

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


  ø[5,1]=(0.5*(2*omat[4,1]-omat[2,1]*omat[4,2])*(1-omat[4,2]^2-omat[2,1]^2-omat[4,1]^2)+omat[4,1]^3)/nv
  ø[1,5]=ø[5,1]

  
  ø[5,2]=omat[2,1]^2+omat[4,1]^2+omat[3,2]^2+omat[4,3]^2
  ø[5,2]=ø[5,2]*0.5*omat[3,1]*omat[4,2]
  ø[5,2]=ø[5,2]+omat[2,1]*omat[4,3]+omat[4,1]*omat[3,2]
  ø[5,2]=ø[5,2]-omat[3,1]*omat[2,1]*omat[4,1]
  ø[5,2]=ø[5,2]-omat[3,1]*omat[3,2]*omat[4,2]
  ø[5,2]=ø[5,2]-omat[2,1]*omat[3,2]*omat[4,2]
  ø[5,2]=ø[5,2]-omat[4,1]*omat[4,3]*omat[4,2]
  ø[5,2]=ø[5,2]/nv 
  ø[2,5]=ø[5,2] 


  ø[5,3]=(0.5*(2*omat[2,1]-omat[4,1]*omat[4,2])*(1-omat[4,2]^2-omat[2,1]^2-omat[4,1]^2)+omat[2,1]^3)/nv
  ø[3,5]=ø[5,3]

  ø[5,4]=(0.5*(2*omat[4,3]-omat[3,2]*omat[4,2])*(1-omat[4,2]^2-omat[4,3]^2-omat[3,2]^2)+omat[4,3]^3)/nv
  ø[4,5]=ø[5,4]



  #variance of the difference
  aøa=t(av)%*%ø%*%(av)
  #return(aøa) 

  }
