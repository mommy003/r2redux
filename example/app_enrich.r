
  #dat=read.table("example.dat")
  #dat=read.table("example.dat2")
  #dat=read.table("PRS_UKBB_JBB_REGU_NON2")
  #dat=read.table("test.dat")
  #dat=read.table("stand_height_reg_nonreg.txt")
  #dat=read.table("PRS_UKBB_BMI_ALL_THRES_R2_scaled")
  dat=read.table("PRS_JBB_BMI_ALL_THRES_R2_scaled")
  nv=length(dat$V1)

  source("r2_diff2.r")
  source("r2_enrich2.r")
  source("r2_enrich_beta.r")

  #to get variance between 1 and 2 where 1=R2(lm(y~x2)) and 2=R2(lm(y~x1)) 
  v1=c(2)
  v2=c(1)
  res1_2=r2_diff(dat,v1,v2,nv)

  #to get variance between 1 and 2 where 1=R2(lm(y~x1+x2)) and 2=R2(lm(y~x1)) 
  v1=c(1,2)
  v2=c(1)
  res12_1=r2_diff(dat,v1,v2,nv)
  
  #to get variance between 1 and 2 where 1=R2(lm(y~x1+x2)) and 2=R2(lm(y~x2)) 
  v1=c(1,2)
  v2=c(2)
  res12_2=r2_diff(dat,v1,v2,nv)

  #enrichment analysis, 
  #e.g. R2(lm(y~x1+x2)) - R2(lm(y~x1)) vs. R2(lm(y~x1+x2)) - R2(lm(y~x2)) 
  v1=1;v2=2
  #v1=2;v2=1
  expected_ratio=0.07 #0.5
  #expected_ratio=0.93 #0.5
  enrich1_2=r2_enrich(dat,v1,v2,nv,expected_ratio)
  enrich_beta1_2=r2_enrich_beta(dat,v1,v2,nv,expected_ratio)





