rm(list=ls())
dataf=read.csv("acitLucy.csv")
library(data.table)
library(SuperLearner)
library(sl3)
library(origami)
library(R6)
library(tidyverse)
library(glmnet)
library(randomForest)
### List of Outcomes
out=c("mt10rbcffp24h", "rbc_by_6h_bin", "ffp_by_6h_bin","mortalityat24h", "mortalityatdisch")	

stand.lab=c("hr0_inr","hr0_ptt")
#compared to the following groupings

#For TEG CK parameters:
teg.ck=c("sample0h_ck_sp","sample0h_ck_r","sample0h_ck_k","sample0h_ck_alpha","sample0h_ck_ma","sample0h_ck_g",
         "sample0h_ck_ly30","sample0h_ck_ly60")

#For TEG CRT parameters:
teg.crt=c("sample0h_crt_act","sample0h_crt_sp","sample0h_crt_r","sample0h_crt_k","sample0h_crt_alpha","sample0h_crt_ma","sample0h_crt_g","sample0h_crt_ly30")

#For TEG FUNCTIONAL FIBRINOGEN parameters:
teg.ff=c("sample0h_ff_flev","sample0h_ff_ma")

#For ROTEM EXTEM parameters:
rotem.ex=c("sample0h_ex_ct","sample0h_ex_cft","sample0h_ex_alpha","sample0h_ex_a10","sample0h_ex_a20",
           "sample0h_ex_mcf","sample0h_ex_ml")

#FOR ROTEM INTEM parameters:
rotem.intem=c("sample0h_in_ct","sample0h_in_cft","sample0h_in_alpha","sample0h_in_a10","sample0h_in_a20",
              "sample0h_in_mcf","sample0h_in_ml")

#For ROTEM APTEM parameters:
rotem.aptem=c("sample0h_ap_ct","sample0h_ap_cft","sample0h_ap_alpha","sample0h_ap_a10","sample0h_ap_a20",
              "sample0h_ap_mcf","sample0h_ap_ml")

#For ROTEM FIBTEM parameters:
rotem.fibtem=c("sample0h_fib_ct","sample0h_fib_alpha","sample0h_fib_a10","sample0h_fib_a20","sample0h_fib_mcf",
               "sample0h_fib_ml")
rotem.all=c(rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem)

# ==============
# Covariates
pred.sets=list(stand.lab,teg.ck,teg.crt,teg.ff,rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem,rotem.all)
prednmes=c("standLab","teg_ck","teg_crt","teg_ff","rotem_ex","rotem_intem","rotem_aptem","rotem_fibtem","rotem_all")
ns=length(pred.sets)


### Make outcomes binary (1st outcome is already binary)
Yt=dataf[,out[2]]
n=length(Yt)
Y=rep(NA,n)
Y[Yt=="Yes"]=1
Y[Yt=="No"]=0
dataf[,out[2]]=Y

Yt=dataf[,out[3]]
n=length(Yt)
Y=rep(NA,n)
Y[Yt=="Yes"]=1
Y[Yt=="No"]=0
dataf[,out[3]]=Y

Yt=dataf[,out[4]]
n=length(Yt)
Y=rep(NA,n)
Y[Yt=="Dead"]=1
Y[Yt=="Alive"]=0
dataf[,out[4]]=Y

Yt=dataf[,out[5]]
n=length(Yt)
Y=rep(NA,n)
Y[Yt=="Dead"]=1
Y[Yt=="Alive"]=0
dataf[,out[5]]=Y

# ========process missing data 
# { confusion? -> } should be at the end of ROC curves
# each round feeds one outcome 

predname=unlist(as.vector(pred.sets))
no=length(out)
predn="Everything"
for( i in 1:no) {
  outc=out[i]
 
  cat(" i = ",i," out of  ",no,"\n")
  
  Y=dataf[,outc]
  include=is.na(Y)==F
  dat.t=dataf[include,]
  
  Y=dat.t[,outc]
  X=dat.t[,predname]
  filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")
}  
  xp=dim(X)[2]
  sum.na=function(x){sum(is.na(x))}
  sna=apply(X,2,sum.na)
  nmesX=names(X)
  
  for(k in 1:xp){
    if(sna[k] > 0) {
      ix=as.numeric(is.na(X[,k])==F)
      X[is.na(X[,k]),k]=0
      X=data.frame(X,ix)
      nmesX=c(nmesX,paste("Imiss_",nmesX[k],sep="")) }
  } 
  names(X)=nmesX    
  

# fit with super learners
  acit_task <-make_sl3_Task(data = dataf,
                            outcome = outc,
                            covariates = predname[1:10] )
  # choose learners
lrn1 = Lrnr_mean$new()
lrn2 = Lrnr_glmnet$new()
lrn3 = Lrnr_randomForest$new()
sl_lrn = Lrnr_sl$new(learners = list(lrn2),metalearner = Lrnr_nnls$new())


# fit
acit_sl = sl_lrn$train(acit_task)


# === try w/o super learner 
lrn_glm = make_learner(Lrnr_glm)
lrn_fit = lrn_glm$train(acit_task)
# failed, try screening

# ===== try screening
scr_cov = Lrnr_pkg_SuperLearner_screener$new("screen.corP")
scr_fit = scr_cov$train(acit_task)




