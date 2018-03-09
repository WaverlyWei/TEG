rm(list=ls())
library(MASS)
library(SuperLearner)
library(origami)
library(glmnet)
library(randomForest)
library(ROCR)
library(cvAUC)
library(doParallel)
library(readstata13)
library(tmle)
library(xtable)
library(sl3)
library(data.table)
library(R6)
library(tidyverse)

set.seed(12313197)
#setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Trauma and Coagulation (White Space Conflict)/ACIT/Data")
dataf=read.dta13("ACIT_6Jan17_deidentified1.dta")
plt0to6=as.numeric(dataf$icu_0to6h_plt_units>0)
bloodUnit0to6=as.numeric(dataf$icu_0to6h_blood_units>0)
ffp0to6=as.numeric(dataf$icu_0to6h_ffp_units>0)
dataf=data.frame(dataf,plt0to6,bloodUnit0to6,ffp0to6)
rm(plt0to6,bloodUnit0to6,ffp0to6)

V=10

### List of Outcomes
out=c("bloodUnit0to6","ffp0to6","plt0to6","mortalityat24h", "mortalityatdisch")	
no=length(out)

stand.lab=c("hr0_inr","hr0_ptt")
standlabnmes=gsub("^.*?hr0_","",stand.lab)
#compared to the following groupings

#For TEG CK parameters:
#### Removed sample0h_ck_ly60 because almost no observations
teg.ck=c("sample0h_ck_sp","sample0h_ck_r","sample0h_ck_k","sample0h_ck_alpha","sample0h_ck_ma","sample0h_ck_g","sample0h_ck_ly30")
tegcknmes=gsub("^.*?sample0h_ck_","",teg.ck)

#For TEG CRT parameters:
teg.crt=c("sample0h_crt_act","sample0h_crt_sp","sample0h_crt_r","sample0h_crt_k","sample0h_crt_alpha","sample0h_crt_ma","sample0h_crt_g","sample0h_crt_ly30")
tegcrtnmes=gsub("^.*?sample0h_crt_","",teg.crt)

#For TEG FUNCTIONAL FIBRINOGEN parameters:
teg.ff=c("sample0h_ff_flev","sample0h_ff_ma")
tegffnmes=gsub("^.*?sample0h_ff_","",teg.ff)

#For ROTEM EXTEM parameters:
rotem.ex=c("sample0h_ex_ct","sample0h_ex_cft","sample0h_ex_alpha","sample0h_ex_a10","sample0h_ex_a20",
           "sample0h_ex_mcf","sample0h_ex_ml")
rotemexnmes=gsub("^.*?sample0h_ex_","",rotem.ex)

#FOR ROTEM INTEM parameters:
rotem.intem=c("sample0h_in_ct","sample0h_in_cft","sample0h_in_alpha","sample0h_in_a10","sample0h_in_a20","sample0h_in_mcf","sample0h_in_ml")
rotemintemnmes=gsub("^.*?sample0h_in_","",rotem.intem)

#For ROTEM APTEM parameters:
rotem.aptem=c("sample0h_ap_ct","sample0h_ap_cft","sample0h_ap_alpha","sample0h_ap_a10","sample0h_ap_a20",
              "sample0h_ap_mcf","sample0h_ap_ml")
rotemaptemnmes=gsub("^.*?sample0h_ap_","",rotem.aptem)

#For ROTEM FIBTEM parameters:
rotem.fibtem=c("sample0h_fib_ct","sample0h_fib_alpha","sample0h_fib_a10","sample0h_fib_a20","sample0h_fib_mcf",
               "sample0h_fib_ml")
rotemfibtemnmes=gsub("^.*?sample0h_fib_","",rotem.fibtem)

rotem.all=c(rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem)
rotemnmes=c(rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes)
pred.sets=list(stand.lab,teg.ck,teg.crt,teg.ff,rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem,rotem.all)
pred.sets.nmes=list(standlabnmes,tegcknmes,tegcrtnmes,tegffnmes,rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes,rotemnmes)

prednmes=c("standLab","teg_ck","teg_crt","teg_ff","rotem_ex","rotem_intem","rotem_aptem","rotem_fibtem","rotem_all")
ns=length(pred.sets)

#change outcome to binary 
Yt=dataf[,out[1]]
n=length(Yt)
Y=rep(NA,n)
Y[Yt=="Dead"]=1
Y[Yt=="Alive"]=0
dataf[,out[1]]=Y


Y=dataf[,out[4]]
include=is.na(Y)==F
dat.t=dataf[include,]

X = dat.t[,pred.sets[[1]]]
Y = dat.t[,out[4]]


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


#this part works 
acit_task <-make_sl3_Task(data = dat.t, outcome = out[4],covariates = pred.sets[[1]] )
lrn1 = Lrnr_mean$new()
lrn2 = Lrnr_glm$new()
lrn3 = Lrnr_randomForest$new()
sl_lrn = Lrnr_sl$new(learners = list(lrn1,lrn3),metalearner = Lrnr_nnls$new())
acit_sl = sl_lrn$train(acit_task)
acit_sl_pred=  acit_sl$predict()
head(acit_sl_pred)
risks = acit_sl$cv_risk(loss_squared_error)

# prediction works now do AUC 
ciout = ci.cvAUC(acit_sl_pred, labels = c(1:1645))


## Previous method with super learner 
SL.library <- c("SL.glmnet", "SL.randomForest", "SL.mean","SL.rpartPrune")
cl <- makeCluster(min(V,4))
registerDoParallel(cl)
n=length(Y)
rownames(X)=paste("",1:n,sep="")
fit.test=CV.SuperLearner(Y, X, SL.library = SL.library, method = method.NNLS(), family = binomial(),V=V)
rr=as.numeric(rownames(fit.test$SL.predict))
YY=Y[rr]
predsY=fit.test$SL.predict
ciout = ci.cvAUC(predsY, labels = c(1:1671))
