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
options(digits = 4)
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
Yt=dataf[,out[5]]
n=length(Yt)
Y=rep(NA,n)
Y[Yt=="Dead"]=1
Y[Yt=="Alive"]=0
dataf[,out[5]]=Y


# ======put together 
auc_df = NULL
for (i in 1:5){
  for (j in 1:9){
Y=dataf[,out[i]]
include=is.na(Y)==F
dat.t=dataf[include,]

X = dat.t[,pred.sets[[j]]]
Y = dat.t[,out[i]]


xp=dim(X)[2]
sum.na=function(x){sum(is.na(x))}
sna=apply(X,2,sum.na)
nmesX=names(X)


# Set missing values delta 
for(k in 1:xp){
  if(sna[k] > 0) {
    ix=as.numeric(is.na(X[,k])==F)
    X[is.na(X[,k]),k]=0
    X=data.frame(X,ix)
    nmesX=c(nmesX,paste("Imiss_",nmesX[k],sep="")) }
} 
names(X)=nmesX    

#  try this intermediate df
temp = cbind(Y, X)
acit_task <-make_sl3_Task(data = temp, outcome = "Y",covariates = colnames(temp)[!(colnames(temp) %in% "Y")],outcome_type = "binomial" )
lrn1 = Lrnr_mean$new()
lrn2 = Lrnr_glmnet$new()
lrn3 = Lrnr_randomForest$new()
sl_lrn = Lrnr_sl$new(learners = list(lrn1,lrn2,lrn3),metalearner = Lrnr_nnls$new())
acit_sl = sl_lrn$train(acit_task)
acit_sl_pred=  acit_sl$predict()
head(acit_sl_pred)
risks = acit_sl$cv_risk(loss_squared_error)
# try this =================AUC after (WORKS )
predsY = acit_sl_pred
YY = dat.t[,out[i]]
ciout=ci.cvAUC(predsY, YY)
auc = ciout$cvAUC
ci_lo = ciout$ci[1]
ci_hi = ciout$ci[2]
se = ciout$se
auc_df = rbind(auc_df,c(prednmes[j],out[i],auc,ci_lo,ci_hi,se))
  }
}
auc_df = as.tibble(auc_df)
names(auc_df) = c("PredSet","Outcome","AUC","CI_low","CI_high","SE")
#  NEXT: PLOT AUC
save(auc_df,file = "AUC.Rda")
# ===================================================
# ===============================================
# ===============end of try =============
is(iris)



# ============ sl3  works / weird output =============
acit_task <-make_sl3_Task(data = dat.t, outcome = out[4],covariates = pred.sets[[1]],outcome_type = "binomial" )
lrn1 = Lrnr_mean$new()
lrn2 = Lrnr_glmnet$new()
lrn3 = Lrnr_randomForest$new()
sl_lrn = Lrnr_sl$new(learners = list(lrn1,lrn2,lrn3),metalearner = Lrnr_nnls$new())
acit_sl = sl_lrn$train(acit_task)
acit_sl_pred=  acit_sl$predict()
head(acit_sl_pred)
risks = acit_sl$cv_risk(loss_squared_error)
# =====================





# prediction works now do AUC 
ciout = ci.cvAUC(acit_sl_pred, labels = c(1:1645))



## Previous method with super learner 
SL.library <- c("SL.glmnet", "SL.randomForest", "SL.mean","SL.rpartPrune")
#cl <- makeCluster(min(V,4))
#registerDoParallel(cl)
#n=length(Y)
#rownames(X)=paste("",1:n,sep="")

fit.test=CV.SuperLearner(Y, X, SL.library = SL.library, method = method.NNLS(), family = binomial(),V=V)
##### ROC Curve
fld=fit.test$fold
# YY / binary / labels 
YY=fit.test$Y
predsY=fit.test$SL.predict
n=length(predsY)
fold=rep(NA,n)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}
ciout=ci.cvAUC(predsY, YY, folds = fold)
txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")


#### Compare via ROC plots
preds=fit.test$SL.predict
# make prefiction object 
pred <- prediction(preds,Y)
perf1 <- performance(pred, "sens", "spec")
filn=paste(filepre,".pdf",sep="")
pdf(filn)
titl=paste("ROC with CV AUC:  ",outc," vs ",predn,sep="")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity"
     )
text(0.75,0.4,txt)
abline(0,1)
dev.off()


dat.tmp2 = data.frame(Y,X)
fit.rpart=rpart(Y ~ ., data =dat.tmp2, control = rpart.control(cp = 0.001, minsplit = 20, 
                                                               xval = 20, maxdepth = 1, minbucket = 10), method = "class")
plot(fit.rpart)

printcp(fit.rpart)
plotcp(fit.rpart)
plot(fit.rpart, uniform=TRUE, 
     main="Classification Tree for Kyphosis")





# ==========================Haven't figured out yet 

xsplt=fit.rpart$splits
varsplt = NULL
if(is.null(xsplt)) {
  xs=c(xs,NA)
  oper=c(oper,"")
}
if(is.null(xsplt)==F) {
  xname=rownames(xsplt)[1]  
  xsplt=xsplt[1,"index"]	
  x1=dat.tmp2[,xname]
  ix1=as.numeric(x1>xsplt)
  tt=table(dat.tmp2[,1],ix1)
  y1=tt[2,1]/sum(tt[,1])
  y2=tt[2,2]/sum(tt[,2])
  optmp = ifelse(y1>y2,"<",">")
  xs=c(xs,xsplt)
  oper=c(oper,optmp)
  varsplt=c(varsplt,xname)
  pp=predict(fit.rpart,type="class")
  dat.tmp2=dat.tmp2[pp==0,names(dat.tmp2)!=xname]
}
nv=length(varsplt)
if(nv = 0) {
  setsx=c(setsx,predn)
  outnames=c(outnames,outc)
  varn=c(varn,"none")
}

dd=data.frame(outnames,setsx=setsx,varnames=varn,cut=xs,oper=oper)
