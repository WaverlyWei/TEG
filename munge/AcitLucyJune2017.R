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

#inc=agrep(out[5],names(dataf),max.distance = 0.5,value=T)

### Re-do outcome

 Yt=dataf[,out[5]]
 n=length(Yt)
 Y=rep(NA,n)
 Y[Yt=="Dead"]=1
 Y[Yt=="Alive"]=0
 dataf[,out[5]]=Y

#    SL.library <- c("SL.glmnet", "SL.randomForest", "SL.mean")

SL.library <- c("SL.glmnet", "SL.randomForest", "SL.mean","SL.rpartPrune","SL.gbm")

sum.na=function(x){sum(is.na(x))}

######### 
PredNames=NULL
OutNames=NULL
AUCout=NULL
CI.out = NULL
for(i in 1:no) {
  for(j in 1:ns) {
outc=out[i]
preds=pred.sets[[j]]
predn=prednmes[j]
predictorn=pred.sets.nmes[[j]]
cat(" i = ",i," out of  ",no,":  j = ",j," out of  ",ns,"\n")
Y=dataf[,outc]
include=is.na(Y)==F
dat.t=dataf[include,]

Y=dat.t[,outc]
X=dat.t[,preds]

filepref=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")

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
    #### Fit CV.SuperLearner
      cl <- makeCluster(min(V,4))
      registerDoParallel(cl)
      n=length(Y)
      rownames(X)=paste("",1:n,sep="")
     fit.test=origami_CV.SuperLearner(Y=Y, X=X, SL.library = SL.library, method = method.NNLS(), family = binomial(),.parallel=T,V=V)
stopCluster(cl)
filen=paste("Latest",filepref,"CVAUC.pdf",sep="")
     flds=folds2foldvec(fit.test$folds)
     predsY=fit.test$SL.predict[,"SuperLearner"]
     rr=as.numeric(rownames(fit.test$SL.predict))
     YY=Y[rr]
     n=length(predsY)
     fold=rep(NA,n)
    for(k in 1:V) {
      ii=unlist(flds[k])
      fold[ii]=k
    }
    fold=fold[rr]
      ciout=ci.cvAUC(predsY, YY,folds=fold)
    txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")
    PredNames=c(PredNames,predn)
    OutNames=c(OutNames,outc)
    AUCout=c(AUCout,ciout$cvAUC)
    CI.out = rbind(CI.out,ciout$ci[1:2])
    pred <- prediction(predsY,YY)
    perf1 <- performance(pred, "sens", "spec")
    pdf(filen)
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main=paste("Results for ",outc," vs ",predn,sep=""))
text(0.75,0.4,txt,cex=0.75)
abline(0,1)
dev.off()
}
}

data.OUT=data.frame(PredNames,OutNames,AUC=AUCout,CILower=CI.out[,1],CIUpeer=CI.out[,2])

print(xtable(data.OUT,caption="AUC for Lucy/Alan Project",label="LA1",digits=4),type="latex",file="AUCOut.tex",caption.placement="top",include.rownames=T)

