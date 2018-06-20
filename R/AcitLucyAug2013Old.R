rm(list=ls())
setwd("~/Dropbox/hubbardlap/Mitch Cohen/Lucy/")
dataf=read.csv("acitLucy.csv")
#install.packages("BayesTree")
#library(BayesTree)
library(MASS)
library(SuperLearner)
library(glmnet)
library(randomForest)
library(ROCR)
library(cvAUC)

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
pred.sets=list(stand.lab,teg.ck,teg.crt,teg.ff,rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem,rotem.all)
prednmes=c("standLab","teg_ck","teg_crt","teg_ff","rotem_ex","rotem_intem","rotem_aptem","rotem_fibtem","rotem_all")
ns=length(pred.sets)
V=10
#SL.library <- c( "SL.stepAIC", "SL.glmnet", "SL.nnet", "SL.mean","SL.randomForest")
SL.library <- c( "SL.stepAIC", "SL.glmnet",  "SL.mean","SL.randomForest")
### Re-do outcomes
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

outy=outp=outauc=NULL
######### 
no=length(out)
for( i in 1:no) {
for(j in 1:ns) {
outc=out[i]
preds=pred.sets[[j]]
predn=prednmes[j]
cat(" i = ",i," out of  ",no,":  j = ",j," out of  ",ns,"\n")
### Standard measurements - mt10rbcffp24h


Y=dataf[,outc]
include=is.na(Y)==F
dat.t=dataf[include,]

Y=dat.t[,outc]
X=dat.t[,preds]
filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")

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




fit.test=CV.SuperLearner(Y,X,family=binomial(),SL.library=SL.library,V=V,verbose=F)

##### ROC Curve
fld=fit.test$fold
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
outy = c(outy,outc)
outp = c(outp,predn)
outauc = rbind(outauc,c(ciout$cvAUC,ciout$ci[1],ciout$ci[2]))
preds=fit.test$SL.predict

#### Compare via ROC plots
pred <- prediction(preds,Y)
perf1 <- performance(pred, "sens", "spec")
filn=paste(filepre,".pdf",sep="")
pdf(filn)
titl=paste("ROC with CV AUC:  ",outc," vs ",predn,sep="")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
main=titl)
text(0.75,0.4,txt)
abline(0,1)
dev.off()

}
}

#####################################################################################
### Run Everything
#####################################################################################
preds=unlist(as.vector(pred.sets))
no=length(out)
predn="Everything"
for( i in 1:no) {
outc=out[i]

cat(" i = ",i," out of  ",no,"\n")

Y=dataf[,outc]
include=is.na(Y)==F
dat.t=dataf[include,]

Y=dat.t[,outc]
X=dat.t[,preds]
filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")

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




fit.test=CV.SuperLearner(Y,X,family=binomial(),SL.library=SL.library,V=V,verbose=F)

##### ROC Curve
fld=fit.test$fold
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
outy = c(outy,outc)
outp = c(outp,"Everything")
outauc = rbind(outauc,c(ciout$cvAUC,ciout$ci[1],ciout$ci[2]))


preds=fit.test$SL.predict

#### Compare via ROC plots
pred <- prediction(preds,Y)
perf1 <- performance(pred, "sens", "spec")
filn=paste(filepre,".pdf",sep="")
pdf(filn)
titl=paste("ROC with CV AUC:  ",outc," vs ",predn,sep="")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
main=titl)
text(0.75,0.4,txt)
abline(0,1)
dev.off()
}
outAUC = data.frame(outy,outp,outauc)
names(outAUC) = c("Y","X","AUC","lower95","upperAUC")
save.image("SL.rdata")