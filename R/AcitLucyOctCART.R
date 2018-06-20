rm(list=ls())
setwd("~/Dropbox/hubbardlap/Mitch Cohen/Lucy/Attachments/")
dataf=read.csv("acitLucy.csv")
#install.packages("BayesTree")
#library(BayesTree)
library(SuperLearner)
library(ROCR)
library(cvAUC)
library(rpart)
library(rattle)
library(RColorBrewer)

### List of Outcomes
out=c("mt10rbcffp24h", "rbc_by_6h_bin", "ffp_by_6h_bin","mortalityat24h", "mortalityatdisch")	

stand.lab=c("hr0_inr","hr0_ptt")
standlabnmes=gsub("^.*?hr0_","",stand.lab)
#compared to the following groupings

#For TEG CK parameters:
#### Removed sample0h_ck_ly60 because almost no observations
teg.ck=c("sample0h_ck_sp","sample0h_ck_r","sample0h_ck_k","sample0h_ck_alpha","sample0h_ck_ma","sample0h_ck_g",
"sample0h_ck_ly30")
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
rotemnmes=c(paste("ex",rotemexnmes,sep="."),paste("intem",rotemintemnmes,sep="."),
paste("aptem",rotemaptemnmes,sep="."),paste("fib",rotemfibtemnmes,sep="."))
pred.sets=list(stand.lab,teg.ck,teg.crt,teg.ff,rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem,rotem.all)
pred.sets.nmes=list(standlabnmes,tegcknmes,tegcrtnmes,tegffnmes,rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes,rotemnmes)
prednmes=c("standLab","teg_ck","teg_crt","teg_ff","rotem_ex","rotem_intem","rotem_aptem","rotem_fibtem","rotem_all")

ns=length(pred.sets)
V=10
SL.library <- c( "SL.rpart")
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


######### Sets to Look at, are j = 1 and j=5
######### Outcome is i=2
no=length(out)

rpart.fun=function (Y, X, family,  cp = 0.01, minsplit = 30, 
    xval = 10, maxdepth = 30, minbucket = round(minsplit/3), 
    ...) 
{
    if (family == "gaussian") {
        fit.rpart <- rpart(Y ~ ., data = data.frame(Y, X), control = rpart.control(cp = cp, 
            minsplit = minsplit, xval = xval, maxdepth = maxdepth, 
            minbucket = minbucket), method = "anova")
    }
    if (family == "binomial") {
        fit.rpart <- rpart(Y ~ ., data = data.frame(Y, X), control = rpart.control(cp = cp, 
            minsplit = minsplit, xval = xval, maxdepth = maxdepth, 
            minbucket = minbucket), method = "class")
    }
       return(fit.rpart)
}

for(i in c(1:3,5)) {
#for(j in c(1,5,9)) {
for(j in 9) {
outc=out[i]
preds=pred.sets[[j]]
coln=pred.sets.nmes[[j]]
predn=prednmes[[j]]

cat(" i = ",i," out of  ",no,":  j = ",j," out of  ",ns,"\n")
### Standard measurements - mt10rbcffp24h


Y=dataf[,outc]
include=is.na(Y)==F
dat.t=dataf[include,]

Y=dat.t[,outc]
X=dat.t[,preds]
names(X)=coln
filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")

xp=dim(X)[2]
sum.na=function(x){sum(is.na(x))}
sna1=apply(X,1,sum.na)
include = sna1 < xp
X=X[include,]
Y=Y[include]
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



titl=paste("ROC with CV AUC:  ",outc," vs ",predn,sep="")
fit.test=CV.SuperLearner(Y,X,family=binomial(),SL.library=SL.library,V=V,verbose=F)
tree.fit=rpart.fun(Y,X,"binomial")
pdf(paste(filepre,"CART.pdf",sep=""))
fancyRpartPlot(tree.fit,main=titl)
dev.off()
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

preds2=fit.test$SL.predict

#### Compare via ROC plots
pred <- prediction(preds2,Y)
perf1 <- performance(pred, "sens", "spec")
filn=paste(filepre,"ROC.pdf",sep="")
pdf(filn)
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




fit.test=CV.SuperLearner(Y,X,family=binomial(),SL.library=SL.library,V=V,verbose=T)

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
#####################################################################################
### Add Standard to the various sets and compare
######### First run data generating steps
#####################################################################################



