rm(list=ls())
setwd("~/Dropbox/hubbardlap/Mitch Cohen/Lucy/")
dataf=read.csv("acitLucy.csv")
library(LogicReg)
library(ROCR)
library(SuperLearner)
library(cvTools)
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
rotemnmes=c(rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes)
pred.sets=list(stand.lab,teg.ck,teg.crt,teg.ff,rotem.ex,rotem.intem,rotem.aptem,rotem.fibtem,rotem.all)
pred.sets.nmes=list(standlabnmes,tegcknmes,tegcrtnmes,tegffnmes,rotemexnmes,rotemintemnmes,rotemaptemnmes,rotemfibtemnmes,rotemnmes)

prednmes=c("standLab","teg_ck","teg_crt","teg_ff","rotem_ex","rotem_intem","rotem_aptem","rotem_fibtem","rotem_all")
ns=length(pred.sets)
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

setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Lucy/Results/Logic Results/")

sum.na=function(x){sum(is.na(x))}

V=10
logic.cv=function(Y,X,V,folds) {
  predY = NULL
  outY = NULL
  foldout=NULL
  for(i in 1:V) {
    Xt = X[folds!=i,]
    Yt = Y[folds!=i]
    Xv = X[folds==i,]
    Yv = Y[folds==i]
    logic.mod=logreg(Yt, Xt,select=1)
    predY=c(predY,predict(logic.mod, newbin = Xv))
    outY = c(outY,Yv)
    foldout=c(foldout,rep(i,rep(length(Yv))))
  }
    return(list(predY=predY,obsY=outY,foldout=foldout))}
######### 
setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Lucy/Results/Logic Results/CV/")
no=length(out)
outy=outp=outauc=NULL
for( i in 1:no) {
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
  Xna=apply(X,1,sum.na)
  X=X[Xna==0,]
  Y=Y[Xna==0]

filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")

xp=dim(X)[2]
qt=apply(na.omit(X),2,quantile,probs=seq(0.1,0.9,0.2))
newX=NULL
coln=NULL
varn=predictorn
for(k in 1:xp) {
  xn=cut(X[,k],breaks=unique(c(min(X[,k])-0.1,qt[,k],max(X[,k]+0.1))))
  inds <- model.matrix(~ factor(xn) - 1)[,-1]
  nmes=colnames(inds)
  nc=nchar(nmes)
  nmes=paste(varn[k],substr(nmes,11,nc),sep="")
  coln=c(coln,nmes)
  newX=cbind(newX,inds)
}
colnames(newX)=coln
ncol=dim(newX)[2]
n=length(Y)
folds=cvFolds(n, K = V, R = 1,type = "random")$which
folds=folds[sample(1:n,n,replace=T)]
cv.pred=logic.cv(Y,newX,V,folds)
##### ROC Curve
fld=cv.pred$foldout
YY=cv.pred$obsY
predsY=cv.pred$predY
ciout=try(ci.cvAUC(predsY, YY, folds = fld),silent=T)
if(length(ciout)>1) {
txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")
outy = c(outy,outc)
outp = c(outp,predn)
outauc = rbind(outauc,c(ciout$cvAUC,ciout$ci[1],ciout$ci[2]))

#### Compare via ROC plots
pred <- prediction(predsY,YY)
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
}
#####################################################################################
### Run Everything
#####################################################################################
predall=unlist(as.vector(pred.sets))
no=length(out)
predn="Everything"
for( i in 1:no) {
  outc=out[i]
  cat(" i = ",i," out of  ",no,"\n")
  
  
  Y=dataf[,outc]
  include=is.na(Y)==F
  dat.t=dataf[include,]
  
  Y=dat.t[,outc]
  X=dat.t[,predall]
  Xna=apply(X,1,sum.na)
  X=X[Xna==0,]
  Y=Y[Xna==0]
  
  filepre=paste("OUTis_",outc,"_PREDSETis_",predn,sep="")
  
  xp=dim(X)[2]
  qt=apply(na.omit(X),2,quantile,probs=seq(0.1,0.9,0.2))
  newX=NULL
  coln=NULL
  varn=predall
  for(k in 1:xp) {
    xn=cut(X[,k],breaks=unique(c(min(X[,k])-0.1,qt[,k],max(X[,k]+0.1))))
    inds <- model.matrix(~ factor(xn) - 1)[,-1]
    nmes=colnames(inds)
    nc=nchar(nmes)
    nmes=paste(varn[k],substr(nmes,11,nc),sep="")    
    coln=c(coln,nmes)
    newX=cbind(newX,inds)
  }
  colnames(newX)=coln
  ncol=dim(newX)[2]
  n=length(Y)
  folds=cvFolds(n, K = V, R = 1,type = "random")$which
  folds=folds[sample(1:n,n,replace=T)]
  cv.pred=logic.cv(Y,newX,V,folds)
  ##### ROC Curve
  fld=cv.pred$foldout
  YY=cv.pred$obsY
  predsY=cv.pred$predY
  ciout=try(ci.cvAUC(predsY, YY, folds = fld),silent=T)
  if(length(ciout)> 1) {
  txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")
  outy = c(outy,outc)
  outp = c(outp,predn)
  outauc = rbind(outauc,c(ciout$cvAUC,ciout$ci[1],ciout$ci[2]))
  
  #### Compare via ROC plots
  pred <- prediction(predsY,YY)
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

outAUC = data.frame(outy,outp,outauc)
names(outAUC) = c("Outcome","Predictors","AUC","lower95","upper95")
save.image("Logic.rdata")
txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")
AUC = format(round(outAUC$AUC,2),nsmall=2)
CI = paste(format(signif(outAUC$lower95,2),nsmall=2),format(signif(outAUC$upper95,2),nsmall=2),sep="-")
outtab = data.frame(outAUC[,1:2],AUC,CI)
write.table(outtab, file="LogitAucRes.csv",sep=",",quote = TRUE,na=".",row.names = FALSE)


