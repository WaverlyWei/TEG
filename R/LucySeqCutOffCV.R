rm(list=ls())
setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/mitch cohen/lucy/")
dataf=read.csv("acitLucy.csv")
#install.packages("BayesTree")
#library(BayesTree)
library(MASS)
library(cvTools)
library(ROCR)
library(cvAUC)
library(rpart)

### List of Outcomes
out=c("mt10rbcffp24h", "rbc_by_6h_bin", "ffp_by_6h_bin",
      "mortalityat24h", "mortalityatdisch")	

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
#####################################################################################
#####  X-validate a tuning parameter alpha which is the improvement in the misclassification rate for adding an additional variable
#####
nalpha=10
alpha=seq(1.1,2,length=nalpha)
##### Start with teg.ck and outcome is rbc_by_6h_bin
nout=length(out)
nsets=length(pred.sets)-1
############################
oper=NULL
xs=NULL
setsx=NULL
varn=NULL
outnames=NULL

for(j in 1:nout) {
	for(k in 1:nsets) {
outc=out[j]
preds=pred.sets[[k]]
predn=prednmes[[k]]
Y=dataf[,outc]
X=dataf[,preds]
nx=length(preds)
for(i in 1:nx) {
	X[,i]=as.numeric(X[,i])
}
dat.tmp=data.frame(Y,X)
dat.tmp=na.omit(dat.tmp)
Y=dat.tmp[,1]
X=dat.tmp[,preds]
nx=length(preds)
n=length(Y)
folds=cvFolds(n, K = V, R = 1,type = "random")$which
for(ii in 1:V) {
  Xt = X[folds!=ii,]
  Yt = Y[folds!=ii]
  Xv = X[folds==ii,]
  Yv = Y[folds==ii]
  oper=NULL
  xs=NULL
  Xn=Xt
#  Ypred=rep(0,length(Y))
### Concatenate a rule this "X1 < 10 | X2 > ... "
dat.tmp2=data.frame(Yt,Xn)
  py=0
  Y.1=NULL
  Yp.1=NULL
  for(i in 1:nx){
    #dat.tmp2=dat.tmp[,c(1,(i+1))]
    fit.rpart=rpart(Yt~.,data=dat.tmp2,control = 
      rpart.control(cp = 0.001, minsplit = 20,xval = 20, 
       maxdepth = 1, minbucket = 10), method = "class")    
    xsplt=fit.rpart$splits
    if(is.null(xsplt)) {
      xs=c(xs,-100)
      oper=c(oper,">")
    }
    if(is.null(xsplt)==F) {
      xsplt=xsplt[1,"index"]  
      x1=dat.tmp2[,2]
      ix1=as.numeric(x1>xsplt)
      tt=table(dat.tmp2[,1],ix1)
      y1=tt[2,1]/sum(tt[,1])
      y2=tt[2,2]/sum(tt[,2])
      optmp = ifelse(y1>y2,"<",">")
      xs=c(xs,xsplt)
      oper=c(oper,optmp)
      pp=predict(fit.rpart,type="class")
      dat.tmp2=dat.tmp2[pp==0,]
### Predict on validation sample
      py=py+predict(fit.rpart,type="class",newdata=Xv)
    }
Yp.1=c(Yp.1,as.numeric(py>0))
Y.1=c(Y.1,Yv)
  }
  

  logic.mod=logreg(Yt, Xt,select=1)
  predY=c(predY,predict(logic.mod, newbin = Xv))
  outY = c(outY,Yv)
  foldout=c(foldout,rep(i,rep(length(Yv))))
}

Xn=Xt
dat.tmp2=data.frame(Yt,Xn)
i=0
stopiter = 0
ppv=NULL
while(i <nx & stopiter==0) {
i=i+1
fit.rpart=rpart(Yt~.,data=dat.tmp2,control = rpart.control(cp = 0.001, minsplit = 20, 
    xval = 20, maxdepth = 1, minbucket = 10), method = "class")
xsplt=fit.rpart$splits
if(is.null(xsplt)) {
  stopiter=1
  ppv = cbind(ppv,as.numeric(predict(fit.rpart,newdata=Xv,type="class"))-1)  
}
if(is.null(xsplt)==F) {
  xsplt=xsplt[1,"index"]	
  xname=rownames(xsplt)[1]  
  pp=predict(fit.rpart,type="class")-1
  dat.tmp2=dat.tmp2[pp==0,names(dat.tmp2)!=xname]
  ppv = cbind(ppv,predict(fit.rpart,newdata=Xv,type="class")-1)
}
}
Yobs =c(Yobs,Yv)
Ypred=c(Ypred,apply(ppv,1,max))
}
}

dd=data.frame(outnames,setsx=setsx,varnames=varn,cut=xs,oper=oper)
setwd("~/Dropbox (UC Berkeley Biostat)/hubbardlap/Mitch Cohen/Lucy/Results/CART Results/CutOffs/")

write.table(dd, file="CutOffs2.csv",sep=",",quote = TRUE,na=".",row.names = FALSE)



fit.rpart=rpart(Yt~.,data=dat.tmp2,control = rpart.control(cp = 0.001, minsplit = 20, 
                                                           xval = 20, maxdepth = 10, minbucket = 10), method = "class")
