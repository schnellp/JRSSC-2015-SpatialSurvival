load("NunnData.RData")

High[High==0]<-NA
M<-is.na(Low+High+apply(X,1:2,sum))

L<-Low
L[M]<-NA
U<-High
U[M]<-NA

#######################################################:
######       SUMMARIZE THE FINAL MODEL          #######:
#######################################################:

load("fit.Robject")

#Extremal coefficient:
num<-1:32
EC<-fit$EC.mn
EC<-cbind(NA,EC[,1:14],NA,NA,EC[,15:28],NA)
EC<-rbind(NA,EC[1:14,],NA,NA,EC[15:28,],NA)

library(fields)
image.plot(1:32,1:32,EC,
           xlab="Tooth number",ylab="Tooth number",
           cex.lab=1,cex.axis=1)


X11()

burn<-40000
iters<-50000

#Table of coefficients
SSS<-fit$beta
for(j in 1:9){
   SSS<-cbind(SSS,fit$alpha*fit$beta[,j])
}
SSS<-cbind(SSS,fit$scale)
SSS<-cbind(SSS,fit$shape[,1])
SSS<-cbind(SSS,fit$alpha*fit$shape[,1])
SSS<-cbind(SSS,fit$alpha)

SSS<-SSS[burn:iters,]
colnames(SSS)<-c(colnames(fit$beta),colnames(fit$beta),
                "scale Molar","scale Pre-Molar","scale Can","scale Inc",
                 "shapeC","shapeM","alpha")
QQQ<-apply(SSS,2,quantile,c(0.5,0.05,0.95))
print(round(QQQ,2))










#######################################################:
######           Plot one sub                   #######:
#######################################################:


tooth<-c(2:15,18:31)
covs<-c("Crown-to-root   ","Prob depth   ","Mobility   ","Miss neighbor   ")
sub<-9
sub<-3
par(mfrow=c(3,2),mar=c(4,5,4,1))
XXX<-X[,sub,5:8]
Y1<-fit$pred.y.mn[,sub]-L[,sub]
#Y2<-pred.y.mn.gauss[,sub]-L[,sub]
S5<-fit$Surv5[,sub]
S10<-fit$Surv10[,sub]

Y1<-ifelse(Y1>100,100,Y1)
#Y2<-ifelse(Y2>100,100,Y2)
M<-is.na(L[,sub])
XXX[M,]<-NA
S5[M]<-S10[M]<-Y1[M]<-Y2[M]<-NA



#Probing depth:
plot(tooth,XXX[,2],xaxt="n",main="Probing depth",
     xlab="",ylab="Probing depth (mm)")
axis(1,tooth,tooth)
for(j in 1:28){
  if(M[j]){abline(v=tooth[j],lty=2)}
  if(!M[j]){if(L[j,sub]==U[j,sub]){
   abline(v=tooth[j],lty=1)
  }}
}

#Mobility:
plot(tooth,XXX[,3],xaxt="n",main="Mobility",
     xlab="",ylab="Mobility")
axis(1,tooth,tooth)
for(j in 1:28){
  if(M[j]){abline(v=tooth[j],lty=2)}
  if(!M[j]){if(L[j,sub]==U[j,sub]){
   abline(v=tooth[j],lty=1)
  }}
}


#Mean residual life
plot(tooth,Y1,ylim=c(0,100),xaxt="n",main="MRL, PS model",
     xlab="",ylab="Mean residual life (years)")
axis(1,tooth,tooth)
for(j in 1:28){
  if(M[j]){
   abline(v=tooth[j],lty=2)
  }
  if(!M[j]){if(L[j,sub]==U[j,sub]){
   abline(v=tooth[j],lty=1)
  }}
}


#Mean residual life
# plot(tooth,Y2,ylim=c(0,100),xaxt="n",main="MRL, Gaussian model",
#      xlab="",ylab="Mean residual life (years)")
# axis(1,tooth,tooth)
# for(j in 1:28){
#   if(M[j]){
#    abline(v=tooth[j],lty=2)
#   }
#   if(!M[j]){if(L[j,sub]==U[j,sub]){
#    abline(v=tooth[j],lty=1)
#   }}
# }


#Survival probs
plot(tooth,S5,ylim=0:1,xaxt="n",main = "5-year survival probability",
     xlab="Tooth Number",ylab="Probability")
axis(1,tooth,tooth)
for(j in 1:28){
  if(M[j]){
   abline(v=tooth[j],lty=2)
  }
  if(!M[j]){if(L[j,sub]==U[j,sub]){
   abline(v=tooth[j],lty=1)
  }}
}


#Survival probs
plot(tooth,S10,ylim=0:1,xaxt="n",main="10-year survival probability",
     xlab="Tooth Number",ylab="Probability")
axis(1,tooth,tooth)
for(j in 1:28){
  if(M[j]){
   abline(v=tooth[j],lty=2)
  }
  if(!M[j]){if(L[j,sub]==U[j,sub]){
   abline(v=tooth[j],lty=1)
  }}
}


S10<-fit$Surv10
S10[is.na(L)]<-NA
apply(S10<0.9,2,sum,na.rm=T)






#######################################################:
######                     PIT                  #######:
#######################################################:
X11()

pcure <- mean(fit$pcure[burn:iters])

betaC<-QQQ[2,1:9]
betaM<-QQQ[2,10:18]
log_theta<-fit$log.theta.mn
XbC<-XbM<-L
for(t in 1:99){
  XbC[,t]<-X[,t,]%*%betaC
  XbM[,t]<-X[,t,]%*%betaM
}
scale<-QQQ[2,19:22]
shapeC<-QQQ[2,23]
shapeM<-QQQ[2,24]



par(mfrow=c(1,2))
PITC<-PITM<-L
nreps<-10

plot(NA,
     xlab="Expected quantiles",
     ylab="Observed quantiles",
     main="Conditional",
     xlim=0:1,ylim=0:1,lty=2,type="l",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25)

for(rep in 1:nreps){
 for(t in 1:99){
  set.seed(0820*t+10000*rep)
  SSS<-exp(-(log_theta[,t]+XbC[t])/shapeC)*scale[type]
  LLL<-(1-pcure)*pweibull(L[,t],shapeC,SSS)
  UUU<-(1-pcure)*pweibull(U[,t],shapeC,SSS)
  PITC[,t]<-runif(28,LLL,UUU)
 }
  
 pitc<-sort(as.vector(PITC))
 pite<-seq(0,1,length=length(pitc))
 lines(pite,pitc)#,lty=2)
}
abline(0,1,col=2,lwd=2)


par(mfrow=c(1,1))
plot(NA,
     xlab="Expected quantiles",
     ylab="Observed quantiles",
     xlim=0:1,ylim=0:1,lty=2,type="l",
     cex.lab=1.25,cex.axis=1.25,cex.main=1.25)

for(rep in 1:nreps){
 for(t in 1:99){
  set.seed(0820*t+10000*rep)
  SSS<-exp(-XbM[t]/shapeM)*scale[type]
  LLL<-(1-pcure)*pweibull(L[,t],shapeM,SSS)
  UUU<-(1-pcure)*pweibull(U[,t],shapeM,SSS)
  PITM[,t]<-runif(28,LLL,UUU)

 }
 pitm<-sort(as.vector(PITM))
 pite<-seq(0,1,length=length(pitm))
 lines(pite,pitm)#,lty=2)
}
lines(0:1,0:1,col=2,lwd=2)







