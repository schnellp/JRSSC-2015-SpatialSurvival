source("MCMCfx_Gauss_PCURE.NUfixed.R")
source("MCMCfx_PS_factor_PCURE.R")
source("MCMCfx_Ind_factor_PCURE.R")


library(fields)
m<-5 # each spatial row has m teeth
L<-ns<-2*m # 2xm grid of teeth, (ns locations)
nsub<-50 # number of subjects

p<-4 # number of regression parameters (beta)
iters<-10000 # MCMC sample size
burn<-2000 # MCMC burn-in period
nthin<-3 # thinning - keep 1 in nthin MCMC steps
update<-100 # update plots every 'update' steps
nsims<-20 # number of simulations to do

alpha<-.7 # PS parameter
rho<-c(1,1) # bandwidth (gauss) and spatial range (PS)
shape<-c(3) # weibull shape param
scale<-c(2) # weibull scale param
pcure <- 0.15 # percentage of teeth NOT susceptible

# set up spatial configuration (positive stable model)
s<-knots<-as.matrix(expand.grid(1:m,1:2))
type<-rep(1,ns)
beta<-c(0, 1, 0, 1)
d2<-array(0,c(ns,L,2))
d2[,,1]<-rdist(s[,1],knots[,1])^2
d2[,,2]<-rdist(s[,2],knots[,2])^2
d2[d2<0.00001]<-0
Wa<-make.W(d2[,,1]+d2[,,2])^(1/alpha)


B<-array(0,c(nsims,3,p,5))
dimnames(B)[[1]]<-paste("Dataset",1:nsims,sep="")
dimnames(B)[[2]]<-c("PS","Gauss","Ind")
dimnames(B)[[3]]<-paste("beta",1:p,sep="")
dimnames(B)[[4]]<-c("mean","sd","q05","q50","q95")

for(sim in 1:nsims){
  
  Y<-matrix(0,ns,nsub)
  X<-array(0,c(ns,nsub,p))
  for(t in 1:nsub){
    X[,t,1] <- rlnorm(1)
    X[,t,2] <- rnorm(ns)
    X[,t,3] <- rbinom(1, 1, 0.7)
    X[,t,4] <- rbinom(ns, 1, 0.7)
    A<-rPS(L,alpha)
    theta<-Wa%*%A
    SSS<-(log(theta)+X[,t,]%*%beta)/shape
    SSS<-exp(-SSS)*scale
    Y[,t]<-rweibull(ns,shape,SSS) + sample(c(0, Inf), ns, TRUE, c((1-pcure), pcure))
  }  
  C<-matrix(FALSE,ns,nsub)
  # 30% censoring rate
  cutoff <- quantile(Y, 1.00 - 0.30)
  LOW<-ifelse(Y<cutoff,Y,cutoff)
  HIGH<-ifelse(Y<cutoff,Y,Inf)

  print(paste("Beginning PS Fit", sim, "at", Sys.time()))
  source("MCMCfx_PS_factor_PCURE.R")
  fit1<-SpatSurvPS(LOW,HIGH,X,s,type,L=ns,
        common_scale=TRUE,
        iters=iters,burn=burn,nthin=nthin,update=update)
  b<-fit1$beta[burn:iters,]
  B[sim,1,,1]<-apply(b,2,mean)
  B[sim,1,,2]<-apply(b,2,sd)
  B[sim,1,,3]<-apply(b,2,quantile,0.05)
  B[sim,1,,4]<-apply(b,2,quantile,0.50)
  B[sim,1,,5]<-apply(b,2,quantile,0.95)
  
  print(paste("Beginning Gauss Fit", sim, "at", Sys.time()))
  source("MCMCfx_Gauss_PCURE.NUfixed.R")
  fit2<-SpatSurvGauss(LOW,HIGH,X,s,type,nblocks=5,
       iters=iters,burn=burn,update=update,nthin=nthin)
  b<-fit2$beta[burn:iters,]
  B[sim,2,,1]<-apply(b,2,mean)
  B[sim,2,,2]<-apply(b,2,sd)
  B[sim,2,,3]<-apply(b,2,quantile,0.05)
  B[sim,2,,4]<-apply(b,2,quantile,0.50)
  B[sim,2,,5]<-apply(b,2,quantile,0.95)
  
  print(paste("Beginning Ind Fit", sim, "at", Sys.time()))
  source("MCMCfx_Ind_factor_PCURE.R")
  fit3<-SpatSurvInd(LOW, HIGH, X, s, type,
                    iters=iters, burn=burn, update=update,
                    nthin=nthin)
  b<-fit3$beta[burn:iters,]
  B[sim,3,,1]<-apply(b,2,mean)
  B[sim,3,,2]<-apply(b,2,sd)
  B[sim,3,,3]<-apply(b,2,quantile,0.05)
  B[sim,3,,4]<-apply(b,2,quantile,0.50)
  B[sim,3,,5]<-apply(b,2,quantile,0.95)

  save(B, file=paste("results/partialResult.", sim, ".RObject", sep=""))
  
  if(sim == 1) {
    fit <- list(fit1=fit1, fit2=fit2, fit3=fit3)
    save(fit, file="results/fits.RObject")
  }
}

MAD<-matrix(0,p,3)
dimnames(MAD)[1]<-dimnames(B)[3]
dimnames(MAD)[2]<-dimnames(B)[2]

MSE<-SD<-COV<-BIAS<-MAD
for(j in 1:p){for(k in 1:3){
  MSE[j,k]  <- mean((beta[j]-B[,k,j,1])^2)
  MAD[j,k]  <- mean(abs(beta[j]-B[,k,j,4]))
  SD[j,k]   <- mean(B[,k,j,2])
  BIAS[j,k] <- mean(beta[j]-B[,k,j,1])
  COV[j,k]  <- mean((beta[j]>B[,k,j,3]) & (beta[j]<B[,k,j,5]))
}}

testSimResults.pcure <- list(MAD=MAD,
                       MSE=MSE,
                       SD=SD,
                       BIAS=BIAS,
                       COV=COV)

save(testSimResults.pcure, file="results/testSimResults.pcure.RObject")
