load("NunnData.RData")
source("MCMCfx_PS_factor_PCURE.R")
 
iters <- 100000 #Number of MCMC samples
burn  <- 20000 #Length of burn in
M     <- 20    #Number of latent factors

L<-Low
L[is.na(L)]<-0
U<-High
U[is.na(U)]<-1000
U[U==Inf]<-1000

nnn<-is.na(apply(X,1:2,sum)) | U==0
L[nnn]<-0
U[nnn]<-1000
for(j in 1:dim(X)[3]){
   XX<-X[,,j]
   XX[nnn]<-0
   X[,,j]<-XX
}

#Fit the model:

fit<-SpatSurvPS(L,U,X,S,type,L=M,
                iters=iters,burn=burn,update=100,
                verbose=TRUE)

save(fit, file="fit.Robject")
