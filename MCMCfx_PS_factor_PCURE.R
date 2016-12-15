############################################################
### MCMC for positive-stable model with cured proportion ###
############################################################

# Log likelihood for Weibull model
loglike <- function(LOW, HIGH, theta, Xb, shape, scale, pcure){
   SSS <- (log(theta) + Xb) / shape
   SSS <- exp(-SSS) * scale
   D   <- dweibull(LOW, shape, SSS)
   L   <- pweibull(LOW, shape, SSS)
   H   <- pweibull(HIGH, shape, SSS)
   LL  <- log(ifelse(LOW == HIGH, (1-pcure) * D, pcure + (1-pcure) * (H - L)))
   return(LL)
}

# Simulates draw from Weibull model for censored observations
samp_censored<-function(LOW,HIGH,theta,Xb,shape,scale,pcure){
   SSS   <- (log(theta) + Xb) / shape
   SSS   <- exp(-SSS) * scale
   L     <- pweibull(LOW, shape, SSS)
   H     <- pweibull(HIGH, shape, SSS)
   cured <- sample(c(0, Inf), length(L), replace=TRUE, prob=c(1-pcure, pcure))
   U     <- runif(length(L), L, H)
   Y     <- qweibull(U, shape, SSS) + cured
return(Y)}

# Constructs random effects
# Not used
make.theta<-function(FAC,logs,alpha){
    #theta is nxnF
    #s is nFxnt
    #alpha in (0,1)
    logs <- as.vector(logs)
    if(length(logs) == 1) {
     xxx <- (FAC ^ (1 / alpha)) * exp(logs)
    }
    if(length(logs) > 1) {
     xxx <- (FAC ^ (1 / alpha)) %*% exp(logs)
    }
    return(xxx)
}  

# Constructs spatial correlation matrix for Gaussian process
make.W <- function(d2, logrho){
   rho <- exp(logrho) ^ 2
   W   <- d2[, , 1] / rho[1] + d2[, , 2] / rho[2]
   W   <- exp(-0.5 * W)
   W   <-sweep(W, 1, rowSums(W), "/")   
   return(W)
}

# Alternate parameterization
make.W<-function(factor){
   W<-exp(factor)
   W<-sweep(W,1,rowSums(W),"/")   
return(W)}

# Simulates draws from positive stable distribution
rPS<-function(n,alpha,iters=1000){
   logA<-rnorm(n)
   B<-runif(n)
   curll<-dPS(logA,B,alpha)

   for(i in 1:iters){
     canlogA<-rnorm(n,logA,1/alpha)
     canll<-dPS(canlogA,B,alpha,log=T)

     R<-exp(canll-curll)
     acc<-runif(n)<R    
     logA<-ifelse(acc,canlogA,logA)
     curll<-ifelse(acc,canll,curll)
   }     
exp(logA)}

# Positive stable density
dPS<-function(logs,u,alpha,log=T){
   s<-exp(logs)
   psi<-pi*u
   a<-(sin(alpha*psi)/sin(psi))^(1/(1-alpha))
   a<-a*sin((1-alpha)*psi)/sin(alpha*psi)

   logd<-log(alpha)-log(1-alpha)+(1/(1-alpha))*log(1/s)+
         log(a)-a*(1/s^(alpha/(1-alpha)))+
         logs
   if(!log){logd<-exp(logd)}
logd}

# Inverts distance matrix
inv<-function(X,thresh=exp(-10)){
  E<-eigen(X)
  V<-E$vectors
  D<-E$values
  D<-ifelse(D<thresh,thresh,D)
  PREC<-V%*%diag(1/D)%*%t(V)
  LOGDET<-sum(-log(D))
list(PREC=PREC,LOGDET=LOGDET)}

# Y = [number sites] x [nsubs] (Follow-up time)
# C = [number sites] x [nsubs] (TRUE if censored, FALSE ow)
# X = [number sites] x [nsubs] x p
# s = [n sites] x 2
# knots = [n knots] x 2  

SpatSurvPS<-function(LOW,HIGH,X,s,type,
    spatial=TRUE,init.range=2,common_scale=F,
    common_shape=T,common_beta=T,L=5,nthin=5,
    iters=5000,burn=1000,update=10, verbose=FALSE){

    library(fields)
    library(emulator)

    #BOOKKEEPING
    ns<-nrow(s)            # number of teeth per subject
    nsub<-nsubs<-ncol(LOW) # number of subjects
    p<-dim(X)[3]           # number of covariates
    ntype<-max(type)       # number of tooth types

    d<-rdist(s,s)          # distances between teeth
    
    #INITIAL VALUES:
    beta<-rep(0,p)             # regression parameters
    scale<-rep(max(LOW),ntype) # Weibull scale
    shape<-rep(1,ntype)        # Weibull shape
    alpha<-0.3                 # positive stable parameter
    pcure<-0.5                 # cured proportion
    logrange<-rep(log(init.range),2) # spatial range for Gaussian random effects

    logA<-matrix(5,L,nsubs)
    B<-matrix(0.5,L,nsubs)
    factor<-matrix(0,ns,L)
    tauf<-1
    E<-inv(exp(-d/1))


    W<-make.W(factor)
    Wa<-W^(1/alpha)
    theta<-Wa%*%exp(logA)
    XB<-matrix(0,ns,nsubs)
    for(j in 1:nsubs){XB[,j]<-X[,j,]%*%beta}
    SCALE<-scale[type]
    SHAPE<-shape[type]

    curll<-matrix(0,ns,nsubs)
    for(t in 1:nsubs){
      curll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,SCALE,pcure)
    }
      

    keep.alpha<-rep(0,iters)
    keep.range<-matrix(0,iters,2)
    keep.shape<-keep.scale<-matrix(0,iters,ntype)
    keep.pcure<-rep(0,iters)
    keep.beta<-matrix(0,iters,p)
    colnames(keep.beta)<-dimnames(X)[[3]]

    attb<-accb<-MHb<-c(0.5,.05,.01,.5,0.05,0.01,0.005,1)
    attA<-accA<-MHA<-c(5,2,1,1,.5,.5,.5)
    attB<-accB<-MHB<-.5
    att<-acc<-MH<-rep(0.05,p+5)
    cuts<-c(1,5,10,12,15,20)

    Y1<-YS<-Y2<-Y3<-Y4<-invCPO<-S5<-S10<-S20<-S5C<-S10C<-S20C<-matrix(0,ns,nsubs)
    EC1<-EC2<-matrix(0,ns,ns)
    dev<-rep(0,iters)

    for(i in 1:iters){
     if (verbose) {
       print(paste(i, "/", iters))
     }

     for(ttt in 1:nthin){

      ##########################################################
      ##############      Random effects A and B    ############
      ##########################################################
      level<-oldA<-logA
      if(spatial){for(l in 1:L){for(t in 1:nsubs){
         MH1<-MH2<-MHA[1]
         level[l,t]<-1
         for(j in 1:length(cuts)){if(oldA[l,t]>cuts[j]){
            MH1<-MHA[j+1];level[l,t]<-j+1
         }}
         canlogA<-rnorm(1,logA[l,t],MH1)
         for(j in 1:length(cuts)){if(canlogA>cuts[j]){
            MH2<-MHA[j+1]
         }}

         diff<-exp(canlogA)-exp(logA[l,t])
         cantheta<-theta[,t]+Wa[,l]*diff
         canll<-loglike(LOW[,t],HIGH[,t],cantheta,XB[,t],SHAPE,SCALE,pcure)

         R<-sum(canll-curll[,t])+
            dnorm(logA[l,t],canlogA,MH2,log=T)-
            dnorm(canlogA,logA[l,t],MH1,log=T)+
            dPS(canlogA,B[l,t],alpha,log=T)-
            dPS(logA[l,t],B[l,t],alpha,log=T)
         if(!is.na(exp(R))){if(runif(1)<exp(R)){
            logA[l,t]<-canlogA
            curll[,t]<-canll
            theta[,t]<-cantheta
         }}
      }}}
      for(j in 1:length(MHA)){
         accA[j]<-accA[j]+sum(oldA[level==j]!=logA[level==j])
         attA[j]<-attA[j]+sum(level==j)
      }


      oldB<-B
      if(spatial){for(t in 1:nsubs){for(l in 1:L){             
         attB<-attB+1
         canB<-rnorm(1,B[l,t],MHB)
         if(canB>0.001 & canB<0.999){
          R<-dPS(logA[l,t],canB,alpha,log=T)-
             dPS(logA[l,t],B[l,t],alpha,log=T)
          if(!is.na(exp(R))){if(runif(1)<exp(R)){
             accB<-accB+1
             B[l,t]<-canB
          }}
         }
      }}}


      ##########################################################
      ##############              alpha             ############
      ##########################################################

      if(spatial & i>500){
       att[p+4]<-att[p+4]+1
       eta<-log(alpha/(1-alpha))
       caneta<-rnorm(1,eta,MH[p+4])
       canalpha<-exp(caneta)/(1+exp(caneta))
       canWa<-W^(1/canalpha)
       cantheta<-canWa%*%exp(logA)
       canll<-curll
       for(t in 1:nsubs){
        canll[,t]<-loglike(LOW[,t],HIGH[,t],cantheta[,t],XB[,t],SHAPE,SCALE,pcure)
       }
       R<-sum(dPS(logA,B,canalpha))-
          sum(dPS(logA,B,alpha))+
          sum(canll-curll)+
          dnorm(caneta,log=T)-dnorm(eta,log=T)
       if(!is.na(exp(R))){if(runif(1)<exp(R)){
          alpha<-canalpha;curll<-canll;
          Wa<-canWa;theta<-cantheta
          acc[p+4]<-acc[p+4]+1
       }}           
      }


      ##########################################################
      ##############        latent factors          ############
      ##########################################################
      for(l in 1:L){
        QQQ<- -0.5*tauf*quad.form(E$PREC,factor[,l])
        for(j in 1:ns){
         att[p+3]<-att[p+3]+1
         canfactor<-factor
         canfactor[j,l]<-rnorm(1,factor[j,l],MH[p+3])
         canW<-exp(canfactor[j,])/sum(exp(canfactor[j,]))
         canWa<-canW^(1/alpha)
         cantheta<-canWa%*%exp(logA)
         canll<-loglike(LOW[j,],HIGH[j,],cantheta,XB[j,],SHAPE[j],SCALE[j],pcure)
         canQQQ<- -0.5*tauf*quad.form(E$PREC,canfactor[,l])
         R<-canQQQ-QQQ+
            sum(canll-curll[j,])
         if(!is.na(exp(R))){if(runif(1)<exp(R)){
           factor[j,l]<-canfactor[j,l];curll[j,]<-canll;
           W[j,]<-canW;Wa[j,]<-canWa;
           QQQ<-canQQQ;theta[j,]<-cantheta
           acc[p+3]<-acc[p+3]+1
         }}            
        }
       }

       SSS<-0
       for(l in 1:L){
          SSS<-SSS+quad.form(E$PREC,factor[,l])
       }
       tauf<-rgamma(1,ns*L/2+.1,SSS/2+.1)
      
      ##########################################################
      ##############              pcure             ############
      ##########################################################
      att[p+5] <- att[p+5]+1
      canpcure <- runif(1)
      canll <- curll
      for(t in 1:nsubs) {
        canll[,t] <- loglike(LOW[,t], HIGH[,t], theta[,t], XB[,t], SHAPE, SCALE, canpcure)
      }
      R <- sum(canll - curll)
      if(!is.na(exp(R))){if(runif(1)<exp(R)){
        pcure<-canpcure;curll<-canll;
        acc[p+5]<-acc[p+5]+1
      }}

      ##########################################################
      ##############              shape             ############
      ##########################################################

      for(j in 1:ntype){if(j==1 | !common_shape){
        att[1]<-att[1]+1  #mh attempts
        canshape<-shape   #Candidate shape value vector
        canshape[j]<-exp(rnorm(1,log(shape[j]),MH[1]))  #MH: metropolis-hastings
        if(common_shape){canshape[1:ntype]<-canshape[1]}  #ntype: number of tooth types
        canSHAPE<-canshape[type]
        canll<-curll      #Set candidate log likelihood to current log likelihood
        for(t in 1:nsubs){  #nsubs: number of subjects
          canll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],canSHAPE,SCALE,pcure)
        }
        R<-sum(canll-curll)+
           dnorm(log(canshape[j]),0,10,log=T)-
           dnorm(log(shape[j]),0,10,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
           shape<-canshape;SHAPE<-canSHAPE;curll<-canll;
           acc[1]<-acc[1]+1 #mh acceptances
        }}           
      }}



      ##########################################################
      ##############              scale             ############
      ##########################################################

      for(j in 1:ntype){if(j==1 | !common_scale){
        att[2]<-att[2]+1
        canscale<-scale
        canscale[j]<-exp(rnorm(1,log(scale[j]),MH[2]))
        if(common_scale){canscale[1:ntype]<-canscale[1]}
        canSCALE<-canscale[type]
        canll<-curll
        for(t in 1:nsubs){
          canll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,canSCALE,pcure)
        }
        R<-sum(canll-curll)+
           dnorm(log(canscale[j]),0,10,log=T)-
           dnorm(log(scale[j]),0,10,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
           scale<-canscale;SCALE<-canSCALE;curll<-canll;
           acc[2]<-acc[2]+1
        }}           
      }}


      ##########################################################
      ##############        regression coefs        ############
      ##########################################################

      for(l in 1:p){
        att[l+2]<-att[l+2]+1
        canbeta<-beta
        canbeta[l]<-rnorm(1,beta[l],MH[l+2])
        canXB<-XB
        for(j in 1:nsubs){canXB[,j]<-X[,j,]%*%canbeta}
        canll<-curll
        for(t in 1:nsubs){
          canll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],canXB[,t],SHAPE,SCALE,pcure)
        }
        R<-sum(canll-curll)+
           dnorm(canbeta[l],0,10,log=T)-
           dnorm(beta[l],0,10,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
           beta<-canbeta;XB<-canXB;curll<-canll;
           acc[l+2]<-acc[l+2]+1
        }}           
      }


     }


      ##########################################################
      ##############      TUNE THE MH STEPS         ############
      ##########################################################

      for(j in 1:length(MHA)){if(i<burn/2 & attA[j]>100){
        if(accA[j]/attA[j]<0.3){MHA[j]<-MHA[j]*0.9}
        if(accA[j]/attA[j]>0.6){MHA[j]<-MHA[j]*1.1}
        accA[j]<-attA[j]<-0
      }}
      if(i<burn/2 & attB>100){
        if(accB/attB<0.3){MHB<-MHB*0.9}
        if(accB/attB>0.6){MHB<-MHB*1.1}
         accB<-attB<-0
      }
      for(j in 1:length(MH)){if(i<burn/2 & att[j]>25){
        if(acc[j]/att[j]<0.3){MH[j]<-MH[j]*0.8}
        if(acc[j]/att[j]>0.6){MH[j]<-MH[j]*1.2}
        acc[j]<-att[j]<-0
      }}



      #KEEP TRACK OF STUFF:
      keep.range[i,]<-exp(logrange)
      keep.alpha[i]<-alpha
      keep.shape[i,]<-shape
      keep.scale[i,]<-scale
      keep.beta[i,]<-beta
      keep.pcure[i]<-pcure
      dev[i]<--2*sum(curll)

      if(i>burn){
        nnn<-iters-burn
        for(t in 1:nsubs){
          Ynew<-samp_censored(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,SCALE,pcure)
          YC <- ifelse(Ynew == Inf, 0, Ynew)
          Y1[,t]<-Y1[,t]+YC
          YS[,t]<-YS[,t]+ifelse(Ynew == Inf, 0, 1)
          Y2[,t]<-Y2[,t]+YC*YC
          S5C[,t] <-S5C[,t]+(YC>LOW[,t]+5)
          S10C[,t] <-S10C[,t]+(YC>LOW[,t]+10)
          S20C[,t] <-S20C[,t]+(YC>LOW[,t]+20)
          S5[,t] <-S5[,t] +(Ynew>LOW[,t]+5)/nnn
          S10[,t]<-S10[,t]+(Ynew>LOW[,t]+10)/nnn
          S20[,t]<-S20[,t]+(Ynew>LOW[,t]+20)/nnn
        }
        Y3<-Y3+log(theta)/nnn
        Y4<-Y4+log(theta)*log(theta)/nnn
        invCPO<-invCPO+exp(-curll)/(iters-burn)

        EC<-matrix(0,ns,ns)
        for(j1 in 1:ns){for(j2 in 1:ns){
           EC[j1,j2]<-sum((Wa[j1,]+Wa[j2,])^alpha)
        }}
        EC1<-EC1+EC/(iters-burn)
        EC2<-EC2+EC*EC/(iters-burn)
      }

      #DISPLAY CURRENT VALUE:
      if(i%%update==0 && FALSE){
         par(mfrow=c(3,3))
         plot(keep.alpha[1:i],type="l",main="alpha")
         image.plot(W,main="W",xlab="s",ylab="l")
         image.plot(factor,main="factor",xlab="s",ylab="l")
         plot(keep.scale[1:i,1],type="l",main="scale[1]")
         plot(keep.shape[1:i,1],type="l",main="shape[1]")
         plot(keep.beta[1:i,p],type="l",main="beta[p,1]")
         for(t in 1:3){if(t<=nsubs){
            plot(s[,1],LOW[,t],pch=ifelse(LOW[,t]==HIGH[,t],1,2),col=s[,2],ylim=c(0,200))
            SSS<-(log(theta[,t])+XB[,t])/SHAPE
            SSS<-exp(-SSS)*SCALE
            lines(s[,1],SSS*(log(2))^(1/SHAPE),col=s[,2])
         }}
      }


    }


    beta<-colMeans(keep.beta[burn:iters,])
    if(ntype==1){
     shape<-mean(keep.shape[burn:iters,])[type]
     scale<-mean(keep.scale[burn:iters,])[type]
    }
    if(ntype>1){
     shape<-colMeans(keep.shape[burn:iters,])[type]
     scale<-colMeans(keep.scale[burn:iters,])[type]
    }
    for(j in 1:nsubs){XB[,j]<-X[,j,]%*%beta}
    theta<-exp(Y3)
    for(t in 1:nsubs){
      curll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,SCALE,pcure)
    }
    dhat<--2*sum(curll)
    dbar<-mean(dev[burn:iters])
    pD<-dbar-dhat
    DIC<-dbar+pD

    CPO<-1/invCPO
    LPML<-sum(log(CPO))

list(pred.y.mn=Y1/YS,
     pred.y.var=(Y2/YS)-(Y1/YS)^2,
     pred.y.susc=YS / (iters - burn),
     log.theta.mn=Y3,
     log.theta.var=Y4-Y3^2,
     Surv5 =S5,
     Surv10=S10,
     Surv20=S20,
     Surv5C=S5C/YS,
     Surv10C=S10C/YS,
     Surv20C=S20C/YS,
     EC.mn=EC1,
     EC.var=EC2-EC1^2,
     range=keep.range,
     alpha=keep.alpha,
     shape=keep.shape,
     scale=keep.scale,
     pcure=keep.pcure,
     beta=keep.beta,
     dev=dev,DIC=DIC,pD=pD,dbar=dbar,
     CPO=CPO,LPML=LPML)}

