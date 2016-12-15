
loglike<-function(LOW,HIGH,Xb,shape,scale,pcure){
  SSS<-(Xb)/shape
  SSS<-exp(-SSS)*scale
  D<-dweibull(LOW,shape,SSS)
  L<-pweibull(LOW,shape,SSS)
  H<-pweibull(HIGH,shape,SSS)
  LL<-log(ifelse(LOW==HIGH, (1-pcure)*D, pcure+(1-pcure)*(H-L)))
return(LL)}

samp_censored<-function(LOW,HIGH,Xb,shape,scale,pcure){
   SSS<-(Xb)/shape
   SSS<-exp(-SSS)*scale
   L<-pweibull(LOW,shape,SSS)
   H<-pweibull(HIGH,shape,SSS)
   cured<-sample(c(0,Inf), length(L), replace=TRUE, prob=c(1-pcure, pcure))
   U<-runif(length(L),L,H)
   Y<-qweibull(U,shape,SSS)+cured
return(Y)} 

make.W<-function(d2,logrho){
   rho<-exp(logrho)^2
   W<-d2[,,1]/rho[1]+d2[,,2]/rho[2]
   W<-exp(-0.5*W)
   W<-sweep(W,1,rowSums(W),"/")   
return(W)}


make.W<-function(factor){
   W<-exp(factor)
   W<-sweep(W,1,rowSums(W),"/")   
return(W)}

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

SpatSurvInd<-function(LOW,HIGH,X,s,type,
    spatial=TRUE,init.range=2,common_scale=F,
    common_shape=T,common_beta=T,L=5,nthin=5,
    iters=5000,burn=1000,update=10){

    library(fields)
    library(emulator)

    #BOOKKEEPING
    ns<-nrow(s)
    nsub<-nsubs<-ncol(LOW)
    p<-dim(X)[3]
    ntype<-max(type)

    d<-rdist(s,s)
    
    #INITIAL VALUES:
    beta<-rep(0,p)
    scale<-rep(max(LOW),ntype)
    shape<-rep(1,ntype)
    logrange<-rep(log(init.range),2)
    pcure<-0.5
    
    XB<-matrix(0,ns,nsubs)
    for(j in 1:nsubs){XB[,j]<-X[,j,]%*%beta}
    SCALE<-scale[type]
    SHAPE<-shape[type]

    curll<-matrix(0,ns,nsubs)
    for(t in 1:nsubs){
      curll[,t]<-loglike(LOW[,t],HIGH[,t],XB[,t],SHAPE,SCALE,pcure)
    }
      

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

    Y1<-YS<-Y2<-invCPO<-matrix(0,ns,nsubs)
    dev<-rep(0,iters)

    for(i in 1:iters){

     for(ttt in 1:nthin){
      
      
      ##########################################################
      ##############              pcure             ############
      ##########################################################
      att[p+5] <- att[p+5]+1
      canpcure <- runif(1)
      canll <- curll
      for(t in 1:nsubs) {
        canll[,t] <- loglike(LOW[,t], HIGH[,t], XB[,t], SHAPE, SCALE, canpcure)
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
          canll[,t]<-loglike(LOW[,t],HIGH[,t],XB[,t],canSHAPE,SCALE,pcure)
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
          canll[,t]<-loglike(LOW[,t],HIGH[,t],XB[,t],SHAPE,canSCALE,pcure)
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
          canll[,t]<-loglike(LOW[,t],HIGH[,t],canXB[,t],SHAPE,SCALE,pcure)
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
      keep.shape[i,]<-shape
      keep.scale[i,]<-scale
      keep.beta[i,]<-beta
      keep.pcure[i]<-pcure
      dev[i]<--2*sum(curll)

      if(i>burn){
        nnn<-iters-burn
        for(t in 1:nsubs){
          Ynew<-samp_censored(LOW[,t],HIGH[,t],XB[,t],SHAPE,SCALE,pcure)
          Y1[,t]<-Y1[,t]+ifelse(Ynew == Inf, 0, Ynew)
          YS[,t]<-YS[,t]+ifelse(Ynew == Inf, 0, 1)
          Y2[,t]<-Y2[,t]+ifelse(Ynew == Inf, 0, Ynew*Ynew)
        }
        invCPO<-invCPO+exp(-curll)/(iters-burn)
      }

      #DISPLAY CURRENT VALUE:
      if(i%%update==0 && FALSE){
         par(mfrow=c(3,3))
         image.plot(W,main="W",xlab="s",ylab="l")
         image.plot(factor,main="factor",xlab="s",ylab="l")
         plot(keep.scale[1:i,1],type="l",main="scale[1]")
         plot(keep.shape[1:i,1],type="l",main="shape[1]")
         plot(keep.beta[1:i,p],type="l",main="beta[p,1]")
         for(t in 1:3){if(t<=nsubs){
            plot(s[,1],LOW[,t],pch=ifelse(LOW[,t]==HIGH[,t],1,2),col=s[,2],ylim=c(0,200))
            SSS<-(XB[,t])/SHAPE
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
    for(t in 1:nsubs){
      curll[,t]<-loglike(LOW[,t],HIGH[,t],XB[,t],SHAPE,SCALE,pcure)
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
     range=keep.range,
     shape=keep.shape,
     scale=keep.scale,
     pcure=keep.pcure,
     beta=keep.beta,
     dev=dev,DIC=DIC,pD=pD,dbar=dbar,
     CPO=CPO,LPML=LPML)}

