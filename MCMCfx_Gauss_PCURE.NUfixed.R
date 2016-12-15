
loglike<-function(LOW,HIGH,theta,Xb,shape,scale,pcure){
  ## SSS<-(log(theta)+Xb)/shape
  
  SSS<-(theta+Xb)/shape
  SSS<-exp(-SSS)*scale
  
  # SSS <- ifelse(SSS <=0, 0.0001, SSS)
  # SSS <- ifelse(SSS == "NaN", 0.0001, SSS)
  
  D<-dweibull(LOW,shape,SSS)
  L<-pweibull(LOW,shape,SSS)
  H<-pweibull(HIGH,shape,SSS)
  LL<-log(ifelse(LOW==HIGH, (1-pcure)*D, pcure+(1-pcure)*(H-L)))
  return(LL)
}




samp_censored<-function(LOW,HIGH,theta,Xb,shape,scale,pcure){
  
  ## SSS<-(log(theta)+Xb)/shape
  
  SSS<-(theta+Xb)/shape
  SSS<-exp(-SSS)*scale
  
  # SSS <- ifelse(SSS <=0, 0.0001, SSS)
  # SSS <- ifelse(SSS == "NaN", 0.0001, SSS)
  
  L<-pweibull(LOW,shape,SSS)
  H<-pweibull(HIGH,shape,SSS)
  cured<-sample(c(0,Inf), length(L), replace=TRUE, prob=c(1-pcure, pcure))
  U<-runif(length(L),L,H)
  Y<-qweibull(U,shape,SSS)+cured
  return(Y)
}



inv<-function(X,thresh=exp(-10)){
  E<-eigen(X)
  V<-E$vectors
  D<-E$values
  D<-ifelse(D<thresh,thresh,D)
  PREC<-V%*%diag(1/D)%*%t(V)
  LOGDET<-sum(-log(D))
list(PREC=PREC,LOGDET=LOGDET)}




# Y = [numbers sites] x [nsubs] (Follow-up time)
# C = [numbers sites] x [nsubs] (TRUE if censored, FALSE ow)
# X = [nsubs] x p
# s = [n sites] x 2
# knots = [n knots] x 2  

SpatSurvGauss<-function(LOW,HIGH,X,s,type,
    common_scale=F,common_shape=T,common_beta=T,
    iters=5000,burn=1000,update=10,nthin=5,nblocks=NULL){

    library(fields)
    library(emulator)
    library(geoR)


    #BOOKKEEPING
    ns<-nrow(s)
    nsub<-nsubs<-ncol(LOW)
    L<-nrow(knots)
    p<-dim(X)[3]
    ntype<-max(type)

    D1<-rdist(s[,1],s[,1])^2
    D2<-rdist(s[,2],s[,2])^2
    D1[D1<0.00001]<-0
    D2[D2<0.00001]<-0
    
    
    
    #INITIAL VALUES:
    beta<-rep(0,p)
    scale<-rep(max(LOW),ntype)
    shape<-rep(1,ntype)
    range<-rep(1,2)
    nu<- 0.5
    tau<-1
    pcure<-0.15

    E<-  inv(matern(sqrt(range[1]*D1+range[2]*D2),0.5,nu))
    
    ## print(E)
    

    theta<-matrix(0,ns,nsubs)
    XB<-matrix(0,ns,nsubs)
    for(j in 1:nsubs){XB[,j]<-X[,j,]%*%beta}
    SCALE<-scale[type]
    SHAPE<-shape[type]

    curll<-LOW
    for(t in 1:nsubs){
      curll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,SCALE, pcure)
    }
      
    ## keep.nu<-
    keep.sd<-rep(0,iters)
    keep.range<-matrix(0,iters,2)
    keep.shape<-keep.scale<-matrix(0,iters,ntype)
    keep.beta<-matrix(0,iters,p)
    keep.pcure<-rep(0,iters)

    colnames(keep.beta)<-dimnames(X)[[3]]
    att<-acc<-MH<-rep(0.05,p+10)
    MH[p+5]<-1
    dev<-rep(0,iters)
    Y1<-Y2<-Y3<-Y4<-invCPO<-matrix(0,ns,nsubs)


    SSS<-s;SSS[,2]<-100*SSS[,2]

    if(!is.null(nblocks)){
       block<-kmeans(SSS,centers=nblocks)$cluster
    }
    if(is.null(nblocks)){
      nblocks<-ns
      block<-1:ns
    }
  
    for(i in 1:iters){

     for(ttt in 1:nthin){
      ##########################################################
      ##############      Random effects theta    ############
      ##########################################################
      QQQ<-rep(0,nsubs)
      for(t in 1:nsubs){
        QQQ[t]<--0.5*tau*quad.form(E$PREC,theta[,t])
        for(l in 1:nblocks){
          bbb<-block==l
          nnn<-sum(bbb)
          att[p+5]<-att[p+5]+1
          cantheta<-theta[,t]
          cantheta[bbb]<-theta[bbb,t]+
                         MH[p+5]*(rnorm(1)+rnorm(nnn))/2
          canQQQ<--0.5*tau*quad.form(E$PREC,cantheta)
          canll<-loglike(LOW[bbb,t],HIGH[bbb,t],cantheta[bbb],XB[bbb,t],SHAPE[bbb],SCALE[bbb], pcure[bbb])
          R<-sum(canll-curll[bbb,t])+canQQQ-QQQ[t]
          if(!is.na(exp(R))){if(runif(1)<exp(R)){
            theta[,t]<-cantheta
            curll[bbb,t]<-canll
            QQQ[t]<-canQQQ
            acc[p+5]<-acc[p+5]+1
          }}
         }
       }


      ##########################################################
      #########     spatial hyperparameters             ########
      ##########################################################

      #Precision
      QQQ<-sum(QQQ)
      if(i > 50){
        QQQ<- -2*QQQ/tau
        tau<-rgamma(1,ns*nsubs/2+.1,QQQ/2+.1)
        QQQ<- -0.5*tau*QQQ
      }

      #Range
      for(j in 1:2){
        att[p+4]<-att[p+4]+1
        canrange<-range
        canrange[j]<-exp(rnorm(1,log(range[j]),MH[p+4]))
        canE<-  inv(matern(sqrt(canrange[1]*D1+canrange[2]*D2),1,nu))
        
     
        canQQQ<-0
        for(t in 1:nsubs){
          canQQQ<-canQQQ-0.5*tau*quad.form(canE$PREC,theta[,t])
        }
        R<-0.5*nsubs*(canE$LOGDET-E$LOGDET)+
           canQQQ-QQQ+
           dnorm(log(canrange[j]),0,1,log=T)-
           dnorm(log(range[j]),0,1,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
           range<-canrange;QQQ<-canQQQ;E<-canE
           acc[p+4]<-acc[p+4]+1
        }}           
      }
      

      #Smoothness
      #att[p+3]<-att[p+3]+1
      ## cannu<-exp(rnorm(1,log(nu),MH[p+3]))
      #canE<- inv(matern(sqrt(range[1]*D1+range[2]*D2),1,nu))
      
      
      # canQQQ<-0
      #for(t in 1:nsubs){
      #  canQQQ<-canQQQ-0.5*tau*quad.form(canE$PREC,theta[,t])
      #}
      #R<-0.5*nsubs*(canE$LOGDET-E$LOGDET)+
      #   canQQQ-QQQ+
      #   dnorm(log(nu),log(.5),1,log=T)-
      #   dnorm(log(nu),log(.5),1,log=T)
         
      #if(!is.na(exp(R))){if(runif(1)<exp(R)){
      #   QQQ<-canQQQ;E<-canE
      #   acc[p+3]<-acc[p+3]+1
      #}}           
      
      
      
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
        att[1]<-att[1]+1
        canshape<-shape
        canshape[j]<-exp(rnorm(1,log(shape[j]),MH[1]))
        if(common_shape){canshape[1:ntype]<-canshape[1]}
        canSHAPE<-canshape[type]
        canll<-curll
        for(t in 1:nsubs){
          canll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],canSHAPE,SCALE, pcure)
        }
        R<-sum(canll-curll)+
           dnorm(log(canshape[j]),0,10,log=T)-
           dnorm(log(shape[j]),0,10,log=T)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
           shape<-canshape;SHAPE<-canSHAPE;curll<-canll;
           acc[1]<-acc[1]+1
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
          canll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,canSCALE, pcure)
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
          canll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],canXB[,t],SHAPE,SCALE, pcure)
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

      for(j in 1:length(MH)){if(i<burn/2 & att[j]>25){
        if(acc[j]/att[j]<0.3){MH[j]<-MH[j]*0.8}
        if(acc[j]/att[j]>0.6){MH[j]<-MH[j]*1.2}
        acc[j]<-att[j]<-0
      }}


      #KEEP TRACK OF STUFF:
      keep.sd[i]<-1/sqrt(tau)
      keep.range[i,]<-1/sqrt(range)
      #keep.nu[i]<-nu
      keep.shape[i,]<-shape
      keep.scale[i,]<-scale
      keep.beta[i,]<-beta
      dev[i]<--2*sum(curll)
      keep.pcure[i]<-pcure

      if(i>burn){
        nnn<-iters-burn
        for(t in 1:nsubs){
          Ynew<-samp_censored(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,SCALE, pcure)
          Y1[,t]<-Y1[,t]+Ynew/nnn
          Y2[,t]<-Y2[,t]+Ynew*Ynew/nnn
        }
        Y3<-Y3+theta/nnn
        Y4<-Y4+theta*theta/nnn
        invCPO<-invCPO+exp(-curll)/(iters-burn)
      }

      #DISPLAY CURRENT VALUE:
      if(i%%update==0 && FALSE){
         par(mfrow=c(3,3))
         plot(keep.range[1:i,1],type="l",main="range1")
         plot(keep.range[1:i,2],type="l",main="range2")
#         plot(keep.nu[1:i],type="l",main="nu")
         plot(keep.sd[1:i],type="l",main="sd")
         plot(keep.scale[1:i,1],type="l",main="scale[1]")
         plot(keep.shape[1:i,1],type="l",main="shape[1]")
         plot(keep.beta[1:i,p],type="l",main="beta[p]")
         for(t in 1:3){if(t<=nsubs){
            plot(s[,1],LOW[,t],pch=ifelse(LOW[,t]==HIGH[,t],1,2),col=s[,2],ylim=c(0,200))
            SSS<-(theta[,t]+XB[,t])/SHAPE
            SSS<-exp(-SSS)*SCALE
            lines(s[,1],SSS*(log(2))^(1/SHAPE),col=s[,2])
         }}
      }


    }

    if(ntype==1){
     shape<-mean(keep.shape[burn:iters,])[type]
     scale<-mean(keep.scale[burn:iters,])[type]
    }
    if(ntype>1){
     shape<-colMeans(keep.shape[burn:iters,])[type]
     scale<-colMeans(keep.scale[burn:iters,])[type]
    }
    beta<-colMeans(keep.beta[burn:iters,])
    for(j in 1:nsubs){XB[,j]<-X[,j,]%*%beta}
    theta<-Y3
    for(t in 1:nsubs){
      curll[,t]<-loglike(LOW[,t],HIGH[,t],theta[,t],XB[,t],SHAPE,SCALE, pcure)
    }
    dhat<--2*sum(curll)
    dbar<-mean(dev[burn:iters])
    pD<-dbar-dhat
    DIC<-dbar+pD

    CPO<-1/invCPO
    LPML<-sum(log(CPO))


list(pred.y.mn=Y1,
     pred.y.var=Y2-Y1^2,
     theta.mn=Y3,
     theta.var=Y4-Y3^3,
     sd=keep.sd,
     range=keep.range,
     #nu=keep.nu,
     shape=keep.shape,
     scale=keep.scale,
     beta=keep.beta,
     pcure=keep.pcure,
     dev=dev,DIC=DIC,pD=pD,dbar=dbar,
     CPO=CPO,LPML=LPML)}





