dat <- read.csv("mcguire.csv")

n<-102

Low <-matrix(NA,32,n)
High<-matrix(NA,32,n)
HyPoor<-HyGood<-smoke<-diabetic<-hiv<-
        cr0<-probe0<-mobile0<-base_age<-matrix(0,32,n)

#Baseline crown-to-root ratio is cr0
#Baseline probing depth is probe0
#Baseline mobility is mobile0

ID<-dat[,1]
for(i in 1:n){

   d<-dat[ID==i,]

   HyPoor[,i]<-ifelse(d[1,5]=="P",1,0)
   HyGood[,i]<-ifelse(d[1,5]=="G",1,0)
   smoke[,i]<-d[1,101]
   diabetic[,i]<-d[1,102]
   hiv[,i]<-d[1,103]
   base_age[,i]<-age<-d[1,38]

   for(j in 1:32){
     dd<-d[d[,56]==j,] #column 56 is "tooth"
     if(nrow(dd)==1){
        dd<-as.numeric(dd)
        Low[j,i]<- ifelse(dd[169]==1,dd[168],dd[168])       #column 169 is "censor"
        High[j,i]<-ifelse(dd[169]==1,Inf,    dd[168]) #column 168 is "time"
        cr0[j,i]    <-dd[110] 
        probe0[j,i]  <-dd[28] 
        mobile0[j,i]<-dd[104] 
     }
   }
}

probe0<-ifelse(probe0<3,3,probe0)


rownames(Low)<-rownames(High)<-paste("tooth",1:32)
Low<-Low[-c(1,16,17,32),]
High<-High[-c(1,16,17,32),]

library(abind)
X<-abind(smoke, base_age,along=3)
X<-abind(X,HyPoor,along=3)
X<-abind(X,HyGood,along=3)
X<-abind(X,cr0,along=3)
X<-abind(X,probe0,along=3)
X<-abind(X,mobile0,along=3)

X<-X[-c(1,16,17,32),,]
rm(smoke,base_age,HyPoor,HyGood,cr0,probe0,mobile0)
nmiss<-apply(is.na(X),2,sum)
junk<-nmiss>10

X<-X[,!junk,]
Low<-Low[,!junk]
High<-High[,!junk]
n<-nrow(Low)

s1<-c(1:14,14:1)
s2<-c(rep(1,14),rep(2,14))

Xmiss1<-Xmiss2<-matrix(0,28,99)
for(j in 1:28){
  neigh<-abs(s1-s1[j])==1 & s2==s2[j]
  for(k in 1:99){
    Xmiss1[j,k]<-sum(is.na(High[neigh,k]))
  }
  neigh<-abs(s2-s2[j])==1 & s1==s1[j]
  for(k in 1:99){
    Xmiss2[j,k]<-sum(is.na(High[neigh,k]))
  }
}
X<-abind(X,Xmiss1,along=3)
X<-abind(X,Xmiss2,along=3)
rm(Xmiss1,Xmiss2,k,j,neigh)
dimnames(X)[[3]]<-c("Smoke","Age","HyPoor","HyGood","cr0","probe0","mobile0","Xmiss1","Xmiss2")

junk<-ls()
S<-cbind(s1,s2)

type<-rep(c(1,1,2,2,3,4,4,4,4,3,2,2,1,1),2)
type_labels<-c("Molar","Pre-Molar","Canine","Insicor")

rm(d,dd,dat,i,ID,diabetic,hiv,s1,s2,junk,age,nmiss)

save.image("NunnData.RData")
 