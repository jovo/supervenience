
library(MASS)
library(lpSolve)

getAk = function(n,p,m,q,labeled=T)
 ### generates adjacency matrix A for K(n,p,m,q)
 {
 A=matrix(sample(0:1,n^2,T,c(1-p,p)),nrow=n,ncol=n)
 if(labeled)  thism = 1:m
 else         thism = sample(1:n,m,F)
 A[thism,thism]=sample(0:1,m^2,T,c(1-q,q))
 A=A*upper.tri(A)
 A=A+t(A)
 return(A)
 }

Lstar = function( n=20,p=0.5,m=8,q=0.9 , q0vec ,pi0=0.5)
 ### class-conditional distributions are K(n,p,m,q0) and K(n,p,m,q)
 ### pi0 = 0.5  => class priors pi0=pi1=1/2
 ### Bayes optimal error rates L* ... assuming my math & coding are correct ?!?
 {
 pi1 = 1-pi0
 nc2 = choose(n,2)
 mc2 = choose(m,2)
 q1=q
 LstarL = LstarU = rep(0,length(q0vec))
 for(i in 1:length(q0vec))
  {
  q0=q0vec[i]
  for(e in 0:nc2)
   {
   for(eA in 0:min(mc2,e))
    {
    x = choose(mc2,eA)*choose(nc2-mc2,e-eA)
    thing0 = pi0*q0^eA*(1-q0)^(mc2-eA)*p^(e-eA)*(1-p)^(nc2-mc2-(e-eA))
    thing1 = pi1*q1^eA*(1-q1)^(mc2-eA)*p^(e-eA)*(1-p)^(nc2-mc2-(e-eA))
    minj = thing0 < thing1
    LstarL[i] = LstarL[i] + x*thing0*minj + x*(1-minj)*thing1
    }
   ttt0=ttt1=0
   for(eA in 0:min(mc2,e))
    {
    x = choose(mc2,eA)*choose(nc2-mc2,e-eA)
    ttt0 = ttt0 + x*pi0*q0^eA*(1-q0)^(mc2-eA)*p^(e-eA)*(1-p)^(nc2-mc2-(e-eA))
    ttt1 = ttt1 + x*pi1*q1^eA*(1-q1)^(mc2-eA)*p^(e-eA)*(1-p)^(nc2-mc2-(e-eA))
    }
   minj = ttt0 < ttt1
   LstarU[i] = LstarU[i] + (ttt0*minj + ttt1*(1-minj))
   }
  }
 return(list( LstarL , LstarU ))
 }

x = Lstar( q0vec = c(0.50,0.51,0.52,0.53,0.54,0.55,0.60,0.65,0.70,0.80,0.90) )
LstarL = x[[1]]
LstarU = x[[2]]




rung = function( q0=0.5,n=20,p=0.5,m=8,q=0.9,s=40,pi0=0.5,nmc=100,maxk=19 ,UUU=0)
# here we run monte carlos, with s training observations each,
# and generate various cross-val error rates Lhat ...
{
lcresL = lcresU = rep(1,nmc)
nnresLall = nnresLnndim = nnresUall = nnresUnndim = matrix(0,nrow=nmc,ncol=maxk)
for(mc in 1:nmc)
{
set.seed(mc)
s0 = rbinom(1,s,pi0)  ### s training observations
s1 = s-s0
A=list()  ### training data: first s0 are class0==K(n,p,m,q0) ; next s1 are class1==K(n,p,m,q)

### assuming ***labeled***
for(k in 1:s0)
 A[[k]]=getAk(n,p,m,q0)
for(k in (s0+1):(s0+s1))
 A[[k]]=getAk(n,p,m,q)

Delta=matrix(0,nrow=s,ncol=s)
for(i in 1:(s-1))
 for(j in (i+1):s)
  Delta[i,j] = sum((A[[i]]-A[[j]])^2)
Delta=Delta+t(Delta)

#plot(prcomp(Delta)$x[,1:2],pch=c(rep("*",s0),rep("+",s1)),col=c(rep("red",s0),rep("blue",s1)))

lcdim = 2
ldaclassif = lda(prcomp(Delta)$x[,1:lcdim],c(rep(0,s0),rep(1,s1)),CV=T)$class
lcresL[mc] = (sum(ldaclassif[1:s0]==0) + sum(ldaclassif[(s0+1):(s0+s1)]==1))/s

Deltap=Delta
for(k in 1:maxk)
 {
 knnresults = apply(Deltap,1,function(x){sum(order(x)[2:(k+1)]<=s0)>k/2})
 nnresLall[mc,k] = sum(knnresults[1:s0]) + sum(1-knnresults[(s0+1):(s0+s1)])
 }

nndim = 1
ppp=prcomp(Delta)$x[,nndim]
for(i in 1:s)
 for(j in 1:s)
  Deltap[i,j] = sqrt( ( ppp[i] - ppp[j] )^2 )
for(k in 1:maxk)
 {
 knnresults = apply(Deltap,1,function(x){sum(order(x)[2:(k+1)]<=s0)>k/2})
 nnresLnndim[mc,k] = sum(knnresults[1:s0]) + sum(1-knnresults[(s0+1):(s0+s1)])
 }

if(UUU){  # assuming ***UNlabeled***
for(k in 1:s0)
 A[[k]]=getAk(n,p,m,q0,labeled=F)
for(k in (s0+1):(s0+s1))
 A[[k]]=getAk(n,p,m,q,labeled=F)

Delta=matrix(0,nrow=s,ncol=s)
lpcost=A[[1]]
for(i in 1:(s-1))
 for(j in (i+1):s)
  {
  for(k in 1:n)
  for(l in 1:n)
   {
   lpcost[k,l] = sum( ( A[[i]][k,] - A[[j]][l,] )^2 )
   }
  lpsol = lp.assign(lpcost)$solution
  #numright[i,j] = sum(lpsol[1:m,1:m])

  Delta[i,j] = sum((lpsol%*%A[[i]]-A[[j]])^2)
  }
Delta=Delta+t(Delta)

#plot(prcomp(Delta)$x[,1:2],pch=c(rep("*",s0),rep("+",s1)),col=c(rep("red",s0),rep("blue",s1)))

ldaclassif = lda(prcomp(Delta)$x[,1:2],c(rep(0,s0),rep(1,s1)),CV=T)$class
#ldaclassif = lda(prcomp(Delta)$x[,1:(s0+s1-1)],c(rep(0,s0),rep(1,s1)),CV=T)$class
lcresU[mc] = (sum(ldaclassif[1:s0]==0) + sum(ldaclassif[(s0+1):(s0+s1)]==1))/s

Deltap=Delta
for(k in 1:maxk)
 {
 knnresults = apply(Deltap,1,function(x){sum(order(x)[2:(k+1)]<=s0)>k/2})
 nnresUall[mc,k] = sum(knnresults[1:s0]) + sum(1-knnresults[(s0+1):(s0+s1)])
 }

nndim = 1
ppp=prcomp(Delta)$x[,1]
for(i in 1:s)
 for(j in 1:s)
  Deltap[i,j] = sqrt( ( ppp[i] - ppp[j] )^2 )
for(k in 1:maxk)
 {
 knnresults = apply(Deltap,1,function(x){sum(order(x)[2:(k+1)]<=s0)>k/2})
 nnresUnndim[mc,k] = sum(knnresults[1:s0]) + sum(1-knnresults[(s0+1):(s0+s1)])
 }
}

}
return(list( lcresL , lcresU , nnresLall , nnresLnndim , nnresUall , nnresUnndim ))
}



 q0vec = c(0.50,0.51,0.52,0.53,0.54,0.55,0.60,0.65,0.70,0.80,0.90) 

date()
 res40 = res80 = list()
 for(i in 1:length(q0vec))
  {
  res40[[i]] = rung(q0=q0vec[i],nmc=100,s=40,UUU=1) # n=20,p=0.5,m=8,q=0.9,pi0=0.5, maxk=19
  res80[[i]] = rung(q0=q0vec[i],nmc=100,s=80,UUU=1) # n=20,p=0.5,m=8,q=0.9,pi0=0.5, maxk=19
  # length(q0)=11 , nmc=100 , s=40&80 , UUU=0 takes ~ 15 minutes (on my laptop)
  # length(q0)=11 , s=40&80 , s=40&80 , UUU=1 takes ~ 30 hours (on my laptop)
  }
date()

 res40meanL = res80meanL = rep(1,length(q0vec))
 res40meanU = res80meanU = rep(1,length(q0vec))
 for(i in 1:length(q0vec))
  {
  whicherror=4
  k=1
  # whicherror=3 , k=1 => LABELED 1nn in full dimensionality
  # whicherror=4 , k=1 => LABELED 1nn in 1st PC dim
  # whicherror=5 , k=1 => UNLABELED 1nn in full dimensionality
  # whicherror=6 , k=1 => UNLABELED 1nn in 1st PC dim
  res40meanL[i] = mean(1-res40[[i]][[whicherror]][,k]/40)
  res80meanL[i] = mean(1-res80[[i]][[whicherror]][,k]/80)
  whicherror=6
  k=1
  # whicherror=5 , k=1 => UNLABELED 1nn in full dimensionality
  # whicherror=6 , k=1 => UNLABELED 1nn in 1st PC dim
  res40meanU[i] = mean(1-res40[[i]][[whicherror]][,k]/40)
  res80meanU[i] = mean(1-res80[[i]][[whicherror]][,k]/80)
  }

#####postscript("supervenienceres1.eps",width=8)
plot(q0vec,res40meanL,xlab="q0",ylab="L",ylim=c(0,.5),type="b",pch="40",col=2)
points(q0vec,res80meanL,type="b",pch="80",col=2)
points(q0vec,LstarL,type="b",pch="*",col=2)
points(q0vec,res40meanU,type="b",pch="40",col=3)
points(q0vec,res80meanU,type="b",pch="80",col=3)
points(q0vec,LstarU,type="b",pch="*",col=3)
#####dev.off()

whichq=6
whicherror=4
# whichq=6 , whicherror=4 => LABELED k-nrstnbrs in 1st PC dim for q0=0.55
for(k in 1:19)
 print(mean(1-res80[[whichq]][[whicherror]][,k]/80))
# compare:  LstarL[[6]] = 0.01325911
