######################################################
#0 load package
library(MASS)   #multiple random normal   
library(gtools) #dirichlet distribution
library(actuar) #rinvgamma random variable

#1 generate data
HMM.data=function(X,pi.0,pi.t,Beta.t,Sigma.t,information=T){
  n=nrow(X)
  k=length(pi.0)
  Y=matrix(NA,n,dim(pi.t)[3]+1)
  S=sample.int(k,n,replace = T,prob = pi.0)
  Y[,1]=sapply(1:n,function(x)rnorm(1,mean=sum(X[x,]*Beta.t[S[x],,1]),sd=sqrt(Sigma.t[S[x],1])))
  #Y[,1]=rnorm(n,mean=sum(X%*%Beta.t[S,,1]),sd=sqrt(Sigma.t[S,1]))  #faster
  S.old=S
  for(i in 1:dim(pi.t)[3]){
    S.temp=vector()
    for(j in 1:k){
      S.temp[S.old==j]=sample.int(k,sum(S.old==j),replace = T,prob = pi.t[j,,i])
    }
    S=cbind(S,S.temp)
    Y[,i+1]=sapply(1:n,function(x)rnorm(1,mean=sum(X[x,]*Beta.t[S.temp[x],,i+1]),sd=sqrt(Sigma.t[S.temp[x],i+1])))
    #Y[,i+1]=rnorm(n,mean=sum(X%*%Beta.t[S.temp,,i+1]),sd=sqrt(Sigma.t[S.temp,i+1]))   #faster
    S.old=S.temp
  }
  if(information){
    True.Value=list()  #True Value
    True.Value$beta=beta.t
    True.Value$sigma=sigma.t
    True.Value$init.distribution=Prob
    True.Value$transfer.matrix=Prob.t
    
    return(list(Y=Y,S=S,TrueValue=True.Value))
  }
  else
    return(Y)
}

HMM.data.missing=function(Y,S,beta_r,na=NA){
  nr=dim(Y)[1]
  nc=dim(Y)[2]
  Y.missing=Y
  R=array(dim=dim(Y))
  for(i in 1:nr)
    for(j in 1:nc){
      s=S[i,j]
      part=exp(beta_r[s,1]+beta_r[s,2]*Y[i,j])
      p0=1/(1+part)
      p1=1-p0
      R[i,j]=sample(c(0,1),1,prob=c(p0,p1))
      if(R[i,j]==0)
        Y.missing[i,j]=na
      #r=0 missing
    }
  return(list(Y.missing=Y.missing,R=R))
}

#test parameter to generate missing data
sz=700  #sample size
X=matrix(c(rep(1,sz),runif(sz,0,5)),sz,2)
beta=rbind(c(-1,0.5),c(1,-0.5))
beta.t=rep(beta,3)
dim(beta.t)=c(2,2,3)
sigma=c(0.3,0.3)
sigma.t=rep(sigma,3)
dim(sigma.t)=c(2,3)
Prob=c(0.4,0.6)
Prob.t=c(0.8,0.2,0.2,0.8,0.8,0.2,0.2,0.8)
dim(Prob.t)=c(2,2,2)

temp=HMM.data(X,Prob,Prob.t,beta.t,sigma.t,T)
Y=temp$Y
S=temp$S
beta_r=array(c(2,2,-1.5,-1.5),dim=c(2,2))
temp.missing=HMM.data.missing(Y,S,beta_r)
Y.missing=temp.missing$Y.missing
R=temp.missing$R
table(R)
R[Y > 2]
