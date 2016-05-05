h######################################################
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



#2 MCMC
gibbs.hmm=function(Y,X,b0,B0,c0,C0,E0,Et,S,m,rep=1){
  #initial value
  d=ncol(X)    #numbers of beta
  p=length(E0) #levels of hidden variable 
  Ti=ncol(Y)   #numbers of periods
  NY=nrow(Y)
  B0.inv=solve(B0)
  Bk=array(dim=c(d,d,p,Ti))
  bk=array(dim=c(p,d,Ti))

  ck=array(dim=c(p,Ti))
  Ck=array(dim=c(p,Ti))
  betak=array(dim=c(p,d,Ti))
  sigmak=array(dim=c(p,Ti))
  Dt=array(dim=c(p,p,Ti-1)) #transfer matrix [f(u_j|u_hat_i)]
  NS=array(dim=c(p,Ti))
  NSt=array(dim=c(p,p,Ti-1))

  p.betak=array(dim=c(p,d,Ti,m))
  p.sigmak=array(dim=c(p,Ti,m))
  p.D=array(dim=c(p,m))
  p.Dt=array(dim=c(p,p,Ti-1,m))

  
  temp.beta=mvrnorm(p,b0,B0)
  temp.sigma=rinvgamma(p,c0,scale=C0)
  for(t in 1:Ti){
    betak[,,t]=temp.beta
    sigmak[,t]=temp.sigma
  }

  
   for(i in 1:m){
    for(t in 1:Ti){
      NS[,t]=sapply(1:p,function(x)sum(S[,t]==x))
    }
    D=rdirichlet(1,E0+NS[,1])
    #generate init distribution

    for(t in 1:(Ti-1)){
      for(k in 1:p){
        NSt[k,,t]=sapply(1:p,function(x)sum(S[S[,t]==k,t+1]==x))
        Dt[k,,t]=rdirichlet(1,Et[k,,t]+NSt[k,,t])
      }
    }
    #generate transfer matrix
      
    for(t in 1:Ti){
      for(j in 1:p){
        Bk[,,j,t]=solve(B0.inv+t(X[S[,t]==j,])%*%X[S[,t]==j,]/sigmak[j,t])
        bk[j,,t]=Bk[,,j,t]%*%(B0.inv%*%b0+t(X[S[,t]==j,])%*%Y[S[,t]==j,t]/sigmak[j,t])
        
        betak[j,,t]=mvrnorm(1,bk[j,,t],Bk[,,j,t])
        
        ck[j,t]=c0+NS[j,t]/2
        Ck[j,t]=C0+sum((Y[S[,t]==j,t]-X[S[,t]==j,]%*%betak[j,,t])^2)/2
        
        sigmak[j,t]=rinvgamma(1,ck[j,t],scale=Ck[j,t])
        
      }
    }
    #renew beta and sigma

      f.density=function(y,x,beta,sigma){
      m=as.vector(beta%*%x)
      s=sqrt(sigma)
      res=dnorm(y,m,s)
      return(res)
    }
    P.S.Update=array(dim=c(Ti,p,NY))
    for(nr in 1:NY){
      P.S.Update[1,,nr]=D*f.density(Y[nr,1],X[nr,],betak[,,1],sigmak[,1])
      for(t in 2:Ti){
        P.S.Update[t,,nr]=t(P.S.Update[t-1,,nr])%*%Dt[,,t-1]*f.density(Y[nr,t],X[nr,],betak[,,t],sigmak[,t])
      }
    }
    for(nr in 1:NY){
      S[nr,Ti]=sample.int(p,1,prob = P.S.Update[Ti,,nr])
      for(t in (Ti-1):1){
        s.temp=P.S.Update[t,,nr]*Dt[,S[nr,t+1],t]
        S[nr,t]=sample.int(p,1,prob = s.temp)
      }
    }
    #FFBS renew lacation parameter


    #label switch
    base=betak[,1,]
    for(t in 1:Ti){
      betak[,,t]=betak[order(base[,t]),,t]
      sigmak[,t]=sigmak[order(base[,t]),t]
      S[,t]=rank(base[,t])[S[,t]]
    }
    D=D[order(base[,1])]
    for(t in 1:(Ti-1)){
      #switch Dt
      Dt[,,t]=Dt[order(base[,t]),order(base[,t+1]),t]

    }

    
    #store parameters
    p.betak[,,,i]=betak
    p.sigmak[,,i]=sigmak
    p.D[,i]=D
    p.Dt[,,,i]=Dt
   }
   ParaMeters=list(p.betak=p.betak,p.sigmak=p.sigmak,p.D=p.D,p.Dt=p.Dt)
   class(ParaMeters)="hmm"
   if(rep==1)
     return(ParaMeters)
   else{
     groups=list()
     groups[[1]]=ParaMeters
     for(r in 2:rep)
       groups[[r]]=gibbs.hmm(Y,X,b0,B0,c0,C0,E0,Et,S,m,rep=1)
     class(groups)='hmmgroups'
   }
   return(groups)
}

#method to summary hmm class
summary.hmm=function(data,n,TrueValue=NULL){
    len=dim(data$p.D)[2]
    if(len<n){
      print("error! need more data")
    }

    beta=apply(data$p.betak[,,,(len-n+1):len],c(1,2,3),mean)
    beta.sd=apply(data$p.betak[,,,(len-n+1):len],c(1,2,3),sd)
    sigma=apply(data$p.sigmak[,,(len-n+1):len],c(1,2),mean)
    sigma.sd=apply(data$p.sigmak[,,(len-n+1):len],c(1,2),sd)
    init.distribution=apply(data$p.D[,(len-n+1):len],1,mean)
    init.distribution.sd=apply(data$p.D[,(len-n+1):len],1,sd)
    transfer.matrix=apply(data$p.Dt[,,,(len-n+1):len],c(1,2,3),mean)
    transfer.matrix.sd=apply(data$p.Dt[,,,(len-n+1):len],c(1,2,3),sd)
    
    summary.hmm=list(beta=beta,beta.sd=beta.sd,sigma=sigma,sigma.sd=sigma.sd,
                     init.distribution=init.distribution,
                     init.distribution.sd=init.distribution.sd,
                     transfer.matrix=transfer.matrix,
                     transfer.matrix.sd=transfer.matrix.sd)
    if(! is.null(TrueValue))
        summary.hmm=list(beta.TrueValue=TrueValue$beta,
                         beta=beta,beta.sd=beta.sd,
                         sigma.TrueValue=TrueValue$sigma,
                         sigma=sigma,sigma.sd=sigma.sd,
                         init.distribution..TrueValue=TrueValue$init.distribution,
                         init.distribution=init.distribution,
                         init.distribution.sd=init.distribution.sd,
                         transfer.matrix.TrueValue=TrueValue$transfer.matrix,
                         transfer.matrix=transfer.matrix,
                         transfer.matrix.sd=transfer.matrix.sd)
    return(summary.hmm)
}

summary.hmmgroups=function(data,n,TrueValue=NULL){
  r=length(data)
  Sample=summary(data[[1]],n)
  
  beta=array(dim=c(dim(Sample$beta),r))
  beta.sd=array(dim=c(dim(Sample$beta.sd),r))
  sigma=array(dim=c(dim(Sample$sigma),r))
  sigma.sd=array(dim=c(dim(Sample$sigma.sd),r))
  init.distribution=array(dim=c(length(Sample$init.distribution),r))
  init.distribution.sd=array(dim=c(length(Sample$init.distribution.sd),r))
  transfer.matrix=array(dim=c(dim(Sample$transfer.matrix),r))
  transfer.matrix.sd=array(dim=c(dim(Sample$transfer.matrix.sd),r))
  
  for(i in 1:r){
    tem=summary(data[[i]],n)
    
    beta[,,,i]=tem$beta
    beta.sd[,,,i]=tem$beta.sd
    sigma[,,i]=tem$sigma
    sigma.sd[,,i]=tem$sigma.sd
    init.distribution[,i]=tem$init.distribution
    init.distribution.sd[,i]=tem$init.distribution.sd
    transfer.matrix[,,,i]=tem$transfer.matrix
    transfer.matrix.sd[,,,i]=tem$transfer.matrix.sd
  }
  if(is.null(TrueValue))
    res=list(
      beta=apply(beta,1:3,mean),beta.sd=apply(beta.sd,1:3,mean),
      sigma=apply(sigma,1:2,mean),sigma.sd=apply(sigma.sd,1:2,mean),
      init.distribution=apply(init.distribution,1,mean),
      init.distribution.sd=apply(init.distribution.sd,1,mean),
      transfer.matrix=apply(transfer.matrix,1:3,mean),
      transfer.matrix.sd=apply(transfer.matrix.sd,1:3,mean)
    )
  else
    res=list(
      beta.TrueValue=TrueValue$beta,
      beta.estimate=apply(beta,1:3,mean),
      beta.sd=apply(beta.sd,1:3,mean),
      
      sigma.TrueValue=TrueValue$sigma,
      sigma.estimate=apply(sigma,1:2,mean),
      sigma.sd=apply(sigma.sd,1:2,mean),
      
      init.distribution.TrueValue=TrueValue$init.distribution,
      init.distribution.estimate=apply(init.distribution,1,mean),
      init.distribution.sd=apply(init.distribution.sd,1,mean),
      
      transfer.matrix.TrueValue=TrueValue$transfer.matrix,
      transfer.matrix.estimate=apply(transfer.matrix,1:3,mean),
      transfer.matrix.sd=apply(transfer.matrix.sd,1:3,mean)
    )
  
  return(res)
}

#test result
sz=500  #sample size
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
TV=temp$TrueValue  #True value
S=HMM.data(X,Prob,Prob.t,beta.t,sigma.t,T)$S
#generate half random S
S=sample.int(length(Prob),nrow(Y)*ncol(Y),replace = T)
dim(S)=dim(Y)
#generate completely random S

b0=solve(t(X)%*%X)%*%t(X)%*%Y[,1]
B0=1.5*diag(2)
c0=1.28
C0=0.36*var(as.vector(Y))
E0=c(1,1)
Et=Prob.t*2
m=700
system.time(g<-gibbs.hmm(Y,X,b0,B0,c0,C0,E0,Et,S,m,rep=1))
Sys.time()
summary(g,m-500,TV)

