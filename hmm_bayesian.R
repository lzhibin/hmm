######################################################
#0 load package
library(MASS)   #multiple random normal   
library(gtools) #dirichlet distribution
library(actuar) #rinvgamma random variable

#1 generate data
HMM.data=function(X,pi.0,pi.t,Beta.t,Sigma.t,loc=F){
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
  if(loc)
    return(list(Y=Y,S=S))
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
  
   # D=Prob #fix initial distribution and transfer matrix  ###code for test
   # Dt=Prob.t
   # dim(Dt)=c(2,2,2)
    
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

    
    #FFBS renew hidden variable
    # P.S.pred=array(dim=c(NY,p,Ti))
    # P.S.Update=array(dim=c(NY,p,Ti))
    # P.S.pred[,,1]=matrix(rep(D,NY),NY,p,byrow = T)
    # for(i in 1:NY){
    #   for(j in 1:p){
    #     P.S.Update[i,j,1]=P.S.pred[i,j,1]*dnorm(Y[i,1],sum(X[i,]*betak[j,,1]),sqrt(sigmak[j,1]))
    #   }
    # }
    # for(i in 1:NY){
    # }
    # #T=1ʱ
    # for(t in 2:Ti){
    #   P.S.pred[,,t]=P.S.Update[,,t-1]%*%Dt[,,t-1]
    #   for(i in 1:NY){
    #     for(j in 1:p){
    #       P.S.Update[i,j,t]=P.S.pred[i,j,t]*dnorm(Y[i,t],sum(X[i,]*betak[j,,t]),sqrt(sigmak[j,t]))
    #     }
    #   }
    # }
    # #T>1ʱ
    
    # #renew S
    # Post.S=array(dim = c(NY,p,Ti))
    # for(i in 1:NY){
    #   S[i,Ti]=sample.int(p,1,prob = P.S.Update[i,,Ti])
    #   for(t in (Ti-1):1){
    #     tem=vector()
    #     for(j in 1:p){
    #       tem[j]=P.S.Update[i,j,t]*Dt[j,S[i,t+1],t]
    #     }
    #     S[i,t]=sample.int(p,1,prob = tem)
    #   }
    # }
    # 
    
    #FFBS renew hidden variable altnative code
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


    #print(betak)
    #print(sigmak)
    #print(D)
    #print(Dt)
    
    #store parameters
    p.betak[,,,i]=betak
    p.sigmak[,,i]=sigmak
    p.D[,i]=D
    p.Dt[,,,i]=Dt
   }
   ParaMeters=list(p.betak=p.betak,p.sigmak=p.sigmak,p.D=p.D,p.Dt=p.Dt)
   class(ParaMeters)="hmm"
   return(ParaMeters)
}

#method to summary hmm class
summary.hmm=function(data,n){
    len=dim(data$p.D)[2]
    if(len<n){
      print("error! need more data")
    }

    beta=apply(data$p.betak[,,,(len-n+1):len],c(1,2,3),mean)
    sigma=apply(data$p.sigmak[,,(len-n+1):len],c(1,2),mean)
    init.distribution=apply(data$p.D[,(len-n+1):len],1,mean)
    transfer.matrix=apply(data$p.Dt[,,,(len-n+1):len],c(1,2,3),mean)
    summary.hmm=list(beta=beta,sigma=sigma,
                     init.distribution=init.distribution,
                     transfer.matrix=transfer.matrix)
    return(summary.hmm)
}

#test result
X=matrix(c(rep(1,1000),runif(1000,0,5)),1000,2)
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

b0=solve(t(X)%*%X)%*%t(X)%*%Y[,1]
B0=1.5*diag(2)
c0=1.28
C0=0.36*var(as.vector(Y))
E0=c(1,1)
Et=Prob.t*2
m=1500
system.time(g<-gibbs.hmm(Y,X,b0,B0,c0,C0,E0,Et,S,m,rep=1))
summary(g,m-500)
