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


#2 MCMC
mcmc.hmm=function(Y,R,X,b0,B0,c0,C0,E0,Et,S,m,beta0_r=NULL,lambda_y=1,lambda_r=1,rep=1){
  # accept_y=vector()
  # accept_r=vector()
  # count_y=1
  # count_r=1

  #initial value
  d=ncol(X)    #numbers of beta
  p=length(E0) #levels of hidden variable 
  Ti=ncol(Y)   #numbers of periods
  NY=nrow(Y)   #numbers of samples
  B0.inv=solve(B0)
  if(is.null(beta0_r))
    beta0_r=array(0,dim=c(p,2))
  
  #variables for every step
  Bk=array(dim=c(d,d,p,Ti))
  bk=array(dim=c(p,d,Ti))
  ck=array(dim=c(p,Ti))
  Ck=array(dim=c(p,Ti))
  
  betak=array(dim=c(p,d,Ti))
  sigmak=array(dim=c(p,Ti))
  betak_r=beta0_r
  NS=array(dim=c(p,Ti))
  NSt=array(dim=c(p,p,Ti-1))
  Dt=array(dim=c(p,p,Ti-1)) #transfer matrix [f(u_j|u_hat_i)]
  
  #variables to save parameters
  p.betak=array(dim=c(p,d,Ti,m))
  p.sigmak=array(dim=c(p,Ti,m))
  p.D=array(dim=c(p,m))
  p.Dt=array(dim=c(p,p,Ti-1,m))
  p.betak_r=array(dim=c(p,2,m))
  
  #init betak and sigmak for hmm
  temp.beta=mvrnorm(p,b0,B0)
  temp.sigma=rinvgamma(p,c0,scale=C0)
  for(t in 1:Ti){
    betak[,,t]=temp.beta
    sigmak[,t]=temp.sigma
  }
  
  #init missing Y base on prior distribution
  for(t in 1:Ti)
    for(ny in 1:NY){
      if(is.na(Y[ny,t]))
        Y[ny,t]=rnorm(1,sum(X[ny,]*betak[S[ny,t],,t]))
    }
  
  #density function for FFBS to renew S
  f.density=function(y,x,beta,sigma){
    m=as.vector(beta%*%x)
    s=sqrt(sigma)
    res=dnorm(y,m,s)
    return(res)
  }
  
  #likelihood ratio for beta_r_new / beta_r_old
  ln_beta_r_ratio=function(r,y,beta1,beta0){
    #beta1 and beta0 are the new value and old value of beta_r
    temp1=sum(r)*(beta1[1]-beta0[1])+sum(r*y)*(beta1[2]-beta0[2])
    temp2=sum(log(1+exp(beta0[1]+beta0[2]*y)))-sum(log(1+exp(beta1[1]+beta1[2]*y)))
    temp=exp(temp1+temp2)
    return(temp)
    #exp to diff of log likelihood
  }
  
  #information matrix estimator for beta_r ;return inverse
  im_beta_r=function(y,beta=NULL){
    #the beta here should be the estimator of beta_r ;using zero vector as default
    temp=array(dim=c(2,2))
    if(is.null(beta)){
      temp[1,1]=length(y)
      temp[1,2]=sum(y)
      temp[2,1]=temp[1,2]
      temp[2,2]=sum(y^2)
      temp=temp/4
      #it means beta=c(0,0)
    }
    else{
      d1=exp(beta[1]+beta[2]*y)
      d=d1/((1+d1)^2)
      temp[1,1]=sum(d)
      temp[1,2]=sum(y*d)
      temp[2,1]=temp[1,2]
      temp[2,2]=sum(d*(y^2))
    }
    return(solve(temp))
  }
  
  #likelihood ratio for yi_new / yi_old
  ln_yi_ratio=function(beta,x,sigma,y1,y0){
    #the beta here is betak[k,,t];sigma is sigmak[k,t] y1 and y0 is the new value and old value
    temp1=((y0-sum(x*beta))^2-(y1-sum(x*beta))^2)/(2*sigma) #+beta[2]*r*(y1-y0)=0 for r=0
    temp2=log(1+exp(beta[1]+beta[2]*y0))-log(1+exp(beta[1]+beta[2]*y1))
    temp=exp(temp1+temp2)
    return(temp)
    #take exp to diff of log likelihood
  }
  
  #information matrix estimator for yi ;return inverse ;in this sample it is number
  im_yi=function(beta,sigma,y=0){
    d1=(beta[2]^2)
    d2=exp(beta[1]+beta[2]*y)
    temp=(d1*d2)/((1+d2)^2)+1/sigma
    return(1/temp)
  }
  
  
  #main iteration
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
    
    
    #FFBS renew hidden variable altnative code
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
    
    #sample beta_r 
    for(k in 1:p){
      yk=Y[S==k]
      rk=R[S==k]
      temp=mvrnorm(1,betak_r[k,],lambda_r*im_beta_r(yk))  #lambda_r to control accept rate
      alpha=min(1,ln_beta_r_ratio(rk,yk,temp,betak_r[k,]))
      #cat("acceptrate_beta:",alpha,"\n")  #to see the accept rate
      U01=runif(1)
      
      # accept_r[count_r]=0 #calc accept rate for beta
      
      if(U01<alpha){
        betak_r[k,]=temp
        # accept_r[count_r]=1
      }
      # count_r=count_r+1 
    }
    
    #sample missing yi
    for(ny in 1:NY)
      for(t in 1:Ti)
        if(R[ny,t]==0){
          nyt.s=S[ny,t]
          temp=rnorm(1,Y[ny,t],sqrt(lambda_y*im_yi(betak_r[nyt.s,],sigmak[nyt.s,t])))
          alpha=min(1,ln_yi_ratio(beta=betak[nyt.s,,t],x=X[ny,],sigma=sigmak[nyt.s,t],y1=temp,y0=Y[ny,t]))
          #cat("acceptrate_yi:",alpha,"\n")  #to see the accept rate
          U01=runif(1)
          
          # accept_y[count_y]=0 #calc accept rate for y
          
          if(U01<alpha){
            Y[ny,t]=temp
            # accept_y[count_y]=1
          }
          # count_y=count_y+1
        }
    
    
    #label switch
    base=betak[,1,]
    for(t in 1:Ti){
      betak[,,t]=betak[order(base[,t]),,t]
      sigmak[,t]=sigmak[order(base[,t]),t]
      S[,t]=rank(base[,t])[S[,t]]
    }
    
    betak_r=betak_r[order(base[,1]),]
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
    p.betak_r[,,i]=betak_r
  }
  ParaMeters=list(p.betak=p.betak,p.sigmak=p.sigmak,p.D=p.D,p.Dt=p.Dt,p.betak_r=p.betak_r)
  class(ParaMeters)="hmm.missing"
  if(rep==1){
    # cat("beta accept rate:",sum(accept_r)/length(accept_r),"\n")
    # cat("Y accept rate:",sum(accept_y)/length(accept_y),"\n")
    return(ParaMeters)
  }
  else{
    groups=list()
    groups[[1]]=ParaMeters
    for(r in 2:rep)
      groups[[r]]=mcmc.hmm(Y,R,X,b0,B0,c0,C0,E0,Et,S,m,beta0_r,lambda_y,lambda_r,rep=1)
    class(groups)='hmmgroups.missing'
  }
  # cat("beta accept rate:",sum(accept_r)/length(accept_r),"\n")
  # cat("Y accept rate:",sum(accept_y)/length(accept_y),"\n")
  return(groups)
}

#method to summary hmm.missing class
summary.hmm.missing=function(data,n,TrueValue=NULL){
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
  beta_r=apply(data$p.betak_r[,,(len-n+1):len],c(1,2),mean)
  beta_r.sd=apply(data$p.betak_r[,,(len-n+1):len],c(1,2),sd)
  if(is.null(TrueValue))
    summary.hmm=list(beta=beta,beta.sd=beta.sd,sigma=sigma,sigma.sd=sigma.sd,
                     init.distribution=init.distribution,
                     init.distribution.sd=init.distribution.sd,
                     transfer.matrix=transfer.matrix,
                     transfer.matrix.sd=transfer.matrix.sd,
                     beta_r=beta_r,
                     beta_r.sd=beta_r.sd)
  else
    summary.hmm=list(beta.True=TrueValue$beta,beta=beta,beta.sd=beta.sd,
                     sigama.True=TrueValue$sigma,sigma=sigma,sigma.sd=sigma.sd,
                     init.distribution.True=TrueValue$init.distribution,
                     init.distribution=init.distribution,
                     init.distribution.sd=init.distribution.sd,
                     transfer.matrix.True=TrueValue$transfer.matrix,
                     transfer.matrix=transfer.matrix,
                     transfer.matrix.sd=transfer.matrix.sd,
                     beta_r.True=TrueValue$beta_r,
                     beta_r=beta_r,
                     beta_r.sd=beta_r.sd)
  return(summary.hmm)
}

summary.hmmgroups.missing=function(data,n,TrueValue=NULL){
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
  beta_r=array(dim=c(dim(Sample$beta_r),r))
  beta_r.sd=array(dim=c(dim(Sample$beta_r.sd),r))
  
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
    beta_r[,,i]=tem$beta_r
    beta_r.sd[,,i]=tem$beta_r.sd
    
  }
  if(is.null(TrueValue))
    res=list(
      beta=apply(beta,1:3,mean),beta.sd=apply(beta.sd,1:3,mean),
      sigma=apply(sigma,1:2,mean),sigma.sd=apply(sigma.sd,1:2,mean),
      init.distribution=apply(init.distribution,1,mean),
      init.distribution.sd=apply(init.distribution.sd,1,mean),
      transfer.matrix=apply(transfer.matrix,1:3,mean),
      transfer.matrix.sd=apply(transfer.matrix.sd,1:3,mean),
      beta_r=apply(beta_r,1:2,mean),
      beta_r.sd=apply(beta_r.sd,1:2,mean)
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
      transfer.matrix.sd=apply(transfer.matrix.sd,1:3,mean),
      
      beta_r.TrueValue=TrueValue$beta_r,
      beta_r=apply(beta_r,1:2,mean),
      beta_r.sd=apply(beta_r.sd,1:2,mean)
    )
  
  return(res)
}


#test parameter to generate missing data
sz=1000  #sample size
X=matrix(c(rep(1,sz),runif(sz,0,5),rnorm(sz,0,5)),sz)
beta=rbind(c(-1,0.5,1),c(0,0.5,0.5),c(1,-0.5,-1))
beta.t=rep(beta,3)
dim(beta.t)=c(3,3,3)
sigma=c(0.3,0.25,0.3)
sigma.t=rep(sigma,3)
dim(sigma.t)=c(3,3)
Prob=c(0.3,0.3,0.4)
Prob.t=rep(c(0.8,0.1,0.05,0.15,0.8,0.15,0.05,0.1,0.8),2)
dim(Prob.t)=c(3,3,2)

temp=HMM.data(X,Prob,Prob.t,beta.t,sigma.t,T)
Y=temp$Y
S=temp$S
TV=temp$TrueValue
beta_r=array(c(2,2.5,3,-0.5,-0.5,-1),dim=c(3,2))
TV$beta_r=beta_r
temp.missing=HMM.data.missing(Y,S,beta_r)
Y.missing=temp.missing$Y.missing
R=temp.missing$R
table(R)

b0=solve(t(X)%*%X)%*%t(X)%*%Y[,1]
B0=1.5*diag(3)
c0=1.28
C0=0.36*var(as.vector(Y))
E0=c(1,1,1)
Et=Prob.t*2
m=4000
system.time(g1<-mcmc.hmm(Y,R,X,b0,B0,c0,C0,E0,Et,S,m,beta_r,lambda_r = 8,lambda_y = 8,rep=1))
summary(g1,m/2,TV)
# beta_r0=array(rep(c(log(length(R)/sum(R==0)-1),0),each=dim(beta_r)[1]),dim=dim(beta_r))
# system.time(g2<-mcmc.hmm(Y.missing,R,X,b0,B0,c0,C0,E0,Et,S,m=10000,beta_r0,lambda_r = 8,lambda_y = 5,rep=10))
# summary(g2,m/2,TV)
