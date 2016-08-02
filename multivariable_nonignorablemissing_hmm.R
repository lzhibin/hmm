library(MASS)
library(dplyr)
library(gtools) #dirichlet distribution
library(actuar) #rinvgamma random variable

HMMdata_missing=function(X,Beta,Sigma,Pi0,Pit,logit_beta,na=NA){
  ### X is N*P ; N is sample size P is the covariate size
  ### Beta is H*P*Ti*M ; M is dimension of Y ; H is hidden variable levels
  # t is the time periods
  ### Sigma is H*Ti*M 
  ### Y is N*Ti*M
  ### S is N*Ti
  ### logit_beta H*(M+1) matrix
  N=nrow(X)
  M=dim(Beta)[4]
  P=dim(Beta)[2]
  H=dim(Beta)[1]
  Ti=dim(Beta)[3]
  Y=array(dim=c(N,Ti,M))
  R=array(0,dim=c(N,Ti,M))
  S=array(dim=c(N,Ti))
  S[,1]=sample.int(H,size = N,replace=T,prob = Pi0)
  for(ni in 1:N){
    for (ti in 2:Ti) {
      S[ni,ti]=sample.int(H,1,prob = Pit[S[ni,ti-1],,ti-1]) 
    }
  }
  for(ni in 1:N){
    for (ti in 1:Ti) {
      Y[ni,ti,]=t(X[ni,]) %*% Beta[S[ni,ti],,ti,] %>%
        as.vector %>%
        mvrnorm(n=1,mu=.,Sigma=diag(Sigma[S[ni,ti],ti,]))
    }
  }
  Y.missing=Y
  for (ni in 1:N) {
    for (ti in 1:Ti) {
      s_yit=S[ni,ti]
      p0=1/(1+exp(sum(logit_beta[s_yit,]*c(1,Y[ni,ti,]))))
      R[ni,ti,]=sample(c(0,1),M,replace = T,prob = c(p0,1-p0))
    }
  }
  Y.missing[R==1]=na
  result=list(Y=Y,Y.missing=Y.missing,Beta=Beta,Sigma=Sigma,Pi0=Pi0,Pit=Pit,S=S,R=R)
  return(result)
}



mcmc.hmm=function(Y,X,b0,B0,c0,C0,E0,Et,S,R,logistic_beta,m,lambda_y=1,lambda_logit=1,rep=1){
  #initial value
  d=ncol(X)    #numbers of beta
  p=length(E0) #levels of hidden variable 
  Ti=dim(Y)[2]   #numbers of periods
  NY=dim(Y)[1]   #sample sizes
  MY=dim(Y)[3] #dimension of Y
  Bk=array(dim=c(d,d,p,Ti,MY))
  bk=array(dim=c(p,d,Ti,MY))
  
  ck=array(dim=c(p,Ti,MY))
  Ck=array(dim=c(p,Ti,MY))
  betak=array(dim=c(p,d,Ti,MY))
  sigmak=array(dim=c(p,Ti,MY))
  
  Dt=array(dim=c(p,p,Ti-1)) #transfer matrix [f(u_j|u_hat_i)]
  NS=array(dim=c(p,Ti))
  NSt=array(dim=c(p,p,Ti-1))
  if(sum(is.na(Y)>0))
    Y[R==1]=runif(sum(R),min(Y,na.rm = T),max(Y,na.rm = T))
  #random init Y
  
  p.betak=array(dim=c(p,d,Ti,m,MY))
  p.sigmak=array(dim=c(p,Ti,m,MY))
  p.D=array(dim=c(p,m))
  p.Dt=array(dim=c(p,p,Ti-1,m))
  p.logistic_beta=array(dim=c(dim(logistic_beta),m))
  
  B0.inv=list()
  for (myi in 1:MY) {
    B0.inv[[myi]]=solve(B0[,,myi])
  }
  for (myi in 1:MY) {
    temp.beta=mvrnorm(p,b0[,myi],B0[,,myi])
    temp.sigma=rinvgamma(p,c0[myi],scale=C0[myi])
    for(t in 1:Ti){
      betak[,,t,myi]=temp.beta
      sigmak[,t,myi]=temp.sigma
    }
    
  }
  
  y_ratio=function(j,y_old,y_new,x,beta,sigma2,logistic_beta,r_M){
    #this function is to calc the loglikelihood ratio for y_new/y_old
    #y_old and y_new are diff in the jth element
    #beta and sigma are the parameter of the jth element
    #logistic_beta and r is the corresponding vector
    #M is globle variable indicate the dimension of Y
    log_ratio=((y_old[j+1]-sum(x*beta))^2-(y_new[j+1]-sum(x*beta))^2)/(2*sigma2)+
      sum(r_M)*logistic_beta[j+1]*(y_new[j+1]-y_old[j+1])+
      MY*(log(1+exp(sum(logistic_beta*y_old)))-log(1+exp(sum(logistic_beta*y_new))))
    return(exp(log_ratio))
  }
  
  
  logistic_beta_ratio=function(Y,lg_beta_old,lg_beta_new,r_rowsum){
    #logistic_beta to be renew correspond to the kth hidden class
    #Y is the Nk*M matrix
    #Nk is the total number of samples belong to kth hidden class
    #r_rowsum is a vector (length of Nk) indicate total missing elements of every sample
    # belong to the kth hidden class
    delta=lg_beta_new-lg_beta_old
    log_ratio=sum(r_rowsum*(Y%*%delta))+sum(MY*log(1+exp(Y%*%lg_beta_old))-
                                            MY*log(1+exp(Y%*%lg_beta_new)))
    return(exp(log_ratio))
  }
  
  im_y=function(j,sigma,logistic_beta){
    #This function is to calc information matrix when sample Yijt
    #actually return a numerical object
    im_inv=MY*(logistic_beta[j+1]^2)*exp(logistic_beta[1])/((1+exp(logistic_beta[1]))^2)+1/sigma
    im=1/im_inv
    return(im)
  }
  
  im_logistic_beta=function(Y){
    #This function is to calc information matrix when sample logistic_beta
    #return M+1 * M+1 dimension matrix
    im_inv_vec=apply(Y,1,function(x)x%*%t(x))
    im_inv=rowSums(im_inv_vec)/4
    dim(im_inv)=rep(MY+1,2)
    im=solve(im_inv)
    return(im)
    
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
    
    for(myi in 1:MY){
      for(t in 1:Ti){
        for(j in 1:p){
          Bk[,,j,t,myi]=solve(B0.inv[[myi]]+t(X[S[,t]==j,])%*%X[S[,t]==j,]/sigmak[j,t,myi])
          bk[j,,t,myi]=Bk[,,j,t,myi]%*%(B0.inv[[myi]]%*%b0[,myi]+
                                          t(X[S[,t]==j,])%*%Y[S[,t]==j,t,myi]/sigmak[j,t,myi])
          
          betak[j,,t,myi]=mvrnorm(1,bk[j,,t,myi],Bk[,,j,t,myi])
          
          ck[j,t,myi]=c0[myi]+NS[j,t]/2
          Ck[j,t,myi]=C0[myi]+sum((Y[S[,t]==j,t,myi]-X[S[,t]==j,]%*%betak[j,,t,myi])^2)/2
          
          sigmak[j,t,myi]=rinvgamma(1,ck[j,t,myi],scale=Ck[j,t,myi])
          
        }
      }
    }
    #renew Yijt
    for(ni in 1:NY){
      for(t in 1:Ti){
        for(mi in 1:MY){
          if(R[ni,t,mi]==1){
            #temp_yijt_im=im_y(mi,sigmak[S[ni,t],t,mi],logistic_beta[S[ni,t],])
            temp_yijt_im=0.5
            temp_yijt=rnorm(1,Y[ni,t,mi],lambda_y*sqrt(temp_yijt_im))
            temp_old=c(1,Y[ni,t,])
            temp_new=temp_old
            temp_new[mi+1]=temp_yijt
            alpha=y_ratio(mi,temp_old,temp_new,X[ni,t],betak[S[ni,t],,t,mi],
                          sigmak[S[ni,t],t,mi],logistic_beta[S[ni,t],],R[ni,t,])
            U_y=runif(1)
            if(U_y<alpha)
              Y[ni,t,mi]=temp_yijt
          }
        }
      }
    }
    
    # #renew logistic_beta
    # S_vec=as.vector(S)
    # R_col=R
    # dim(R_col)=c(NY*Ti,MY)
    # Y_col=Y
    # dim(Y_col)=c(NY*Ti,MY)
    # Y_col=cbind(rep(1,NY*Ti),Y_col)
    # for (k in 1:p) {
    #   Yk=Y_col[S_vec==k,]
    #   Rk=R_col[S_vec==k,]
    #   temp=mvrnorm(1,logistic_beta[k,],lambda_logit*im_logistic_beta(Yk))
    #   alpha=logistic_beta_ratio(Yk,logistic_beta[k,],temp,rowSums(Rk))
    #   U_logit=runif(1)
    #   if(U_logit<alpha)
    #     logistic_beta[k,]=temp
    # }
    
    # #FFBS renew hidden variable altnative code
    # f.density=function(y,x,beta,sigma){
    #   res=1
    #   for(myi in 1:MY){
    #     m=as.vector(beta[,,myi]%*%x)
    #     s=sqrt(sigma[,myi])
    #     p=dnorm(y[myi],m,s)
    #     p=p/sum(p)
    #     res=res*p
    #   }
    #   return(res)
    # }
    # P.S.Update=array(dim=c(Ti,p,NY))
    # for(nr in 1:NY){
    #   P.S.Update[1,,nr]=D*f.density(Y[nr,1,],X[nr,],betak[,,1,],sigmak[,1,])
    #   for(t in 2:Ti){
    #     P.S.Update[t,,nr]=t(P.S.Update[t-1,,nr])%*%Dt[,,t-1]*
    #       f.density(Y[nr,t,],X[nr,],betak[,,t,],sigmak[,t,])
    #   }
    # }
    # for(nr in 1:NY){
    #   S[nr,Ti]=sample.int(p,1,prob = P.S.Update[Ti,,nr])
    #   for(t in (Ti-1):1){
    #     s.temp=P.S.Update[t,,nr]*Dt[,S[nr,t+1],t]
    #     S[nr,t]=sample.int(p,1,prob = s.temp)
    #   }
    # }
    
    
    #label switch
    base=betak[,1,,1]
    for(t in 1:Ti){
      for(myi in 1:MY){
        betak[,,t,myi]=betak[order(base[,t]),,t,myi]
        sigmak[,t,myi]=sigmak[order(base[,t]),t,myi]
      }
      S[,t]=rank(base[,t])[S[,t]]
    }
    D=D[order(base[,1])]
    logistic_beta=logistic_beta[order(base[,1]),]
    for(t in 1:(Ti-1)){
      #switch Dt
      Dt[,,t]=Dt[order(base[,t]),order(base[,t+1]),t]
      
    }
    #print(paste0("error rate:",sum(S-S.T != 0))) #for test
    
    #store parameters
    p.betak[,,,i,]=betak
    p.sigmak[,,i,]=sigmak
    p.D[,i]=D
    p.Dt[,,,i]=Dt
    p.logistic_beta[,,i]=logistic_beta
  }
  ParaMeters=list(p.betak=p.betak,p.sigmak=p.sigmak,p.D=p.D,p.Dt=p.Dt,p.logistic_beta=p.logistic_beta)
  class(ParaMeters)="hmm"
  if(rep==1)
    return(ParaMeters)
  else{
    groups=list()
    groups[[1]]=ParaMeters
    for(r in 2:rep)
      groups[[r]]=gibbs.hmm(Y,X,b0,B0,c0,C0,E0,Et,S,R,logistic_beta,m,lambda_y,lambda_logit,rep=1)
    class(groups)='hmmgroups'
  }
  return(groups)
}

#method to summary hmm class
summary.hmm=function(data,n){
  len=dim(data$p.D)[2]
  if(len<n){
    print("error! need more data")
  }
  
  beta=apply(data$p.betak[,,,(len-n+1):len,],c(1,2,3,5),mean)
  beta.sd=apply(data$p.betak[,,,(len-n+1):len,],c(1,2,3,5),sd)
  sigma=apply(data$p.sigmak[,,(len-n+1):len,],c(1,2,4),mean)
  sigma.sd=apply(data$p.sigmak[,,(len-n+1):len,],c(1,2,4),sd)
  init.distribution=apply(data$p.D[,(len-n+1):len],1,mean)
  init.distribution.sd=apply(data$p.D[,(len-n+1):len],1,sd)
  transfer.matrix=apply(data$p.Dt[,,,(len-n+1):len],c(1,2,3),mean)
  transfer.matrix.sd=apply(data$p.Dt[,,,(len-n+1):len],c(1,2,3),sd)
  logit_beta=apply(data$p.logistic_beta[,,(len-n+1):len],c(1,2),mean)
  logit_beta.sd=apply(data$p.logistic_beta[,,(len-n+1):len],c(1,2),sd)
  
  summary.hmm=list(beta=beta,beta.sd=beta.sd,sigma=sigma,sigma.sd=sigma.sd,
                   init.distribution=init.distribution,
                   init.distribution.sd=init.distribution.sd,
                   transfer.matrix=transfer.matrix,
                   transfer.matrix.sd=transfer.matrix.sd,
                   logit_beta=logit_beta,
                   logit_beta.sd=logit_beta.sd)
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

sz=1000
MY=2
X=cbind(rep(1,sz),runif(sz,0,5),rnorm(sz,2))
Sigma=array(0.3,dim=c(3,3,2))
Pi0=c(0.3,0.4,0.3)
Pit=array(rep(c(0.6,0.2,0.15,0.25,0.6,0.25,0.15,0.2,0.6),2),dim=c(3,3,2))
beta_y1s1=c(-1,0.5,-1)
beta_y1s2=c(0,1,-1)
beta_y1s3=c(1,-0.5,1)
beta_y2s1=c(0,1,-1)
beta_y2s2=c(1,-0.5,1)
beta_y2s3=c(-1,0.5,-1)
Beta_t1=array(dim = c(3,3,2))
Beta_t1[1,,1]=beta_y1s1
Beta_t1[2,,1]=beta_y1s2
Beta_t1[3,,1]=beta_y1s3
Beta_t1[1,,2]=beta_y2s1
Beta_t1[2,,2]=beta_y2s2
Beta_t1[3,,2]=beta_y2s3
Beta=array(dim=c(3,3,3,2))
for (ti in 1:3) {
  Beta[,,ti,]=Beta_t1
}
logit_beta=matrix(c(-3,1,1,-2,0.5,0.5,-1,0.25,0.75),3,3,byrow=T)
HD=HMMdata_missing(X,Beta,Sigma,Pi0,Pit,logit_beta)
Y=HD$Y
Y.missing=HD$Y.missing
R=HD$R

# table(test$S[,1])/sz
# pit1=table(test$S[,1],test$S[,2])
# pit2=table(test$S[,2],test$S[,3])
# pit1/colSums(pit1)
# pit2/colSums(pit2)
S.T=HD$S
S=sample.int(length(Pi0),dim(Y)[1]*dim(Y)[2],replace = T)
dim(S)=dim(Y)[1:2]
#generate completely random S
b0=array(dim=c(ncol(X),MY))
for(myi in 1:MY){
  b0[,myi]=solve(t(X)%*%X)%*%t(X)%*%Y[,1,myi]
}
B0=array(rep(1.5*diag(3),MY),dim=c(3,3,2))
c0=rep(1.28,MY)
C0=0.36*apply(Y,3,sd)^2
E0=c(1,1,1)
Et=Pit*2
logit_beta_random=matrix(c(-1,-1,-1,0,0,0,0,0,0),3,3)
m=2000
Sys.time()
system.time(g2<-mcmc.hmm(Y,X,b0,B0,c0,C0,E0,Et,S.T,R,logit_beta,m,lambda_y=1,lambda_logit=1,rep=1))
Sys.time()
summary(g,m/2)

