library(MASS)
library(glmnet)
#library(SIS)
#library(screening)

library(scalreg)
set.seed(1)

ρ1=0.5
ρ2=0.3
#ρ3=0.2
p=300
n=100
M=100
#length=27
length=1.4*round(n/log(n)) 
lambda.univ=(2.025/n*log(p))^0.5
alpha=0.05
theta=diag(1,p,p)   

#setting1
for (i in 1:p){
  for (j in 1:p){
    if(abs(i-j)==1){
      theta[i,j]=ρ1
    }
    if(abs(i-j)==2){
      theta[i,j]=ρ2
    }
    # if(abs(i-j)==3){
    #  theta[i,j]=ρ3
    #}
  }
}

sigma=solve(theta)              
mean=c(rep(0,p))               



gamma0=matrix(0,(p-1),p)           
for (j in 1:(p)){
  gamma0[,j]=-theta[-j,j]/theta[j,j]
}
timestart=Sys.time()


bingxing<-function(n,p,mean,sigma,gamma0,length){
  colsum=c(rep(0,p))
  
  I=diag(1,n,n)
  p_valu=matrix(0,p-1,p)
  #p_valuA= matrix(0,p,p)
  gamma_hat= matrix(0,p-1,p)
  
  X = matrix(0,n,p)   
  
  
  X = mvrnorm(n=n,mean,sigma)     
 
  for (j in 1:p){
    Xjbu = X[,-j]  
    Xj   = X[,j]
    gujilingling=scalreg(Xjbu,Xj,lam0=lambda.univ) 
    # gujilingling=scalreg(Xjbu,Xj)
    gammajini=coef(gujilingling) 
    tausj = gujilingling$hsigma 
    taus = tausj
    
    #SIS 
    #omegaj=t(Xjbu)%*%Xj
    #Holp
    omegaj=t(Xjbu)%*%(solve(Xjbu%*%t(Xjbu)))%*%Xj 
    BICj1 = c(rep(0,length))    
    min = 2
    for(sj in 2:length){
      Sj = sort(order(abs(omegaj),decreasing=TRUE)[1:sj])
      xj = Xjbu[,Sj]%*%solve(t(Xjbu[,Sj])%*%Xjbu[,Sj])%*%t(Xjbu[,Sj])%*%Xj 
      residualj = Xj-xj
      resj = sum(residualj^2)                                    
      
      BICj1[sj] = resj/(n)+ sj*(taus^2)*log(p)*log(log(n))/n
      if(BICj1[sj] <= BICj1[min]){  
        min = sj                                                 
      }
    }
    #Min[k] = min
    Sj.hat = sort(order(abs(omegaj),decreasing=TRUE)[1:min])           
    Xj2 = Xjbu[,Sj.hat]                                         
    Pj = Xj2%*%solve(t(Xj2)%*%Xj2)%*%t(Xj2) 
    
    Sjbu=setdiff(1:(p-1),Sj.hat) 
    lamdazheng=(2/n*log(length(Sjbu)))^0.5 
    
    
    psij= matrix(0,n,p-1) 
    zj = matrix(0,n,p-1) 
    zetaj= c(rep(0,p-1))
    gamma_hatj = c(rep(0,p-1))
    p_valuj = c(rep(0,p-1))
    
    for (i in 1:(p-1)){
      if (i %in% Sj.hat){  
        Xj3 = Xjbu[,setdiff(Sj.hat,i)]
        #psij= (I-Xj3%*%solve(t(Xj3)%*%Xj3)%*%t(Xj3))%*%Xjbu
        psij[,Sjbu]=(I-Xj3%*%solve(t(Xj3)%*%Xj3)%*%t(Xj3))%*%Xjbu[,Sjbu]
        psij[,i]= (I-Xj3%*%solve(t(Xj3)%*%Xj3)%*%t(Xj3))%*%Xjbu[,i]
        fitj=glmnet(psij[,Sjbu], psij[,i],intercept=0) 
        keaij=coef(fitj)[-1,] 
        meihaoj=ncol(keaij) 
        BIC2j = c(rep(0,meihaoj))
        min1 = 1
        for(sj1 in 1:meihaoj){
          residualj=psij[,i]-psij[,Sjbu]%*%keaij[,sj1]
          resj = sum(residualj^2)
          BIC2j[sj1] = log(resj) + length(which(keaij[,sj1]!=0))*log(p)*log(log(n))/n            
          if(BIC2j[sj1] <= BIC2j[min1]){
            min1 = sj1
          }
        }
        
        zj[,i] =psij[,i]- psij[,Sjbu]%*%keaij[,min1]
      }
      else if (i %in% Sjbu){
        
        psij[,Sjbu]=(I-Pj)%*%Xjbu[,Sjbu]
        fitj=glmnet(psij[,setdiff(Sjbu,i)], psij[,i],intercept=0)
        keaij=coef(fitj)[-1,]
        meihaoj=ncol(keaij)
        BIC2j = c(rep(0,meihaoj))
        min1 = 1
        for(sj1 in 1:meihaoj){
          residualj= psij[,i]- psij[,setdiff(Sjbu,i)]%*%keaij[,sj1] 
          resj = sum(residualj^2)
          BIC2j[sj1] = log(resj) + length(which(keaij[,sj1]!=0))*log(p)*log(log(n))/n
          if(BIC2j[sj1] <= BIC2j[min1]){
            min1 = sj1
          }
        }
        
        zj[,i] =psij[,i]- psij[,setdiff(Sjbu,i)]%*%keaij[,min1]
      }
      
      zetaj[i]=((t(zj[,i])%*%zj[,i])^0.5)/abs(t(Xjbu[,i])%*%zj[,i])
    
      gamma_hatj[i] = (t(zj[,i])%*%Xj)/(t(Xjbu[,i])%*%zj[,i])
      
      p_valuj[i] =2*(1-pnorm(abs(gamma_hatj[i])/(zetaj[i]*taus)) )
      
    }
    gamma_hat[,j]= t(gamma_hatj) 
    p_valu[,j]= t(p_valuj)
    
  }
  #return(p_valu)
  obj=list(p_valu)
}
#install.packages('doParallel')  #当运行提示没有包的时候，需要安装一下
#install.packages('parallel')
#install.packages('foreach')
#install.packages('iterators')
#install.packages('doParallel')

require(doParallel) #同样是载入包的操作
require(foreach)

cl = makeCluster(20) #初始化核心集群，加快进程（提升效率）
registerDoParallel(cl) # 这是设置进程的过程

resultmcp = foreach (t = 1:M, .combine = 'rbind',
                     .packages = c('MASS','scalreg','glmnet')) %dopar% {
                       bingxing(n,p,mean,sigma,gamma0,length)
                     } 



stopImplicitCluster()
stopCluster(cl) 

a=array(0,c(p-1,p,M))

for (t in 1:M){
  for(jj in 1:p){
    a[,jj,t]=resultmcp [[t]][,jj]
  }
}
s1=matrix(0,p-1,p) 
for(jj in 1:p){
  for (ii in 1:(p-1)){
    s1[ii,jj]=sum(I(a[ii,jj,]< 0.05))     
  }
}
s=matrix(0,p,p) 
for(jj in 1:p){
  s[-jj,jj]= s1[,jj]    
}
EP=matrix(0,p,p)
for(ii in 1:p){
  for (jj in 1:p){
    if(abs(ii-jj)==1){
      EP[ii,jj]=s[ii,jj]
    }
    if(abs(ii-jj)==2){
      EP[ii,jj]=s[ii,jj]
    }
    # if(abs(ii-jj)==3){
    #  EP[ii,jj]=s[ii,jj]
    #}
  }
}


ET=s-EP


ETs=sum(ET)/((p*p-p-4*p+6)*M) 
ETs
EPs=sum(EP)/((p*4-6)*M)

EPs

timeend=Sys.time()
runningtime=timeend-timestart
print(runningtime)
