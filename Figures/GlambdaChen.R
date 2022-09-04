#===============================================
#                 Chen
#===============================================
lambdaChen<-function(u){
  u=diag(length(u))%*%u
  df=dim(u)[2]
  M=rep(0,df)
  k=0
  gama=1
  tol=1e-6
  dif=1
  # R<-function(lam){R0=sum(log(1+t(lam)%*%t(u)));return(R0)}
  R1=rep(0,df)
  R2=R1%*%t(R1)
  
  while(dif>tol && k<=300){
    # 计算R1、R2
    aa=1+t(M)%*%t(u)
    for(i in 1:df){
      R1[i]=sum(t(u[,i])/aa)
      for(j in 1:df){
        R2[i,j]=-sum(u[,i]*u[,j]/aa^2)
      }
    }
    
    delta=-solve(R2)%*%R1
    dif=c(sqrt(t(delta)%*%delta))
    sigma=gama*delta
    while(min(1+t(M+sigma)%*%t(u))<=0){
      gama=gama/2
      sigma=gama*delta
    }
    
    M=M+sigma
    gama=1/sqrt(k+1)
    k=k+1
  }
  # cat(k,'\n')
  return(c(M))
} 
#===================================================#
