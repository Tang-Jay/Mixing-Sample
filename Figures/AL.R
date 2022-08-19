# =============================================== #
#                   计算置信域长度                #
# =============================================== #
AL<-function(x,X,h){
  source('GlambdaChen.R')
  source('K.R')
  source('f.R')
  
  n = length(X)
  if(1>log(n)/2) an=1 else an=log(n)
  u = (x-X)/h
  cut = exp(-qchisq(0.95,1)/2)
  cut_u = qnorm(0.975)
  tol = 1e-08
  L = 0; aL = 0
  U = 0; aU = 0
  times = 0; atimes = 0
  elrMax_x = 0; aelrMax_x = 0
  elRatio = c(0); aelRatio = c(0)
  
  thetas = seq(0,1,0.0002)
  for(theta in thetas){
    # 计算EL值
    z = K2(u)/h-theta
    lam = lambdaChen(z)
    npi = 1/(1+lam*z)
    elr = prod(npi)
    # el = 2*sum(log(1+t(lam)%*%t(z)))
    # (-2)*log(elr)==el
    if(elr>max(elRatio)){elrMax_theta=theta}
    elRatio = c(elRatio,elr)
    if(elr>=cut && times==0){L=theta;times=1}
    if(elr>=cut && times==1){U=theta}
    
    # 计算AEL值
    az=c(z,-an*mean(z))
    alam=lambdaChen(az)
    anpi = 1/(1+alam*az)
    aelr = prod(anpi)
    # ael = 2*sum(log(1+t(alam)%*%t(az)))
    # (-2)*log(aelr)==ael
    if(aelr>max(aelRatio)){aelrMax_theta=theta}
    aelRatio = c(aelRatio,aelr)
    if(aelr>=cut && atimes==0){aL=theta;atimes=1}
    if(aelr>=cut && atimes==1){aU=theta}
    
  }
  
  ELAL = U-L
  AELAL = aU-aL
  
  fnx = mean(K2(u))/h
  z = K2(u)/h-fnx
  sigma2 = h*mean(z^2)
  naL = fnx-sqrt(sigma2/n/h)*cut_u
  naU = fnx+sqrt(sigma2/n/h)*cut_u
  NAAL = naU-naL
  
  return(list(L=L,U=U,aU=aU,aL=aL,naL=naL,naU=naU,ELAL=ELAL,AELAL=AELAL,NAAL=NAAL))
}



