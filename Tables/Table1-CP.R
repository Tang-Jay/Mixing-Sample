# =============================================== #
#                  Table-1  CP                    #                  #
# =============================================== #
rm(list = ls())
source('GlambdaChen.R')
source('K.R')
source('f.R')
# =============================================== #
#                  Change values                  #                  #
# =============================================== #
nsim = 2000 
x = 0  # x = 0,1
theta = f(x,0,4/3)
cut = qchisq(0.95,1)
cut_u = qnorm(0.975)
# =============================================== #
#             Simulation begin here               #                  #
# =============================================== #
cat('x','n','  AEL ',' EL',' NA','\n')
for(n in c(50,100,150,200,250)){
  h = 4*n^(-0.2)/log(n)
  cp1 = 0
  cp2 = 0
  na = 0
  for(i in 1:nsim){
    X = c(0)
    for(j in 2:(n+100)){
      X[j] = 0.5*X[j-1] + rnorm(1,0,1)
    }
    X = X[101:(n+100)]
    u = (x-X)/h
    
    # EL CP
    z = K2(u)/h-theta
    lam = lambdaChen(as.matrix(z))
    el = 2*sum(log(1+t(lam)%*%t(z)))
    if(el<cut) cp1 = cp1+1
    
    # AEL CP
    if(1>log(n)/2) an=1 else an=log(n)#/2
    az=c(z,-an*mean(z))
    alam=lambdaChen(as.matrix(az))
    ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    if(ael<cut) cp2=cp2+1
    
    # NA CP
    fnx = mean(K2(u))/h
    z = K2(u)/h-fnx
    sigma2 = h*mean(z^2)
    L = fnx-sqrt(sigma2/n/h)*cut_u
    U = fnx+sqrt(sigma2/n/h)*cut_u
    if(L<theta & theta<U) na = na+1
  }
  cat(x,n,cp2/nsim,cp1/nsim,na/nsim,'\n')
}
  




