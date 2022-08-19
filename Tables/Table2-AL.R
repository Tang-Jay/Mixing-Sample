# =============================================== #
#                    Table-2 AL                   #
# =============================================== #
rm(list = ls())
source('GlambdaChen.R')
source('K.R')
source('f.R')
# =============================================== #
#                   Change values                 #
# =============================================== #
x = 0   # x = 0,1
m = 100 # Remove the first m samples
set.seed(10)
# =============================================== #
#               Simulation begin here             #                  #
# =============================================== #
cat('x',' n',' ELAL',' AELAL',' NA','\n')
for(n in c(50,100,150,200,250)){
  h = 4*n^(-0.2)/log(n)
  if(1>log(n)/2) an=1 else an=log(n)
  X = c(0,0)
  for(i in 3:(n+m)){
    X[i] = 0.75*X[i-1] - 0.125*X[i-2] + rnorm(1,0,1)
  }
  X = X[(m+1):(n+m)]
  u = (x-X)/h
  cut = exp(-qchisq(0.95,1)/2)
  cut_u = qnorm(0.975)
  
  tol = 1e-08
  L = 0; aL = 0
  U = 0; aU = 0
  times = 0; atimes = 0
  elrMax_theta = 0; aelrMax_theta = 0
  elRatio=c(0);aelRatio=c(0)
  thetas = seq(0,1,0.0001)
  for(theta in thetas){
    # EL
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
    
    # AEL
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
  
  # figure
  plot(thetas,thetas,ylim=c(0,1),col='white',xlab='theta',ylab='elr',
       main=paste0('x = ',x,', n = ',n),yaxt='n')
  axis(2,las = 1)
  legend('topright',legend=c('AELCP','  ELCP'),
         col=c(2,3),lty=c(1,1),bty='n',lwd=1.5)
  # plot EL
  points(thetas,elRatio[-1],type='o',lty=1,lwd=0.2,pch=20,col=3)
  abline(h=cut)
  abline(v=elrMax_theta,col=3)
  abline(v=L,col=3)
  abline(v=U,col=3)
  text(elrMax_theta,0.6,label = elrMax_theta,col=3)
  text(U+0.05,0.6,label = U,col=3)
  text(L-0.05,0.6,label = L,col=3)
  # plot AEL
  points(thetas,aelRatio[-1],type='o',lty=1,lwd=0.2,pch=20,col=2)
  abline(v=aelrMax_theta,col=2)
  abline(v=aL,col=2)
  abline(v=aU,col=2)
  text(aelrMax_theta,0.5,label = aelrMax_theta,col=2)
  text(aU+0.05,0.5,label = aU,col=2)
  text(aL-0.05,0.5,label = aL,col=2)
  # Drawing completed 
  
  # And then to calculate ALs of all methods
  ELAL = U-L
  AELAL = aU-aL
  fnx = mean(K2(u))/h
  z = K2(u)/h-fnx
  sigma2 = h*mean(z^2)
  NAAL = sqrt(sigma2/n/h)*cut_u*2
  
  cat(x,n,ELAL,AELAL,round(NAAL,4),'\n')
}





