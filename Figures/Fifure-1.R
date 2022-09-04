# =============================================== #
#                   Figure-1                      #
# =============================================== #
rm(list = ls())
source('GlambdaChen.R')
source('K.R')
source('f.R')
source('AL.R')
figure_AL_as_h<-function(n){
  X = c(0)
  for(j in 2:(n+100)){
    X[j] = 0.5*X[j-1] + rnorm(1,0,1)
  }
  X = X[101:(n+100)]
  AL1 = c()
  AL2 = c()
  AL3 = c()
  for(h in H){
    re = AL(x=x,X=X,h=h)
    AL1 = c(AL1,re$ELAL)
    AL2 = c(AL2,re$AELAL)
    AL3 = c(AL3,re$NAAL)
  }
  lines(H,AL2,type='l',col='brown2',lty=1,lwd=1.5)
  lines(H,AL1,type='l',col='aquamarine4',lty=6,lwd=1.5)
  lines(H,AL3,type='l',col='blue1',lty=1,lwd=1.5)
  # text(0.35,AL2[3]+0.02,label = paste0('n=',n))
  cat('完成画图 n =',n,'\n')
}

# =============================================== #
#                  Change values                  #                  #
# =============================================== #
nsim = 1000
n = 50
x = 0
seedNum = 15 
# =============================================== #
#                  Figure-1-Left                  #
# =============================================== #
set.seed(seedNum)
theta = f(x,0,4/3)
cut = qchisq(0.95,1)
cut_u = qnorm(0.975)
an=log(n)
CP1 = c()
CP2 = c()
CP3 = c()
H = seq(0.1,0.6,0.05)
for(h in H){
  cp1 = 0
  cp2 = 0
  cp3 = 0
  for(i in 1:nsim){
    X = c(0)
    for(j in 2:(n+100)){
      X[j] = 0.5*X[j-1] + rnorm(1,0,1)
    }
    X = X[101:(n+100)]
    u = (x-X)/h
    
    # EL CP
    z = K2(u)/h-theta
    lam = lambdaChen(z)
    el = 2*sum(log(1+t(lam)%*%t(z)))
    if(el<cut) cp1 = cp1+1
    
    # AEL CP
    az=c(z,-an*mean(z))
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    if(ael<cut) cp2=cp2+1
    
    # NA CP
    fnx = mean(K2(u))/h
    z = K2(u)/h-fnx
    sigma2 = h*mean(z^2)
    L = fnx-sqrt(sigma2/n/h)*cut_u
    U = fnx+sqrt(sigma2/n/h)*cut_u
    if(L<theta & theta<U) cp3 = cp3+1
  }
  cat(h,cp2/nsim,cp1/nsim,cp3/nsim,'\n')
  CP1 = c(CP1,cp1/nsim)
  CP2 = c(CP2,cp2/nsim)
  CP3 = c(CP3,cp3/nsim)
}
par(mfrow=c(1,2))
plot(H,H,col='white',ylim=c(0.84,0.96),xlab='',ylab='',yaxt='n')
axis(2,las = 1)
grid(nx=NULL,ny=NULL,col='lightgray',
     lty='dotted',lwd=par('lwd'),equilogs=T)
legend('bottomright',legend=c('AELCP',' ELCP',' NACP'),
       col=c('brown2','aquamarine4','blue1'),
       pch=c(1,1,1),cex=0.88,
       lty=c(1,6,6),bty='o',lwd=1.5)
lines(H,CP2,type='o',col='brown2',lty=1,lwd=1.5,cex=0.7)
lines(H,CP1,type='o',col='aquamarine4',lty=6,lwd=1.5,cex=0.7)
lines(H,CP3,type='o',col='blue1',lty=6,lwd=1.5,cex=0.7)
cat('在x =',x,',n=',n,'时cp随h变化趋势图完成','\n')

# =============================================== #
#                  Figure-1-Right                 #
# =============================================== #
plot(H,H,col='white',ylim=c(0.30,1.2),xlab='',ylab='',yaxt='n')
axis(2, at=seq(0.30,1.2,0.15), las = 1)
grid(nx=NULL,ny=NULL,col='lightgray',
     lty='dotted',lwd=par('lwd'),equilogs=T)
legend('bottomleft',legend=c('AELAL','   ELAL','   NAAL'),
       col=c('brown2','aquamarine4','blue1'),
       lty=c(1,6,1),bty='o',lwd=1.5,cex=0.82)
set.seed(seedNum)
figure_AL_as_h(n = n)
cat('在x =',x,'时AL随h变化趋势图完成','\n')
# dev.off()
