# =============================================== #
#                  Table-1  CP                    #
# =============================================== #
rm(list = ls())
source('GlambdaChen.R')
source('K.R')
source('f.R')
# =============================================== #
#                  Change values                  #
# =============================================== #
nsim = 2000 
x = 0  # x = 0,1
size = c(50,100,150,200,250)
theta = f(x,0,4/3)
cut = qchisq(0.95,1)
cut_u = qnorm(0.975)
# =============================================== #
#             Simulation begin here               #
# =============================================== #
tt1 = 0
tt2 = 0
tt3 = 0
cp_el = c()
cp_ael = c()
cp_na = c()
cp_time = c()
cp_atime = c()
cp_ntime = c()
T1 = lubridate::now()
cat('x',' n','  AEL ','  EL ','  NA','\n')
for(n in size){
  h = 4*n^(-0.2)/log(n)
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
    t1 = lubridate::now()
    lam = lambdaChen(as.matrix(z))
    el = 2*sum(log(1+t(lam)%*%t(z)))
    t2 = lubridate::now()
    if(el<cut) cp1 = cp1+1
    tt1=tt1+t2-t1
    
    # AEL CP
    if(1>log(n)/2) an=1 else an=log(n) # or an=log(n)/2
    az=c(z,-an*mean(z))
    at1 = lubridate::now()
    alam = lambdaChen(as.matrix(az))
    ael = 2*sum( log(1+t(alam)%*%t(az) ) )  
    at2 = lubridate::now()
    if(ael<cut) cp2 = cp2+1
    tt2 = tt2+at2-at1
    
    # NA CP
    nt1 = lubridate::now()
    fnx = mean(K2(u))/h
    z = K2(u)/h-fnx
    sigma2 = h*mean(z^2)
    L = fnx-sqrt(sigma2/n/h)*cut_u
    U = fnx+sqrt(sigma2/n/h)*cut_u
    nt2 = lubridate::now()
    if(L<theta & theta<U) cp3 = cp3+1
    tt3 = tt3+nt2-nt1
  }
  cat(x,n,cp2/nsim,cp1/nsim,cp3/nsim,'\n')
  
  cp_el=append(cp_el,cp1/nsim)
  cp_ael=append(cp_ael,cp2/nsim)
  cp_na=append(cp_na,cp3/nsim)
  
  cp_time=append(cp_time,tt1)
  cp_atime=append(cp_atime,tt2)
  cp_ntime=append(cp_ntime,tt3)
  
  tt1=0
  tt2=0
  tt3=0
}
T2 = lubridate::now()
TT = T2-T1
# =============================================== #
#             Display of results                  #      
# =============================================== #
# CP = matrix(NA,nrow=length(size),ncol=3)
# CP[,1]=cp_ael
# CP[,2]=cp_el
# CP[,3]=cp_na
# rownames(CP)  =  size
# colnames(CP)  =  c("AEL_cps","EL_cps","NA_cps")
# print(CP)

Time = matrix(NA,nrow=length(size),ncol=3)
Time[,1]=cp_atime
Time[,2]=cp_time
Time[,3]=cp_ntime
rownames(Time)  =  size
colnames(Time)  =  c("AEL_times","EL_times","NA_times")
print(Time)

cat('\n')
cat('The total running time is',TT,'sec.\n')
cat('AEL beats EL time by',sum(cp_time)-sum(cp_atime),'seconds','\n')
cat(sum(cp_atime),"+",sum(cp_time),"+",sum(cp_ntime),'\n')



