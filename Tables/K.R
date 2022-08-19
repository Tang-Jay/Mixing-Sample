K1 <- function(u, a = 1, sig2 = 1){
  k1 = 3/(4*a)*(1-(u/a)^2)^2*(abs(u)<=a)
  return(k1)
}

K2 <- function(u, a = 1, sig2 = 1){
  k2 = 15/(16*a)*(1-(u/a)^2)^2*(abs(u)<=a)
  return(k2)
}

K3 <- function(u, a = 1, sig2 = 1){
  k3 = 1/a*(1-abs(u/a))*(abs(u)<=a)
  return(k3)
}

K4 <- function(u, a = 1, sig2 = 1){
  k4 = 0.5/a*(abs(u)<=a)
  return(k4)
}

K5 <- function(u, a = 1, sig2 = 1){
  k5 = (2*pi*sig2)^(-0.5)*exp(-u^2/(2*sig2))
  return(k5)
}