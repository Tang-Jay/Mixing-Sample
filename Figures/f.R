f <- function(x, mu = 0, sigma2 = 1){
  fx = (2*pi*sigma2)^(-0.5)*exp(-(x-mu)^2/(2*sigma2))
  return(fx)
}
