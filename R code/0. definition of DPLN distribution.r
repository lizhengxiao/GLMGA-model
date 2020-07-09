# define the DPLN distribution
library(NormalLaplace)


dDPLN <- function(y, lambda1, lambda2, v, tau, log.p = FALSE) {
  p1 <- exp(0.5*tau^2*lambda1^2 - lambda1*(log(y) - v))*pnorm(q = (log(y) - tau^2*lambda1 - v)/tau, 0, 1)
  p2 <- exp(0.5*tau^2*lambda2^2 - lambda2*(v - log(y)))*(1 - pnorm(q = (log(y) + tau^2*lambda2 - v)/tau, 0, 1))
  if(log.p == TRUE){
    out <- log(lambda1) + log(lambda2) - log(lambda1 + lambda2) - log(y) + log(p1 + p2)
  } else {
    out <- lambda1*lambda2/(lambda1 + lambda2)*(1/y)*(p1 + p2) 
  }
  out
}


pDPLN2 <- function(y, lambda1, lambda2, v, tau) {
  myfun <- function(y){dDPLN(y, lambda1 = lambda1, lambda2 = lambda2, v = v, tau = tau, log.p = FALSE)}
  myfun <- Vectorize(myfun)
  out <- pracma::integral(myfun, xmin = 0, xmax = y)
  return(out)
}
pDPLN2 <- Vectorize(pDPLN2)


pDPLN <- function(y, lambda1, lambda2, v, tau) {
  out <- pnl(log(y), mu = v, sigma = tau, alpha = lambda1, beta = lambda2)
  
  if(is.na(out)|out==-Inf|out==Inf){
    out <- pDPLN2(y, lambda1 = lambda1, lambda2 = lambda2, v = v, tau = tau)
  } 
  out
}
pDPLN <- Vectorize(pDPLN)


# pDPLN1(c(10000, 10), lambda1, lambda2, v, tau)
# pDPLN2(c(10000, 10), lambda1, lambda2, v, tau)
# pDPLN3(c(10000, 10), lambda1, lambda2, v, tau)
#pDPLN4(c(10000, 2000), lambda1, lambda2, v, tau)
# 0.4981106


qDPLN <- function(u, lambda1, lambda2, v, tau) {
  z <- qnl(u, mu = v, sigma = tau, alpha = lambda1, beta = lambda2)
  out <- exp(z)
  return(out)
}




rDPLN <- function(n, lambda1, lambda2, v, tau) {
  T1 <- rexp(n, rate = 1)
  T2 <- rexp(n, rate = 1)
  Z <- rnorm(n, 0, 1)
  # logX0 <- rnorm(n, v, sd = tau)
  Y <- v + tau*Z + T1/lambda1 - T2/lambda2
  #Y <- rnl(n, mu = v, sigma = tau, alpha = lambda1, beta = lambda2)
  xsim <- exp(Y)
  return(xsim)
}

