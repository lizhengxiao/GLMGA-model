

# 1. Descriptive analysis -------------------------------------------------
library(gamlss)
library(data.table)
library(actuar)
library(extRemes)
library(POT)
source('0. definition of logmoyal distribution.r')
source('0. definition of logmoyal-gamma distribution.r')
setwd("D:/Github - Lee/GLMGA-model")
dtnew <- fread('/eqdata.csv', header = T)
y <- dtnew$ynew

# difine the hessian matrix output and the stardard error -------------------------------------------------
modout <- function(m) {
  Hessian <- m$hessian
  se <- sqrt(diag(solve(Hessian)))                  ## standard error
  Z <- m$par/se                                             ##  Z statistics
  p <- ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)           ## p value
  summarytable <- round(data.frame(m$par, se, Z, p), 3)
  NLL <- m$value  # - loglikelihood value
  AIC <- 2*length(m$par) + 2*NLL
  BIC <- log(length(y))*length(m$par) + 2*NLL
  list(summary = summarytable, ll =  - m$value, AIC = AIC, BIC = BIC)
}

# base model: GP distribution -------------------------------------------------
umin <- 0
GPD2.MLE <- fevd(x = ynew, 
                 scale.fun = ~ 1, 
                 shape.fun = ~ magnitude + intensity, 
                 threshold = 0, method = 'MLE', use.phi = F,
                 type = "GP", data = dtnew)
GPD2.MLE
GPD2.par <- GPD2.MLE$results$par
pars <- GPD2.MLE$results$par
scale <- pars[1]
shape <- X%*%pars[2:4]
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pgpd(y[i], loc = umin, scale = scale, shape = shape[i])
}
resGPD <- qnorm(ufit)
qqnorm(resGPD, ylim = c(-3, 3), xlim = c(-3, 3))
abline(0, 1, col = 'blue')

GPD2.hessian <- GPD2.MLE$results$hessian
GPD2.se <- sqrt(diag(solve(GPD2.hessian)))   
GPD2.Z = GPD2.par/GPD2.se                                            ##  Z 
GPD2.p = ifelse(GPD2.Z>=0, pnorm(GPD2.Z, lower=F)*2, pnorm(GPD2.Z)*2)    ## p value
data.table(GPD2.par, GPD2.se, GPD2.p)

mlGPD <- GPD2.MLE$results
modout(mlGPD)

# 1.1 log normal distribution -------------------------------------------------
X <- model.matrix(~ magnitude + intensity, data = dtnew) # desgin matrix
LLlognormal <-  function(pars, y, X){
  lengthpar <- dim(X)[2]
  mu <- (X%*%pars[1:lengthpar])
  sigma <- exp(pars[lengthpar+1])  # log link function
  logL <- dLOGNO(y, mu, sigma, log = T)
  ll <- -sum(logL)
  return(ll)
}
mllognormal <- optim(fn = LLlognormal, 
                     X = X, y = y, 
                     par = rep(times = (dim(X)[2]+1), x = 1), 
                     hessian = T)
modout(mllognormal)
mllognormal
pars <- mllognormal$par
lengthpar <- dim(X)[2]
mu <- (X%*%pars[1:lengthpar])
sigma <- exp(pars[lengthpar+1])  
ufit <- pLOGNO(y, mu, sigma)
reslognormal <- qnorm(ufit)
qqnorm(reslognormal)
abline(0, 1, col = 'blue')

# 1.2 Weibull distribution ------------------------------------------------------------- 
LLweibull <- function(pars, y, X){
  lengthpar <- dim(X)[2]
  a <- exp(pars[lengthpar+1])  # shape
  b <- exp(X %*% pars[1:lengthpar]) # scale
  logL <- dweibull(y, shape = a, scale = b, log = T)
  return(-sum(logL))
}
mlweibull <- optim(f = LLweibull, 
                 y = y, X = X, 
                 control = list(maxit = 100000),
                 par = rnorm(n = (dim(X)[2]+1), 0, 0.1),
                 hessian = T)
mlweibull
modout(mlweibull)

pars <- mlweibull$par
lengthpar <- dim(X)[2]
a <- exp(pars[lengthpar+1])  # shape
b <- exp(X %*% pars[1:lengthpar]) # scale
ufit <- pweibull(y, shape = a, scale = b)
resweibull <- qnorm(ufit)
qqnorm(resweibull)
abline(0, 1, col = 'blue')

# 1.3 GBII distribution (1 - only in mu) ------------------------------------------------------------- 
LLGB2 <- function(pars, y, X){
  lengthpar <- dim(X)[2]
  sigma <- exp(pars[lengthpar+3])  # shape
  mu <- exp(X %*% pars[1:lengthpar]) # scale
  nu <- exp(pars[lengthpar+1])
  tau <- exp(pars[lengthpar+2])
  logL <- dGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau, log = T)
  return(-sum(logL))
}
mlGB2 <- optim(fn = LLGB2, y = y, X = X, 
             par = c(-1, 1, 0.1, 0.1, 0.05, 0.1), 
             control = list(maxit = 100000),
             hessian = T)
mlGB2
mlGB2$hessian
modout(mlGB2)

pars <- mlGB2$par
lengthpar <- dim(X)[2]
sigma <- exp(pars[lengthpar+3])  # shape
mu <- exp(X %*% pars[1:lengthpar]) # scale
nu <- exp(pars[lengthpar+1])
tau <- exp(pars[lengthpar+2])
ufit <- pGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
resGB2 <- qnorm(ufit)
qqnorm(resGB2)
abline(0, 1, col = 'blue')

# 1.4 GBII distribution (2) - in mu and sigma ------------------------------------------------------------- 
X1 <- model.matrix(~ 1, data = dtnew); 
X2 <- model.matrix(~ magnitude + intensity, data = dtnew)
pars <-  c(-0.5, 1.5, 0.1, 0.1, 0.1, 0.05, 1)
LLGB2.sigma <- function(pars, y, X1, X2){
  lengthpar1 <- dim(X1)[2]
  lengthpar2 <- dim(X2)[2]
  mu <- exp(X1 %*% pars[1:lengthpar1]) # scale
  sigma <- X2 %*% pars[(lengthpar1 + 1):(lengthpar1 + lengthpar2)]  # shape
  nu <- exp(pars[lengthpar1 + lengthpar2 + 1])
  tau <- exp(pars[lengthpar1 + lengthpar2 + 2])
  logL <- -lbeta(nu, tau) + log(abs(sigma)) + sigma*nu*log(y/mu) - log(y) - (nu + tau)*log(1 + (y/mu)^sigma)
  return(-sum(logL))
}
mlGB2.sigma <- optim(fn = LLGB2.sigma, y = y, X1 = X1, X2 = X2, 
               par = c(5, 0.1, 0.1, 0.1, 0.1, 0.1),
               control = list(maxit = 100000),
               hessian = T)
mlGB2.sigma
mlGB2.sigma$hessian
modout(mlGB2.sigma)

pars <- mlGB2.sigma$par
lengthpar1 <- dim(X1)[2]
lengthpar2 <- dim(X2)[2]
mu <- exp(X1 %*% pars[1:lengthpar1]) # scale
sigma <- X2 %*% pars[(lengthpar1 + 1):(lengthpar1 + lengthpar2)]  # shape
nu <- exp(pars[lengthpar1 + lengthpar2 + 1])
tau <- exp(pars[lengthpar1 + lengthpar2 + 2])
ufit <- pGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
resGB2.sigma <- qnorm(ufit)
qqnorm(resGB2.sigma)
abline(0, 1, col = 'blue')

# 1.5 GBII distribution (3) - in mu and nu ------------------------------------------------------------- 
X1 <- model.matrix(~ 1, data = dtnew); 
X2 <- model.matrix(~ intensity + magnitude, data = dtnew)
pars <-  c(-0.5, 1.5, 0.1, 0.1, 0.1, 0.05)
LLGB2.nu <- function(pars, y, X1, X2){
  lengthpar1 <- dim(X1)[2]
  lengthpar2 <- dim(X2)[2]
  mu <- exp(X1 %*% pars[1:lengthpar1])                     # scale
  sigma <- (pars[lengthpar1 + 1])                # shape
  nu <- exp(X2 %*% pars[(lengthpar1 + 2):(lengthpar1 + lengthpar2 + 1)])
  tau <- exp(pars[lengthpar1 + lengthpar2 + 2])
  logL <- -lbeta(nu, tau) + log(abs(sigma)) + sigma*nu*log(y/mu) - log(y) - (nu + tau)*log(1 + (y/mu)^sigma)
  return(-sum(logL))
}
mlGB2.nu <- optim(fn = LLGB2.nu, 
                  y = y, X1 = X1, X2 = X2, 
                     par = c(1, 0.1, 1, 0.1, 0.1, 0.1), 
                     control = list(maxit = 100000),
                     hessian = T)
mlGB2.nu
mlGB2.nu$hessian
modout(mlGB2.nu)

pars <- mlGB2.nu$par
lengthpar1 <- dim(X1)[2]
lengthpar2 <- dim(X2)[2]
mu <- exp(X1 %*% pars[1:lengthpar1])                     # scale
sigma <- (pars[lengthpar1 + 1])                # shape
nu <- exp(X2 %*% pars[(lengthpar1 + 2):(lengthpar1 + lengthpar2 + 1)])
tau <- exp(pars[lengthpar1 + lengthpar2 + 2])

ufit <- pGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
resGB2.nu <- qnorm(ufit)
qqnorm(resGB2.nu)
abline(0, 1, col = 'blue')

# 1.6 GBII distribution (4) - in mu and tau ------------------------------------------------------------- 
X1 <- model.matrix(~ 1, data = dtnew); 
X2 <- model.matrix(~ intensity + magnitude, data = dtnew)
pars <-  c(-0.5, 1.5, 0.1, 0.1, 0.1, 0.05)
LLGB2.tau <- function(pars, y, X1, X2){
  lengthpar1 <- dim(X1)[2]
  lengthpar2 <- dim(X2)[2]
  mu <- exp(X1 %*% pars[1:lengthpar1])                     # scale
  sigma <- (pars[lengthpar1 + 1])                # shape
  nu <- exp(pars[lengthpar1 + 2])
  tau <- exp(X2 %*% pars[(lengthpar1 + 3):(lengthpar1 + lengthpar2 + 2)])
  logL <- -lbeta(nu, tau) + log(abs(sigma)) + sigma*nu*log(y/mu) - log(y) - (nu + tau)*log(1 + (y/mu)^sigma)
  return(-sum(logL))
}
mlGB2.tau <- optim(fn = LLGB2.tau, 
                  y = y, X1 = X1, X2 = X2, 
                  par = c(1, 0.1, 1, 0.1, 0.1, 0.1), 
                  control = list(maxit = 100000),
                  hessian = T)
mlGB2.tau
mlGB2.tau$hessian
modout(mlGB2.tau)

pars <- mlGB2.tau$par
lengthpar1 <- dim(X1)[2]
lengthpar2 <- dim(X2)[2]
mu <- exp(X1 %*% pars[1:lengthpar1])                     # scale
sigma <- (pars[lengthpar1 + 1])                # shape
nu <- exp(pars[lengthpar1 + 2])
tau <- exp(X2 %*% pars[(lengthpar1 + 3):(lengthpar1 + lengthpar2 + 2)])

ufit <- pGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
resGB2.tau <- qnorm(ufit)
qqnorm(resGB2.tau)
abline(0, 1, col = 'blue')

# 1.7 log moyal - gamma distribution (1) ------------------------------------------------------------- 
LLlogmoyalGA <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  sigma <- exp(X %*% pars[1:lengthpar])
  b <- exp(pars[lengthpar+1])
  kappa <- exp(pars[lengthpar+2])
  ll <- -0.5*log(2*pi) - log(sigma) + kappa*log(b) + lgamma(kappa + 0.5) - lgamma(kappa) - (1/(2*sigma)+1)*log(y) - (kappa + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  
  loglike <- -sum(ll)  
  return(loglike)
}
pars <- c(1,1)
mlogmoyalGA <- optim(fn = LLlogmoyalGA, X = X,
                     par = c(-1, 0.3,0.1, 1, -1.3),
                     control = list(maxit = 100000),
                     y = y, hessian = T)
modout(mlogmoyalGA)
mlogmoyalGA

pars <- mlogmoyalGA$par
lengthpar <- dim(X)[2]
sigma <- exp(X %*% pars[1:lengthpar])
b <- exp(pars[lengthpar+1])
a <- exp(pars[lengthpar+2])
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pLMGA(y[i], sigma = sigma[i], a = a, b = b)
}
resLMGA <- qnorm(ufit)
qqnorm(resLMGA)
abline(0, 1, col = 'blue')


# 1.8 log moyal - gamma distribution (2) ------------------------------------------------------------- 
LLlogmoyalGA2 <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  b <- exp(X %*%pars[1:lengthpar])
  sigma <- exp(pars[lengthpar+1])
  a <- exp(pars[lengthpar+2])
  ll <- -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  
  loglike <- -sum(ll)  
  return(loglike)
}
pars <- c(1,1)
mlogmoyalGA2 <- optim(par = rep(times = (dim(X)[2]+2), x = 0.1),
                   fn = LLlogmoyalGA2, X = X,
                   control = list(maxit = 100000),
                   y = y, hessian = T)

modout(mlogmoyalGA2)
mlogmoyalGA2

pars <- mlogmoyalGA2$par
lengthpar <- dim(X)[2]
b <- exp(X %*%pars[1:lengthpar])
sigma <- exp(pars[lengthpar+1])
a <- exp(pars[lengthpar+2])
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pLMGA(y[i], sigma = sigma, a = a, b = b[i])
}
resLMGA2 <- qnorm(ufit)
qqnorm(resLMGA2)
abline(0, 1, col = 'blue')


# 1.9 log moyal - gamma distribution (3) ------------------------------------------------------------- 
Xsigma <- model.matrix(~ magnitude, data = dtnew); 
Xb <- model.matrix(~ intensity, data = dtnew)
LLlogmoyalGA3 <- function(y, pars, Xsigma, Xb) {
  sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
  b <- exp(Xb %*%pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])])
  a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
  ll <- -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  loglike <- -sum(ll)  
  return(loglike)
}
mlogmoyalGA3 <- optim(fn = LLlogmoyalGA3, Xsigma = Xsigma,
                     Xb = Xb,
                     y = y, hessian = T,
                     control = list(maxit = 50000),
                     method = 'Nelder-Mead',
                     par = c(-1,-0.1,0.5,-2,-1))
modout(mlogmoyalGA3)
mlogmoyalGA3

pars <- mlogmoyalGA3$par
sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
b <- exp(Xb %*%pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])])
a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pLMGA(y[i], sigma = sigma[i], a = a, b = b[i])
}
resLMGA3 <- qnorm(ufit)
qqnorm(resLMGA3)
abline(0, 1, col = 'blue')

# 2.0 log moyal distribution ------------------------------------------------------------- 
LLlogmoyal <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  sigma <- exp(X %*% pars[1:lengthpar])
  tau <- exp(pars[lengthpar+1])
  ll <- 0.5*log(tau) - 0.5*log(2*pi) - log(sigma) - (1/(2*sigma) + 1)*log(y) - tau/2*(1/y)^(1/sigma)
    
  loglike <- -sum(ll)  
  return(loglike)
}
mlogmoyal <- optim(fn = LLlogmoyal, X = X,
                    par = c(0.1, 0.3, 0.1, 0.1),
                    y = y, hessian = T)
modout(mlogmoyal)
mlogmoyal

pars <- mlogmoyal$par
lengthpar <- dim(X)[2]
sigma <- exp(X %*% pars[1:lengthpar])
tau <- exp(pars[lengthpar+1])
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pGlogM(y[i], sigma = sigma[i], tau = tau)
}
resGlogM <- qnorm(ufit)
qqnorm(resGlogM)
abline(0, 1, col = 'blue')

# 2.1 exponentiated Frechet regression ------------------------------------------------------------- 
LLeFrechet <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  lambda <- exp(X %*% pars[1:lengthpar])
  tau <- exp(pars[lengthpar+1])
  alpha <- exp(pars[lengthpar+2])
  
  ll <- log(alpha) + log(tau) + log(lambda) + (alpha - 1)*log(1 - exp(-tau/y^lambda)) - (1+lambda)*log(y) - tau/y^lambda
  loglike <- -sum(ll)  
  return(loglike)
}
meFrechet <- optim(fn = LLeFrechet, X = X,
                   par = c(-0.8,-0.2,0.1,2.5,0),
                   control = list(maxit = 50000),
                   y = y, hessian = T)
modout(meFrechet)
meFrechet
pars <- meFrechet$par
lengthpar <- dim(X)[2]
lambda <- exp(X %*% pars[1:lengthpar])
tau <- exp(pars[lengthpar+1])
alpha <- exp(pars[lengthpar+2])
ufit <- c()
peFrechet <- function(y, lambda, tau, alpha){
  sigma <- tau^(1/lambda)
  1 - (1 - exp(-(sigma/y)^lambda))^alpha
  }
for(i in 1:length(y)){
  ufit[i] <- peFrechet(y[i], lambda = lambda[i], tau = tau, alpha = alpha)
}
ufit
reseFrechet <- qnorm(ufit)
qqnorm(reseFrechet)
abline(0, 1, col = 'blue')


# 2.2 burr regression ------------------------------------------------------------- 
LLburr <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  tau <- exp(X %*% pars[1:lengthpar])
  lambda <- exp(pars[lengthpar+1])
  beta <- exp(pars[lengthpar+2])
  
  ll <- log(lambda) + lambda*log(beta) + log(tau) + (tau-1)*log(y) - (lambda+1)*log(beta+y^tau)
  loglike <- -sum(ll)  
  return(loglike)
}
mburr <- optim(fn = LLburr, X = X,
               par = c(2, -1, -1, -1, 5),    
               y = y, hessian = T)
modout(mburr)
mburr
pars <- mburr$par
lengthpar <- dim(X)[2]
tau <- exp(X %*% pars[1:lengthpar])
lambda <- exp(pars[lengthpar+1])
beta <- exp(pars[lengthpar+2])
ufit <- c()
pburrnew <- function(y, lambda, beta, tau){1 - (beta/(beta+y^tau))^lambda}
for(i in 1:length(y)){
  ufit[i] <- pburrnew(y[i], lambda = lambda, tau = tau[i], beta = beta)
}
resburr <- qnorm(ufit)
qqnorm(resburr)
abline(0, 1, col = 'blue')

# 2.3 gamma-generalized inverse gaussian regression ------------------------------------------------------------- 

pars <- rep(times = (dim(X)[2] + 4), x = -1)
LLGAGIG <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  t <- exp(X %*% pars[1:lengthpar])
  a <- exp(pars[lengthpar+1])
  lambda <- exp(pars[lengthpar+2])
  mu <- exp(pars[lengthpar+3])
  phi <- exp(pars[lengthpar+4])
  
  munew <- mu*sqrt(phi/(phi+2*mu^2*t*y))
  K1 <- besselK(x = phi/munew, nu = (a+lambda))
  K2 <- besselK(x = phi/mu, nu = (lambda))
  
  ll <- a*(log(y) + log(t) + log(mu)) - log(y) - lgamma(a) + (a+lambda)/2*(log(phi) - log(phi+2*mu^2*t*y)) + log(K1) - log(K2)
  loglike <- -sum(ll)  
  return(loglike)
}

mGAGIG <- optim(fn = LLGAGIG, X = X,
               #par = c(1,-1,-1,1,-1,1,1),
               par = c(1, -1,-1,1,-1,2,2),
               control = list(maxit = 50000),
               y = y, hessian = T)
modout(mGAGIG)
mGAGIG
pars <- mGAGIG$par
lengthpar <- dim(X)[2]
t <- exp(X %*% pars[1:lengthpar])
a <- exp(pars[lengthpar+1])
lambda <- exp(pars[lengthpar+2])
mu <- exp(pars[lengthpar+3])
phi <- exp(pars[lengthpar+4])
pGAGIG <- function(y, t, a, lambda, mu){
  pGAGIG <- function(y, t , a, lambda, mu){
    munew <- mu*sqrt(phi/(phi+2*mu^2*t*y))
    K1 <- besselK(x = phi/munew, nu = (a+lambda))
    K2 <- besselK(x = phi/mu, nu = (lambda))
    ll <- a*(log(y) + log(t) + log(mu)) - log(y) - lgamma(a) + (a+lambda)/2*(log(phi) - log(phi+2*mu^2*t*y)) + log(K1) - log(K2)
    exp(ll)
  }
  z <- c()
  for (i in 1:length(y)){
    pnew <- function(y){
      pGAGIG(y, t = t[i], a = a, lambda = lambda, mu = mu)
    }
    z[i] = integrate(pnew, lower = 0, upper = y[i])$value
  }
  z
}

ufit <- c()
ufit <- pGAGIG(y = y, t = t, a = a, lambda = lambda, mu = mu)
resGAGIG <- qnorm(ufit)
qqnorm(resGAGIG)
abline(0, 1, col = 'blue')

# 2.4 gamma regression distribution ------------------------------------------------------------- 
LLgamma <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  mu <- exp(X %*% pars[1:lengthpar])
  sigma <- exp(pars[lengthpar+1])

  ll <- dGA(y, mu = mu, sigma = sigma, log = T)
  loglike <- -sum(ll)
  return(loglike)
}
mgamma <- optim(fn = LLgamma, X = X,
                par = rep(times = (dim(X)[2]+1), x = 1),
                control = list(maxit = 50000),
                y = y, hessian = T)
modout(mgamma)
mgamma
pars <- mgamma$par
lengthpar <- dim(X)[2]
mu <- exp(X %*% pars[1:lengthpar])
sigma <- exp(pars[lengthpar+1])
ufit <- c()
ufit <- pGA(y, mu = mu, sigma = sigma)
resgamma <- qnorm(ufit)
qqnorm(resgamma)
abline(0, 1, col = 'blue')

# 2.5 inverse guaissn regression distribution ------------------------------------------------------------- 
LLigua <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  mu <- exp(X %*% pars[1:lengthpar])
  sigma <- exp(pars[lengthpar+1])
  
  ll <- dIG(y, mu = mu, sigma = sigma, log = T)
  loglike <- -sum(ll)
  return(loglike)
}
migua <- optim(fn = LLigua, X = X,
                #par = rep(times = (dim(X)[2]+1), x = 0.1),
                par = c(-5,1,1,1),
                control = list(maxit = 50000),
                y = y, hessian = T)
modout(migua)
migua
pars <- migua$par
lengthpar <- dim(X)[2]
mu <- exp(X %*% pars[1:lengthpar])
sigma <- exp(pars[lengthpar+1])
ufit <- c()
ufit <- pIG(y, mu = mu, sigma = sigma)
resigua <- qnorm(ufit)
qqnorm(resigua)
abline(0, 1, col = 'blue')

# 2.6 exp - inverse guaissn regression distribution ------------------------------------------------------------- 
LLexpigua <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  t <- exp(X %*% pars[1:lengthpar])
  delta <- exp(pars[lengthpar+1])
  phi <- (delta^2+2*y/t)^(0.5)
  
  ll <- log(delta)-log(t)-delta*phi-delta^2-3*log(phi)+log(delta*phi+1)
  loglike <- -sum(ll)
  return(loglike)
}
mexpigua <- optim(fn = LLexpigua, X = X,
               #par = rep(times = (dim(X)[2]+1), x = 0.1),
               par = c(-1,1,1,1),
               control = list(maxit = 50000),
               y = y, hessian = T)
modout(mexpigua)
mexpigua
pexpigua <- function(y, delta, t){
  phi <- (delta^2+2*y/t)^(0.5)
  1 - delta*exp(-delta*(phi-delta))/phi 
}
pars <- mexpigua$par
lengthpar <- dim(X)[2]
t <- exp(X %*% pars[1:lengthpar])
delta <- exp(pars[lengthpar+1])
phi <- (delta^2+2*y/t)^(0.5)
ufit <- c()
ufit <- pexpigua(y, delta = delta, t = t)
resexpigua <- qnorm(ufit)
qqnorm(resexpigua)
abline(0, 1, col = 'blue')

# comparsion of all models ------------------------------------------------------------- 
npars.com <- rbind(nrow(modout(mllognormal)$summary),
                   nrow(modout(mlweibull)$summary),
                   nrow(modout(mlGB2)$summary), 
                   nrow(modout(mlogmoyalGA)$summary),
                   nrow(modout(mlogmoyalGA2)$summary),
                   nrow(modout(mlogmoyalGA3)$summary),
                   nrow(modout(mlogmoyal)$summary),
                   nrow(modout(meFrechet)$summary),
                   nrow(modout(mburr)$summary),
                   nrow(modout(mGAGIG)$summary),
                   nrow(modout(mgamma)$summary),
                   nrow(modout(migua)$summary),
                   nrow(modout(mexpigua)$summary),
                   nrow(modout(mlGPD)$summary),
                   nrow(modout(mlGB2.sigma)$summary),
                   nrow(modout(mlGB2.nu)$summary),
                   nrow(modout(mlGB2.tau)$summary))

loglike.com <- rbind(modout(mllognormal)$ll,
                     modout(mlweibull)$ll,
                     modout(mlGB2)$ll, 
                     modout(mlogmoyalGA)$ll,
                     modout(mlogmoyalGA2)$ll,
                     modout(mlogmoyalGA3)$ll,
                     modout(mlogmoyal)$ll,
                     modout(meFrechet)$ll,
                     modout(mburr)$ll,
                     modout(mGAGIG)$ll,
                     modout(mgamma)$ll,
                     modout(migua)$ll,
                     modout(mexpigua)$ll,
                     modout(mlGPD)$ll,
                     modout(mlGB2.sigma)$ll,
                     modout(mlGB2.nu)$ll,
                     modout(mlGB2.tau)$ll)

AIC.com <- rbind(modout(mllognormal)$AIC,
                 modout(mlweibull)$AIC,
                 modout(mlGB2)$AIC, 
                 modout(mlogmoyalGA)$AIC,
                 modout(mlogmoyalGA2)$AIC,
                 modout(mlogmoyalGA3)$AIC,
                 modout(mlogmoyal)$AIC,
                 modout(meFrechet)$AIC,
                 modout(mburr)$AIC,
                 modout(mGAGIG)$AIC,
                 modout(mgamma)$AIC,
                 modout(migua)$AIC,
                 modout(mexpigua)$AIC,
                 modout(mlGPD)$AIC,
                 modout(mlGB2.sigma)$AIC,
                 modout(mlGB2.nu)$AIC,
                 modout(mlGB2.tau)$AIC)

BIC.com <- rbind(modout(mllognormal)$BIC,
                 modout(mlweibull)$BIC,
                 modout(mlGB2)$BIC, 
                 modout(mlogmoyalGA)$BIC,
                 modout(mlogmoyalGA2)$BIC,
                 modout(mlogmoyalGA3)$BIC,
                 modout(mlogmoyal)$BIC,
                 modout(meFrechet)$BIC,
                 modout(mburr)$BIC,
                 modout(mGAGIG)$BIC,
                 modout(mgamma)$BIC,
                 modout(migua)$BIC,
                 modout(mexpigua)$BIC,
                 modout(mlGPD)$BIC,
                 modout(mlGB2.sigma)$BIC,
                 modout(mlGB2.nu)$BIC,
                 modout(mlGB2.tau)$BIC)

convergence.com <- rbind((mllognormal)$convergence,
                         (mlweibull)$convergence,
                         (mlGB2)$convergence, 
                         (mlogmoyalGA)$convergence,
                         (mlogmoyalGA2)$convergence,
                         (mlogmoyalGA3)$convergence,
                         (mlogmoyal)$convergence,
                         (meFrechet)$convergence,
                         (mburr)$convergence,
                         (mGAGIG)$convergence,
                         (mgamma)$convergence,
                         (migua)$convergence,
                         (mexpigua)$convergence,
                         (mlGPD)$convergence,
                         (mlGB2.sigma)$convergence,
                         (mlGB2.nu)$convergence,
                         (mlGB2.tau)$convergence)
  

row.names(AIC.com) <- c('lognormal', 'weibull', 'GB2',
                        'logmoyalGA','logmoyalGA2','logmoyalGA3',
                        'GlogM', 'eFrechet', 'burr', 'GAGIG', 
                        'gamma', 'ig', 'expigua', 'GP',
                        'GB2.sigma','GB2.nu','GB2.tau')
row.names(BIC.com) <- c('lognormal', 'weibull', 'GB2',
                        'logmoyalGA','logmoyalGA2','logmoyalGA3',
                        'GlogM', 'eFrechet', 'burr', 'GAGIG', 
                        'gamma', 'ig', 'expigua', 'GP',
                        'GB2.sigma','GB2.nu','GB2.tau')
colnames(AIC.com) <- 'AIC'
colnames(BIC.com) <- 'BIC'
colnames(loglike.com) <- 'loglike'
colnames(npars.com) <- 'npars'
outtable <- round(cbind(npars.com, loglike.com, AIC.com, BIC.com), 1)
outtable <- round(cbind(npars.com, loglike.com, AIC.com, BIC.com), 1)

outtable[order(outtable[,3], decreasing = F),]  # AIC sorting
outtable[order(outtable[,4], decreasing = F),]  # BIC sorting


# QQ - plot
par(mfrow = c(2, 2))
qqnorm(resLMGA, main = 'GLMGA I', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resLMGA2, main = 'GLMGA II', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resLMGA3, main = 'GLMGA III', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resGB2, main = 'GB2', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')

