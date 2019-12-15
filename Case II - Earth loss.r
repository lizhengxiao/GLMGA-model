

# 1. Descriptive analysis -------------------------------------------------
library(gamlss)
library(data.table)
library(actuar)
setwd("D:/对外经济贸易大学 - 科研/0. Generalizing the log-Moyal distribution and regression models for heavy tailed loss data/合作文件 - 杨亮 - 文稿+代码")
source('0. definition of logmoyal distribution.r')
source('0. definition of logmoyal-gamma distribution.r')
dtraw <- fread('D:/对外经济贸易大学 - 科研/0. Generalizing the log-Moyal distribution and regression models for heavy tailed loss data/composite models - r/eqdata.csv', header = T)
#dtraw <- fread('D:/坚果云/同步文件夹/对外经济贸易大学 - 科研/29. 地震损失金额和死亡人数的精算模型/合作文件 - 杨亮 - 文稿+代码/eqdata.csv', header = T)
names(dtraw)[21] <- 'y'
names(dtraw)[1] <- 'year'
names(dtraw)[10] <- 'magnitude'
names(dtraw)[11] <- 'intensity'
names(dtraw)[9] <- 'location'
names(dtraw)[8] <- 'time'
names(dtraw)[6] <- 'hour'
names(dtraw)[7] <- 'min'

dtraw

dtgpd <- fread('D:/对外经济贸易大学 - 科研/0. Generalizing the log-Moyal distribution and regression models for heavy tailed loss data/composite models - r/gpd.csv', header = T)
#dtgpd <- fread('D:/坚果云/同步文件夹/对外经济贸易大学 - 科研/29. 地震损失金额和死亡人数的精算模型/composite models - r/gpd.csv', header = T)
gpd <- dtgpd[year %in% seq(1990, 2015)]$`国内生产总值(亿元)`
names(gpd) <- rev(seq(1990, 2015))
dtgpd <- data.table(year = rev(seq(1990, 2015)), gpd)
dtgpd[, gpdfactor := rev(gpd/18872.9)]

setkey(dtgpd, year)
setkey(dtraw, year)

dt <- merge(dtraw, dtgpd)
dt[, ynew := y*gpdfactor]
dt  # contains NA's and 0's

dtnew <- dt[!is.na(ynew)&ynew != 0] # omit 0 and NA in the dataset.
dt[ynew ==0]
dtnew <- dtnew[, .(year, ynew, y,
                   death = 死亡人数, 
                   inj = 受伤人数, magnitude, intensity, location, hour, min, time)]

dtnew$ynew <- dtnew$ynew/100
dtnew$y <- dtnew$y/100
dtnew[, i.f := year - 1989]
dtnew$count <- 1
dtnew[is.na(intensity)]$intensity = median(dtnew$intensity, na.rm = T)
dtnew[, night := ifelse(hour>=1&hour<=6, 1, 0)]

y <- dtnew$ynew
summary(dtnew)
summary(y)
y
X <- model.matrix(~ magnitude + intensity, data = dtnew)
dtnew
#write.csv(dtnew, 'dtnew.csv')

#hist(y, xlab = 'Economic losses (0-50000)', main = '',
#     xlim = c(50000, 3600000), breaks = 100, col = 'gray')
library(fitdistrplus)
m0 <- fitdist(y, distr = 'exp', 
              method = 'mme', 
              start = list(rate = 6.862071e-05))
rate <- m0$estimate
ufit <- pexp(y, rate = rate)
resexp <- qnorm(ufit)
resexp <- resexp[!resexp==Inf]
par(mfrow = c(1,2))
hist(y, xlab = 'Economic losses (0-20000)', main = '',
     xlim = c(0, 20000), breaks = 3000, col = 'gray')
qqnorm(resexp, main = '', ylim = c(-4, 4), xlim = c(-4, 4))
abline(0, 1, col = 'blue')


#hist(y)


library(moments)
dtsum <- dtnew[,.(N = sum(count), ymean = mean(ynew), ym = quantile(ynew, 0.5),
                  ysd = sd(ynew), yskew = skewness(ynew),
                  ykur = kurtosis(ynew)
), by = 'year']
mean(y)
median(y)
sd(y)
skewness(y)
kurtosis(y)
#write.csv(dtsum,'C:/Users/Lee/Desktop/summary.csv')


# base model: GP distribution -------------------------------------------------
library(extRemes)
library(POT)
umin <- 0
GPD2.MLE <- fevd(x = ynew, 
                 scale.fun = ~ 1, 
                 shape.fun = ~ magnitude + intensity, 
                 threshold = 0, method = 'MLE', use.phi = F,
                 type = "GP", data = dtnew)
GPD2.MLE
#plot(GPD2.MLE)
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
GPD2.Z = GPD2.par/GPD2.se                                            ##  Z统计量
GPD2.p = ifelse(GPD2.Z>=0, pnorm(GPD2.Z, lower=F)*2, pnorm(GPD2.Z)*2)    ## p值
data.table(GPD2.par, GPD2.se, GPD2.p)

mlGPD <- GPD2.MLE$results
modout(mlGPD)

# difine the hessian matrix output and the stardard error -------------------------------------------------
modout <- function(m) {
  Hessian <- m$hessian
  se <- sqrt(diag(solve(Hessian)))                  ## 标准误
  Z <- m$par/se                                             ##  Z统计量
  p <- ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)           ## p值
  summarytable <- round(data.frame(m$par, se, Z, p), 3)
  NLL <- m$value  # - loglikelihood value
  AIC <- 2*length(m$par) + 2*NLL
  BIC <- log(length(y))*length(m$par) + 2*NLL
  list(summary = summarytable, ll =  - m$value, AIC = AIC, BIC = BIC)
}

# X <- model.matrix(~ magnitude + intensity, data = dtnew)
# X <- model.matrix(~ magnitude + factor(intensity), data = dtnew)

# 1. log normal distribution -------------------------------------------------
# X <- model.matrix(~ 1, data = dtnew)
LLlognormal <-  function(pars, y, X){
  lengthpar <- dim(X)[2]
  mu <- (X%*%pars[1:lengthpar])
  sigma <- exp(pars[lengthpar+1])  # 注意sigma的估计值应该是结果的exp形式
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
sigma <- exp(pars[lengthpar+1])  # 注意sigma的估计值应该是结果的exp形式
ufit <- pLOGNO(y, mu, sigma)
reslognormal <- qnorm(ufit)
qqnorm(reslognormal)
abline(0, 1, col = 'blue')

# 1.2 Pareto distribution -------------------------------------------------

pars <- rep(times = (dim(X)[2]+1), x = 1)
#X <- model.matrix(~ magnitude + intensity, data = dtnew)
LLpareto <- function(pars, y, X){
  lengthpar <- dim(X)[2]
  beta <- exp(pars[lengthpar+1])  # 同样的
  sigma <- exp(X%*%pars[1:lengthpar]) # 同样的
  logL <- log(beta) + beta*log(sigma) - (beta+1)*log(y+sigma) # 对数似然函数
  ll <- -sum(logL)
  return(ll)
}
mlpareto <- optim(fn = LLpareto, 
                  y = y, X = X, control = list(maxit = 100000),
                  par = rep(times = (dim(X)[2]+1), x = 0.1), hessian = T)
mlpareto
modout(mlpareto)
pars <- mlpareto$par
lengthpar <- dim(X)[2]
beta <- exp(pars[lengthpar+1])  # 同样的
sigma <- exp(X%*%pars[1:lengthpar]) # 同样的
ufit <- ppareto2(y, shape = beta, scale = sigma)
respareto <- qnorm(ufit)
qqnorm(respareto)
abline(0, 1, col = 'blue')

# 1.3 Weibull distribution ------------------------------------------------------------- 
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
                 #par = rep(times = (dim(X)[2]+1), x = 0.1), 
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

# 1.4 GBII distribution (1 - only in mu) ------------------------------------------------------------- 
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
             #par = c(1.5, 0.1, 0.1, 0.1, 0.05, 0.7),
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
               #par = c(-1, 0.1, 1, 0.1, 0.1, 0.1, 0.1), 
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

# 1.4 GBII distribution (3) - in mu and nu ------------------------------------------------------------- 
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
                     #par = c(1.5, 0.1, 0.1, 0.1, 0.05, 0.7),
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

# 1.4 GBII distribution (4) - in mu and tau ------------------------------------------------------------- 
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
                  #par = c(1.5, 0.1, 0.1, 0.1, 0.05, 0.7),
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

# 1.5 log moyal - gamma distribution (1) ------------------------------------------------------------- 
LLlogmoyalGA <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  sigma <- exp(X %*% pars[1:lengthpar])
  b <- exp(pars[lengthpar+1])
  kappa <- exp(pars[lengthpar+2])
  #ll <- -0.5*log(2) - log(sigma) + 0.5*log(b) - lbeta(0.5, 0.5) - (1/(2*sigma)+1)*log(y) - (0.5 + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  ll <- -0.5*log(2*pi) - log(sigma) + kappa*log(b) + lgamma(kappa + 0.5) - lgamma(kappa) - (1/(2*sigma)+1)*log(y) - (kappa + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  
  loglike <- -sum(ll)  
  return(loglike)
}
pars <- c(1,1)
mlogmoyalGA <- optim(fn = LLlogmoyalGA, X = X,
                     #par = c(-3, 0.3, -11, -1.3),
                     #par = c(-3, 0.3,0.1, -11, -1.3),
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


# 1.6 log moyal - gamma distribution (2) ------------------------------------------------------------- 
LLlogmoyalGA2 <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  b <- exp(X %*%pars[1:lengthpar])
  sigma <- exp(pars[lengthpar+1])
  a <- exp(pars[lengthpar+2])
  #ll <- -0.5*log(2) - log(sigma) + 0.5*log(b) - lbeta(0.5, 0.5) - (1/(2*sigma)+1)*log(y) - (0.5 + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  ll <- -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  
  loglike <- -sum(ll)  
  return(loglike)
}
pars <- c(1,1)
mlogmoyalGA2 <- optim(par = rep(times = (dim(X)[2]+2), x = 0.1),
                   fn = LLlogmoyalGA2, X = X,
                   control = list(maxit = 100000),
                   y = y, hessian = T)
# mlogmoyalGA2 <- nlm(p = c(-3, -1, -1.3, 0, 0),
#                     f = LLlogmoyalGA2, X = X,
#                     y = y, hessian = T)
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


# 1.7 log moyal - gamma distribution (3) ------------------------------------------------------------- 

# Xsigma <- model.matrix(~ magnitude + intensity, data = dtnew);
# Xb <- model.matrix(~ intensity, data = dtnew)
Xsigma <- model.matrix(~ magnitude, data = dtnew); 
Xb <- model.matrix(~ intensity, data = dtnew)
#Xb <- model.matrix(~ 1, data = dt0)
LLlogmoyalGA3 <- function(y, pars, Xsigma, Xb) {
  sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]])
  b <- exp(Xb %*%pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])])
  a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
  #ll <- -0.5*log(2) - log(sigma) + 0.5*log(b) - lbeta(0.5, 0.5) - (1/(2*sigma)+1)*log(y) - (0.5 + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
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



# 1.8 log moyal distribution ------------------------------------------------------------- 

LLlogmoyal <- function(y, pars, X) {
  lengthpar <- dim(X)[2]
  sigma <- exp(X %*% pars[1:lengthpar])
  tau <- exp(pars[lengthpar+1])
  #ll <- -0.5*log(2) - log(sigma) + 0.5*log(b) - lbeta(0.5, 0.5) - (1/(2*sigma)+1)*log(y) - (0.5 + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
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

# 1.9 exponentiated Frechet regression ------------------------------------------------------------- 
#X <- model.matrix(~ magnitude, data = dtnew)

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
                   #par = rep(times = (dim(X)[2]+2), x = -1),
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


# 2.0 burr regression ------------------------------------------------------------- 

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
                   #par = rep(times = (dim(X)[2]+2), x = 0.1),
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

# 2.1 gamma-generalized inverse gaussian regression ------------------------------------------------------------- 

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

# # 2.2 The Log–Lindley distribution -------------------------------------------------------------  

# LLlogLy <- function(y, pars, X) {
#   lengthpar <- dim(X)[2]
#   theta <- exp(X %*% pars[1:lengthpar])
#   sigma <- exp(pars[lengthpar+1])
#   
#   ll <- log(lambda) + lambda*log(beta) + log(tau) + (tau-1)*log(y) - (lambda+1)*log(beta+y^tau)
#   loglike <- -sum(ll)  
#   return(loglike)
# }


# 2.2 gamma regression distribution ------------------------------------------------------------- 
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
                #par = c(-2,-1,-1,1,-1,2,2),
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

# 2.2 inverse guaissn regression distribution ------------------------------------------------------------- 
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

# 2.3 exp - inverse guaissn regression distribution ------------------------------------------------------------- 
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


# 2.4 GB2 speical case ------------------------------------------------------------- 
LLGB2s <- function(pars, X, y){
  lengthpar <- dim(X)[2]
  mu <- exp(X %*% pars[1:lengthpar])
  m <- -exp(pars[lengthpar+1])
  nu <- 0.5
  tau <- exp(pars[lengthpar+2])
  
  ll <- dGB2(y, mu = mu, sigma = m, nu = nu, tau = tau, log = T)
  loglike <- -sum(ll)
  return(loglike)
}
mGB2s <- optim(fn = LLGB2s, X = X,
                par = c(-1,1,1,1,1),
                control = list(maxit = 50000),
                y = y, hessian = T)
modout(mGB2s)
mGB2s
lengthpar <- dim(X)[2]
pars <- mGB2s$par
mu <- exp(X %*% pars[1:lengthpar])
m <- -exp(pars[lengthpar+1])
nu <- 0.5
tau <- exp(pars[lengthpar+2])
ufit <- c()
ufit <- pGB2(y, mu = mu, sigma = m, nu = nu, tau = tau)
resmGB2s <- qnorm(ufit)
qqnorm(resmGB2s)
abline(0, 1, col = 'blue')

# 结果比较 ------------------------------------------------------------- 
npars.com <- rbind(nrow(modout(mllognormal)$summary),
                   nrow(modout(mlpareto)$summary),
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
                   nrow(modout(mGB2s)$summary),
                   nrow(modout(mlGB2.sigma)$summary),
                   nrow(modout(mlGB2.nu)$summary),
                   nrow(modout(mlGB2.tau)$summary))

loglike.com <- rbind(modout(mllognormal)$ll,
                     modout(mlpareto)$ll,
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
                     modout(mGB2s)$ll,
                     modout(mlGB2.sigma)$ll,
                     modout(mlGB2.nu)$ll,
                     modout(mlGB2.tau)$ll)

AIC.com <- rbind(modout(mllognormal)$AIC,
                 modout(mlpareto)$AIC,
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
                 modout(mGB2s)$AIC,
                 modout(mlGB2.sigma)$AIC,
                 modout(mlGB2.nu)$AIC,
                 modout(mlGB2.tau)$AIC)

BIC.com <- rbind(modout(mllognormal)$BIC,
                 modout(mlpareto)$BIC,
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
                 modout(mGB2s)$BIC,
                 modout(mlGB2.sigma)$BIC,
                 modout(mlGB2.nu)$BIC,
                 modout(mlGB2.tau)$BIC)
convergence.com <- rbind((mllognormal)$convergence,
                         (mlpareto)$convergence,
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
                         (mGB2s)$convergence,
                         (mlGB2.sigma)$convergence,
                         (mlGB2.nu)$convergence,
                         (mlGB2.tau)$convergence)
  

row.names(AIC.com) <- c('lognormal', 'pareto', 'weibull', 'GB2',
                        'logmoyalGA','logmoyalGA2','logmoyalGA3',
                        'GlogM', 'eFrechet', 'burr', 'GAGIG', 
                        'gamma', 'ig', 'expigua', 'GP', 'GB2s',
                        'GB2.sigma','GB2.nu','GB2.tau')
row.names(BIC.com) <- c('lognormal', 'pareto', 'weibull', 'GB2',
                        'logmoyalGA','logmoyalGA2','logmoyalGA3',
                        'GlogM', 'eFrechet', 'burr', 'GAGIG', 
                        'gamma', 'ig','expigua', 'GP', 'GB2s',
                        'GB2.sigma','GB2.nu','GB2.tau')
colnames(AIC.com) <- 'AIC'
colnames(BIC.com) <- 'BIC'
colnames(loglike.com) <- 'loglike'
colnames(npars.com) <- 'npars'
outtable <- round(cbind(npars.com, loglike.com, AIC.com, BIC.com), 1)
outtable <- round(cbind(npars.com, loglike.com, AIC.com, BIC.com), 1)

outtable[order(outtable[,3], decreasing = F),]  # AIC sorting
outtable[order(outtable[,4], decreasing = F),]  # BIC sorting

est.table <- cbind(
                rbind(as.matrix(modout(mllognormal)$summary), matrix(NA, nrow = 7-nrow(modout(mllognormal)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlpareto)$summary), matrix(NA, nrow = 7-nrow(modout(mlpareto)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlweibull)$summary), matrix(NA, nrow = 7-nrow(modout(mlweibull)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlGB2)$summary), matrix(NA, nrow = 7-nrow(modout(mlGB2)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlogmoyalGA)$summary), matrix(NA, nrow = 7-nrow(modout(mlogmoyalGA)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlogmoyalGA2)$summary), matrix(NA, nrow = 7-nrow(modout(mlogmoyalGA2)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlogmoyalGA3)$summary), matrix(NA, nrow = 7-nrow(modout(mlogmoyalGA3)$summary), ncol = 4)),
                rbind(as.matrix(modout(mlogmoyal)$summary), matrix(NA, nrow = 7-nrow(modout(mlogmoyal)$summary), ncol = 4)),
                rbind(as.matrix(modout(meFrechet)$summary), matrix(NA, nrow = 7-nrow(modout(meFrechet)$summary), ncol = 4)),
                rbind(as.matrix(modout(mburr)$summary), matrix(NA, nrow = 7-nrow(modout(mburr)$summary), ncol = 4)),
                rbind(as.matrix(modout(mGAGIG)$summary), matrix(NA, nrow = 7-nrow(modout(mGAGIG)$summary), ncol = 4))
                )
round(est.table, 2)
#write.csv(round(est.table,2), 'C:/Users/Lee/Desktop/esttable.csv')
        
# QQ - plot
par(mfrow = c(2, 2))                                
q1 <- qqnorm(resLMGA, main = 'GLMGA I', xlim = c(-3, 3), ylim = c(-3, 3))
abline(0, 1, col = 'red')
q2 <- qqnorm(resLMGA2, main = 'GLMGA II', xlim = c(-3, 3), ylim = c(-3, 3))
abline(0, 1, col = 'red')
q3 <- qqnorm(resLMGA3, main = 'GLMGA III', xlim = c(-3, 3), ylim = c(-3, 3))
abline(0, 1, col = 'red')
q4 <- qqnorm(resGB2, main = 'GB2', xlim = c(-3, 3), ylim = c(-3, 3))
abline(0, 1, col = 'red')

cor(q1$data, q1$qnorm); # GLMGA I
cor(q2$data, q2$qnorm); # GLMGA II
cor(q3$data, q3$qnorm); # GLMGA III
cor(q4$data, q4$qnorm); # GB2

# QQ - plot
par(mfrow = c(4, 3))
qqnorm(reslognormal, main = 'lognormal', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(respareto, main = 'pareto', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resweibull, main = 'weibull', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resGB2, main = 'GB2', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resLMGA, main = 'logmoyal-gamma', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resLMGA2, main = 'logmoyal-gamma2', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resLMGA3, main = 'logmoyal-gamma3', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resGlogM, main = 'logmoyal', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(reseFrechet, main = 'eFrechet', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resburr, main = 'burr', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resGAGIG, main = 'GAGIG', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')
qqnorm(resGB2.sigma, main = 'GB2.sigma', xlim = c(-4, 4), ylim = c(-4, 4))
abline(0, 1, col = 'red')



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


# qqnorm(resLMIG2, main = 'logmoyal-inverse guassian(2)', xlim = c(-4, 4), ylim = c(-4, 4))
# abline(0, 1, col = 'red')
# qqnorm(resLMGIG, main = 'logmoyal-generalized inverse guassian', xlim = c(-4, 4), ylim = c(-4, 4))
# abline(0, 1, col = 'red')

# # =================================================================
# # 2. 组合分布
# # ==================================================================
# # 2.1 Weibull–Pareto model proposed by Scollnik and Sun (2012).
# remove.naninf <- function(x) {
#   x[!is.nan(x) & is.finite(x)] # 去除 NA 的值 和 无限值 
# }
# # design matrix
# X <- model.matrix(~ magnitude, data = dtnew)
# 
# lg.weibull.pareto <- function(y, par, X){
#   #coef <- par[grep('coef', names(par))]  
#   coef <- par[4:5]  
#   alpha <- exp(par[1])
#   theta <- exp(par[2]) # is theta > 0 or in real?
#   beta <- exp(par[3])
#   sigma <- exp(X %*% coef) 
#   # Bakar et al.(2015) 引用 Scollnik and Sun (2012) 论文
#   # 符号转换的时候出错
#   f <- NA
#   for (i in 1:length(y)) {
#     lambda = theta/((beta*theta-sigma[i])/(alpha*(sigma[i]+theta))+1)^(1/alpha)
#     r = (beta/alpha)/(beta/alpha+(sigma[i]+theta)/(lambda*(exp(theta/lambda)-1)))
#     if (y[i] <= theta) {
#       f[i] <- log(r*dweibull(y[i],alpha,lambda)/pweibull(theta,alpha,lambda))
#     }else{
#       f[i] <- log((1-r)*actuar::dpareto(y[i], shape = beta, scale = sigma[i])/(1 - actuar::ppareto(theta, shape = beta, scale = sigma[i])))
#     }
#   }
#   #f <- remove.naninf(f) # loglikelihood function at each observation
#   loglike <- -sum(f)
#   return(loglike)
# }
# starts <- c(theta = 1,
#             beta = 2,
#             sigma = 1,
#             coef = c(1,1))
# par <- starts
# y <- dtnew$ynew
# 
# mlweibull.pareto <- nlm(lg.weibull.pareto, y = dtnew$ynew, 
#                         X = model.matrix(~ magnitude, data = dtnew),
#                         p = starts, 
#                         #p = starts.mat[index.max, ],
#                         #fscale = 4000,
#                         #typsize = c(10, 10, 1, 1, 1),
#                         print.level = 2,
#                         stepmax = c(1, 1,1,1,1), steptol = 10^-5,
#                         iterlim = 5000, hessian = T)
# mlweibull.pareto #　补充相关结果
# Hessian <- mlweibull.pareto$hessian
# se <- sqrt(diag(solve(Hessian)))                  ## 标准误
# Z <- mlweibull.pareto$estimate/se                                             ##  Z统计量
# p <- ifelse(Z>=0, pnorm(Z, lower=F)*2, pnorm(Z)*2)           ## p值
# (summarytable <- round(data.frame(mlweibull.pareto$estimate, se, Z, p), 3))
# list(summary = summarytable, ll =  - mlweibull.pareto$minimum)
# 
# NLLweibull.pareto <- mlweibull.pareto$minimum  # - loglikelihood value
# AICweibull.pareto <- 2*length(mlweibull.pareto$estimate) + 2*NLLweibull.pareto
# BICweibull.pareto <- log(length(y))*length(mlweibull.pareto$estimate) + 2*NLLweibull.pareto
# c(NLLweibull.pareto, AICweibull.pareto, BICweibull.pareto)
# 
# 
# 
# est <- mlweibull.pareto$estimate
# coef <- est[1:2]  
# alpha <- X %*% coef
# theta <- est[3]
# beta <- exp(est[4]) + 0.0001
# sigma <- exp(est[5]) + 0.0001
# 
# # 组合分布的比例
# lambda = theta/((beta*theta-sigma)/(alpha*(sigma+theta))+1)^(1/alpha)
# r = (beta/alpha)/(beta/alpha+(sigma+theta)/(lambda*(exp(theta/lambda)-1)))
# summary(lambda);hist(lambda)
# summary(r); hist(r)
# 
# # ===================================================================
# # trying a set of the initial values of parameters
# # ===================================================================
# set.seed(122)
# ntry <- 20
# starts.mat <- matrix(NA, nrow = ntry, ncol = 5)
# m.list <- list()
# loglike.vector <- NA
# for(i in 1:ntry){
#   starts <- c(coef = rnorm(2, 0, 1),
#               theta = rnorm(1, 0, 1),
#               beta = rnorm(1, 0, 1),
#               sigma = rnorm(1, 0, 1))
#   starts.mat[i,] <- starts
#   m.list <- try(nlm(lg.weibull.pareto, y = dtnew$ynew, 
#                                              X = model.matrix(~ magnitude, data = dtnew),
#                                              p = starts, 
#                                              iterlim = 5000, hessian = F), silent = T)
# 
#   loglike.vector[i] <- as.numeric(try(- m.list$minimum, silent = T))
# }
# max(loglike.vector[loglike.vector < -1000])
# max.loglike <- max(remove.naninf(loglike.vector))
# index.max <- which(loglike.vector == max.loglike)
# index.max <- which(loglike.vector == max(loglike.vector[loglike.vector < -1000]))
# starts.mat[index.max, ]
# 
# 
# 
# 
# 
# mlweibull.pareto
# 
# 
# 
# 
# 
# # ----------------------------------------------------------------
# # 2.2 lognormal–Pareto model proposed by Scollnik and Sun (2007).
# # lognormal-pareto type II
# par <- c(0.1, 0.2, 1)
# lg.lognormal.pareto <- function(y, par){
#   theta <- par[1]  # 
#   alpha <- exp(par[2])
#   sigma <- exp(par[3])
#   
#   r = ((2*pi)^0.5*alpha*sigma*pnorm(alpha*sigma)*exp(0.5*((alpha*sigma)^2)))/(((2*pi)^0.5*alpha*sigma*pnorm(alpha*sigma)*exp(0.5*((alpha*sigma)^2))) + 1)
#   mu = log(theta) - alpha*sigma^2
#   F1 = pnorm((log(theta) - mu)/sigma)
#   
#   f <- NULL
#   for (i in 1:length(y)) {
#     if(y[i] <= theta){
#       f[i] = log(r) - 0.5*log(2*pi) - log(y[i]) - log(sigma) - 0.5*((log(y[i]) - mu)/sigma)^2 - log(F1);
#       #f[i] = log(r) + dlnorm(y[i], meanlog = mu, sdlog = sigma, log = T) - log(F1);
#     }else{
#       f[i] = log(1-r) + log(alpha) + alpha*log(theta) - (alpha + 1)*log(y[i]);
#       #f[i] = log(1-r) + dpareto1(y[i], shape = alpha, min = theta, log = T);
#     }
#   }
#   #f <- remove.naninf(f) # loglikelihood function at each observation
#   loglike <- - sum(f)
#   return(loglike)
# }
# mllognormal.pareto <- nlm(lg.lognormal.pareto, 
#                           y = y, 
#                           p = c(20, 0.2, 1),
#                           print.level = 2,
#                           stepmax = 10, steptol = 10^-6,
#                           iterlim = 5000, hessian = F)
# 
# mllognormal.pareto <- optim(par = c(100, 4, 1),
#                             fn = lg.lognormal.pareto, 
#                             y = y, 
#                             lower = c(-Inf, 0.0001, 0.0001),
#                             upper = c(Inf, Inf, Inf),
#                             method = 'L-BFGS-B', hessian = F)
# 
# 
# mllognormal.pareto <- optim(par = c(1, 1, 1),
#                             fn = lg.lognormal.pareto, 
#                             y = y, 
#                             #lower = c(-Inf, 0.01, 0.01),
#                             #upper = c(Inf, Inf, Inf),
#                             control = list(trace = 4, REPORT = 1),
#                             method = 'Nelder-Mead', hessian = F)
# 
# 
# mllognormal.pareto
# theta <- mllognormal.pareto$estimate[1]
# alpha <- exp(mllognormal.pareto$estimate[2])
# sigma <- exp(mllognormal.pareto$estimate[3])
# # 组合分布的比例
# r = ((2*pi)^0.5*alpha*sigma*pnorm(alpha*sigma)*exp(0.5*((alpha*sigma)^2)))/(((2*pi)^0.5*alpha*sigma*pnorm(alpha*sigma)*exp(0.5*((alpha*sigma)^2))) + 1)
# mu = log(theta) - alpha*sigma^2
# F1 = pnorm((log(theta) - mu)/sigma)
# # 补充相关的计算结果
# 
# 
# 
# 
# 
# # 2.3 lognormal–bull model 
# par <- c(0.178, 1.039, 0.347, 4.111, 0.841)
# lg.lognormal.bull <- function(y, par){
#   sigma <- exp(par[1])
#   theta <- exp(par[2])
#   alpha <- exp(par[3])
#   beta <- exp(par[4])
#   s <- exp(par[5])
#   
#   mu <- log(theta) - sigma^2*((alpha+1)*beta*theta^beta/(theta^beta + s^beta) - beta)
#   F1 = pnorm((log(theta) - mu)/sigma)
#   F2 <- actuar::pburr(theta, shape1 = alpha, shape2 = beta, scale = s)
#   phi <- ((theta^beta + s^beta)*dnorm((log(theta) - mu)/sigma))/(sigma*alpha*beta*theta^beta*pnorm((log(theta) - mu)/sigma))
#   r <- 1/(1 + phi)
#   
#   f <- NULL
#   for (i in 1:length(y)) {
#     if(y[i] <= theta){
#       #f[i] = log(r) + dlnorm(y[i], meanlog = mu, sdlog = sigma, log = T) - log(F1);
#       f[i] = log(r) - 0.5*log(2*pi) - log(y[i]) - log(sigma) - 0.5*((log(y[i]) - mu)/sigma)^2 - log(F1);
#     }else{
#       f[i] = log(1-r) + log(alpha) + log(beta) + beta*(log(y[i]) - log(s)) - log(y[i]) - (alpha + 1)*log(1 + (y[i]/s)^beta) + (alpha)*log(1 + (theta/s)^beta) 
#       # f[i] = log(1-r) + dburr(y[i], shape1 = alpha, shape2 = beta, scale = s, log = T) - log(1 - F2)
#     }
#   }
#   f <- - sum(f)
#   return(f)
# }
# mllognormal.bull <- nlm(lg.lognormal.bull, y = y, p = c(0.65, .1, .1, .1, .1), iterlim = 1000)
# mllognormal.bull
# # 补充相关的计算结果
