# ======================================================================================
# fire claim data
# ======================================================================================

library(SMPracticals)
library(data.table)
library(ggplot2)
library(gamlss)
setwd("D:/对外经济贸易大学 - 科研/0. Generalizing the log-Moyal distribution and regression models for heavy tailed loss data/合作文件 - 杨亮 - 文稿+代码")
source('0. definition of logmoyal distribution.r')
source('0. definition of logmoyal-gamma distribution.r')
# dt <- fread('fire loss experience of a major university.csv')
# y <- dt$`Total losses`
# loss data set - 2
y <- c(290.40,1248.49,2202.96,3941.30,10560.10,
       537.19,1268.24,2222.80,4017.01,11179.54,
       756.80,1284.56,2255.72,4100.00,11461.39,
       769.19,1363.85,2274.61,4166.98,14538.13,
       787.69,1436.20,2328.64,4355.02,14789.81,
       796.18,1445.96,2384.37,5117.93,17186.09,
       933.62,1469.48,2847.83,5335.96,18582.57,
       967.97,1507.47,2947.04,5453.02,22857.33,
       1010.56,1662.36,2948.35,5568.96,23177.85,
       1017.40,1674.58,3036.51,5761.83,23446.13,
       1033.49,1690.91,3287.68,6161.81,28409.82,
       1034.33,1739.96,3331.62,6348.69,57612.82,
       1056.93,1776.56,3416.67,6859.37,59582.78,
       1124.09,1932.09,3604.66,7972.20,113164.70,
       1165.73,1975.89,3671.16,8028.32,123228.90,
       1217.64,2099.79,3739.30,10047.22,626402.80)
# ---------------------------------------------------------
# difine the hessian matrix output and the stardard error -------------------
# ---------------------------------------------------------
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
# ================================================================
# 1. log moyal distribution
# =================================================================
LLlogmoyal <- function(y, pars) {
  sigma <- exp(pars[1])
  mu <- exp(pars[2])
  log_lik <- 0
  for (i in 1:length(y)) {
    log_lik[i] = -0.5*log(2*3.1415926) - log(sigma) - log(y[i]) + 1/(2*sigma)*(log(mu) - log(y[i])) - 0.5*(mu/y[i])^(1/sigma);
  } 
  loglike <- -sum(log_lik)  
  return(loglike)
}
# mlogmoyal <- optim(par = c(1,-1),
#                       fn = LLlogmoyal,
#                       y = y)
mlogmoyal <- optim(par = c(-1.1345, 0.27),
                   fn = LLlogmoyal,
                   y = y, hessian = T)
mlogmoyal
modout(mlogmoyal)
pars <- mlogmoyal$par
sigma <- exp(pars[1])
mu <- exp(pars[2])
tau <- mu^(1/sigma)

ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pGlogM(y[i], sigma = sigma, tau = tau)
}
resGlogM <- qnorm(ufit)
Var95.GlogM <- qGlogM(0.95, sigma = as.numeric(unique(sigma)), tau = tau)
Var99.GlogM <- qGlogM(0.99, sigma = as.numeric(unique(sigma)), tau = tau)



# ================================================================
# 2. log moyal - gamma distribution
# ================================================================
LLlogmoyalGA <- function(y, pars) {
  sigma <- exp(pars[1])
  b <- exp(pars[2])
  kappa <- exp(pars[3])
  #ll <- -0.5*log(2) - log(sigma) + 0.5*log(b) - lbeta(0.5, 0.5) - (1/(2*sigma)+1)*log(y) - (0.5 + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  ll <- -0.5*log(2*pi) - log(sigma) + kappa*log(b) + lgamma(kappa + 0.5) - lgamma(kappa) - (1/(2*sigma)+1)*log(y) - (kappa + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  
  loglike <- -sum(ll)  
  return(loglike)
}
pars <- c(1,1)
mlogmoyalGA <- optim(par = c(1,-1, -1),
                      fn = LLlogmoyalGA,
                      y = y, hessian = T)

mlogmoyalGA
modout(mlogmoyalGA)

pars <- mlogmoyalGA$par
sigma <- exp(pars[1])
b <- exp(pars[2])
a <- exp(pars[3])
ufit <- c()
for(i in 1:length(y)){
  ufit[i] <- pLMGA(y[i], sigma = sigma, a = a, b = b)
}
resLMGA <- qnorm(ufit)
Var95.GLMGA <- qLMGA(0.95, sigma = sigma, a = a, b = b)
Var99.GLMGA <- qLMGA(0.99, sigma = sigma, a = a, b = b)

# TVar95.GLMGA <- TVaR.GLMGA(0.95, sigma = sigma, a = a, b = b)
# TVar99.GLMGA <- TVaR.GLMGA(0.99, sigma = sigma, a = a, b = b)
I.GLMGA <- as.numeric(y < Var95.GLMGA)
length(I.GLMGA)
length(I.GLMGA[I.GLMGA==1])
pbinom(length(I.GLMGA[I.GLMGA==1]), size = length(I.GLMGA), prob = length(I.GLMGA[I.GLMGA==1])/length(I.GLMGA))
# -------------------------------------------------------------
# 3. GB2 distribution
# -------------------------------------------------------------
LLGB2 <- function(pars, y){
  sigma <- (pars[2])  # shape
  mu <- exp(pars[1]) # scale
  nu <- exp(pars[3])
  tau <- exp(pars[4])
  logL <- dGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau, log = T)
  return(-sum(logL))
}
mlGB2 <- optim(fn = LLGB2, y = y, 
               #par = c(1.5, 1.5, 1, 0.05), 
               par = c(0.5,2, 1, 1),
               hessian = T, 
               control = list(maxit = 50000))
mlGB2
modout(mlGB2)

pars <- mlGB2$par
sigma <- (pars[2])  # shape
mu <- exp(pars[1]) # scale
nu <- exp(pars[3])
tau <- exp(pars[4])
ufit <- pGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
resGB2 <- qnorm(ufit)
Var95.GB2 <- qGB2(0.95, mu = mu, sigma = sigma, nu = nu, tau = tau)
Var99.GB2 <- qGB2(0.99, mu = mu, sigma = sigma, nu = nu, tau = tau)


# -------------------------------------------------------------
# 4. inverse Lindely distribution
# -------------------------------------------------------------
LLGIL <- function(pars, y) {
  lambda <- exp(pars[1])
  alpha <- exp(pars[2])
  
  ll <- log(alpha) + 2*log(lambda) - log(1 + lambda) + log(1 + y^(-alpha)) - (alpha + 1)*log(y) - lambda*y^(-alpha)
  return(-sum(ll))
}
mGIL <- optim(fn = LLGIL, y = y, 
               par = c(-1,1),
               hessian = T, 
               control = list(maxit = 50000))

mGIL
modout(mGIL)

pars <- mGIL$par
lambda <- exp(pars[1])
alpha <- exp(pars[2])
pGIL <- function(y, lambda, alpha){
  (1 + lambda + lambda*y^(-alpha))/(1 + lambda)*exp(-lambda*y^(-alpha))
}
ufit <- pGIL(y, lambda = lambda, alpha = alpha)
resGIL <- qnorm(ufit)
Var95.GIL <- NA
Var99.GIL <- NA
# -------------------------------------------------------------
# 5. Lomax distribution
# -------------------------------------------------------------
X <- model.matrix(~ 1, data = data.table(y))
LLpareto <- function(pars, y, X){
  lengthpar <- dim(X)[2]
  beta <- exp(pars[lengthpar+1])  # 同样的
  sigma <- exp(X%*%pars[1:lengthpar]) # 同样的
  logL <- log(beta) + beta*log(sigma) - (beta+1)*log(y+sigma) # 对数似然函数
  ll <- -sum(logL)
  return(ll)
}
mlpareto <- optim(fn = LLpareto, 
                  y = y, X = X, 
                  par = rep(times = (dim(X)[2]+1), x = 0.1), hessian = T)
mlpareto
modout(mlpareto)
pars <- mlpareto$par
lengthpar <- dim(X)[2]
beta <- exp(pars[lengthpar+1])  # 同样的
sigma <- exp(X%*%pars[1:lengthpar]) # 同样的
ufit <- ppareto2(y, shape = beta, scale = sigma)
respareto <- qnorm(ufit)
Var95.pareto <- qpareto2(0.95, shape = beta, scale = unique(sigma))
Var99.pareto <- qpareto2(0.99, shape = beta, scale = unique(sigma))
# -------------------------------------------------------------
# 6. log gamma distribution
# -------------------------------------------------------------
library(actuar)
pars <- c(0, 0)
LLlogG <- function(pars, y){
  alpha <- exp(pars[1])  # 同样的
  beta <- exp(pars[2]) # 同样的
  #logL <- (alpha - 1)*log(log(y)) - log(y) - alpha*log(beta) - lgamma(alpha) - log(y)/beta
  logL <- dlgamma(y, shapelog = alpha, ratelog = beta, log = T)
  ll <- -sum(logL)
  return(ll)
}
mllogG <- optim(fn = LLlogG, 
                  y = y, 
                  par = c(0, 0), hessian = T)
mllogG
modout(mllogG)
pars <- mllogG$par
alpha <- exp(pars[1])  # 同样的
beta <- exp(pars[2]) # 同样的
ufit <- plgamma(y, shapelog = alpha, ratelog = beta)
reslogG <- qnorm(ufit)
Var95.logG <- qlgamma(0.95, shapelog = alpha, ratelog = beta)
Var99.logG <- qlgamma(0.99, shapelog = alpha, ratelog = beta)
# -------------------------------------------------------------
# 7. Frechet distribution
# -------------------------------------------------------------
LLfrechet <- function(pars, y){
  a <- exp(pars[1])  # 同样的
  b <- exp(pars[2]) # 同样的
  #logL <- (alpha - 1)*log(log(y)) - log(y) - alpha*log(beta) - lgamma(alpha) - log(y)/beta
  logL <- dweibull(y, shape = a, scale = b, log = T)
  ll <- -sum(logL)
  return(ll)
}
mfrechet <- optim(fn = LLfrechet, 
                y = y, 
                par = c(0, 0), hessian = T)
mfrechet
modout(mfrechet)
pars <- mfrechet$par
a <- exp(pars[1])  # 同样的
b <- exp(pars[2]) # 同样的
ufit <- pweibull(y, shape = a, scale = b)
resfrechet <- qnorm(ufit)
Var95.frechet <- qweibull(0.95, shape = a, scale = b)
Var99.frechet <- qweibull(0.99, shape = a, scale = b)



# qq-plot ------------------------------------------------------------------
par(mfrow = c(2, 3))
library(extRemes)
q.GlogM <- qqnorm(resGlogM, xlim = c(-4,4), ylim = c(-4,4), main = 'GlogM')
abline(0, 1, col = 'blue')
q.GLMGA <- qqnorm(resLMGA, xlim = c(-4,4), ylim = c(-4,4), main = 'GLMGA')
abline(0, 1, col = 'blue')
q.GB2 <- qqnorm(resGB2, xlim = c(-3,3), ylim = c(-3,3), main = 'GB2')
abline(0, 1, col = 'blue')
q.Lomax <- qqnorm(respareto, xlim = c(-3,3), ylim = c(-3,3), main = 'Lomax')
abline(0, 1, col = 'blue')
q.logGA <- qqnorm(reslogG, xlim = c(-3,3), ylim = c(-3,3), main = 'Log-gamma')
abline(0, 1, col = 'blue')
q.Frechet <- qqnorm(resfrechet, xlim = c(-3,3), ylim = c(-3,3), main = 'Frechet')
abline(0, 1, col = 'blue')



# ------------------------------------------------------------------
# ks test + ad test + cvm test
# ------------------------------------------------------------------
library(goftest)
pars <- mlogmoyal$par
sigma <- exp(pars[1])
mu <- exp(pars[2])
tau <- mu^(1/sigma)
pGlogMnew <- function(y) pGlogM(y, sigma = sigma, tau = tau)
ksGlogM <- ks.test(x = y, y = pGlogMnew)
adGlogM <- ad.test(x = y, 'pGlogM', sigma = sigma, tau = tau)
cvmGlogM <- cvm.test(x = y, 'pGlogM', sigma = sigma, tau = tau)

pars <- mlogmoyalGA$par
sigma <- exp(pars[1])
b <- exp(pars[2])
a <- exp(pars[3])
pLMGAnew <- function(y) pLMGA(y, sigma = sigma, a = a, b = b)
ksGLMGA <- ks.test(x = y, y = pLMGAnew)
adGLMGA <- ad.test(x = y, 'pLMGA', sigma = sigma, a = a, b = b)
cvmGLMGA <- cvm.test(x = y, 'pLMGA', sigma = sigma, a = a, b = b)


pars <- mlGB2$par
sigma <- (pars[2])
mu <- exp(pars[1])
nu <- exp(pars[3])
tau <- exp(pars[4])
pGB2new <- function(y) pGB2(y, mu = mu, sigma = sigma, nu = nu, tau = tau)
ksGB2 <- ks.test(x = y, y = pGB2new)
adGB2 <- ad.test(x = y, 'pGB2new')
cvmGB2 <- cvm.test(x = y, 'pGB2new')

pars <- mGIL$par
lambda <- exp(pars[1])
alpha <- exp(pars[2])
pGILnew <- function(y) pGIL(y, lambda = lambda, alpha = alpha)
ksGIL <- ks.test(x = y, y = pGILnew)
adGIL <- ad.test(x = y, 'pGILnew')
cvmGIL <- cvm.test(x = y, 'pGILnew')

pars <- mlpareto$par
lengthpar <- dim(X)[2]
beta <- exp(pars[lengthpar+1])  # 同样的
sigma <- exp(X%*%pars[1:lengthpar]) # 同样的
plomaxnew <- function(y) ppareto2(y, shape = beta, scale = sigma)
kslomax <- ks.test(x = y, y = plomaxnew)
adlomax <- ad.test(x = y, 'plomaxnew')
cvmlomax <- cvm.test(x = y, 'plomaxnew')

pars <- mllogG$par
alpha <- exp(pars[1])  # 同样的
beta <- exp(pars[2]) # 同样的
plogMnew <- function(y) plgamma(y, shapelog = alpha, ratelog = beta)
kslogM <- ks.test(x = y, y = plogMnew)
adlogM <- ad.test(x = y, 'plogMnew')
cvmlogM <- cvm.test(x = y, 'plogMnew')

pars <- mfrechet$par
a <- exp(pars[1])  # 同样的
b <- exp(pars[2]) # 同样的
ufit <- pweibull(y, shape = a, scale = b)
#resfrechet <- qnorm(ufit)
pweibullnew <- function(y) pweibull(y, shape = a, scale = b)
ksfre <- ks.test(x = y, y = pweibullnew)
adfre <- ad.test(x = y, 'pweibullnew')
cvmfre <- cvm.test(x = y, 'pweibullnew')

# ---------------------------------------------------------
pars.com <- rbind(modout(mlogmoyal)$summary,
                  modout(mlGB2)$summary, 
                  modout(mlogmoyalGA)$summary,
                  modout(mGIL)$summary,
                  modout(mlpareto)$summary,
                  modout(mllogG)$summary,
                  modout(mfrechet)$summary)
exp(pars.com)
npars.com <- rbind(nrow(modout(mlogmoyal)$summary),
                        nrow(modout(mlGB2)$summary), 
                             nrow(modout(mlogmoyalGA)$summary),
                                  nrow(modout(mGIL)$summary),
                                       nrow(modout(mlpareto)$summary),
                                            nrow(modout(mllogG)$summary),
                                                 nrow(modout(mfrechet)$summary))

LL.com <- rbind(modout(mlogmoyal)$ll,
                  modout(mlGB2)$ll, 
                  modout(mlogmoyalGA)$ll,
                  modout(mGIL)$ll,
                  modout(mlpareto)$ll,
                  modout(mllogG)$ll,
                  modout(mfrechet)$ll)

AIC.com <- rbind(modout(mlogmoyal)$AIC,
                 modout(mlGB2)$AIC, 
                 modout(mlogmoyalGA)$AIC,
                 modout(mGIL)$AIC,
                 modout(mlpareto)$AIC,
                 modout(mllogG)$AIC,
                 modout(mfrechet)$AIC)

BIC.com <- rbind(modout(mlogmoyal)$BIC, modout(mlGB2)$BIC, 
                 modout(mlogmoyalGA)$BIC,
                 modout(mGIL)$BIC,
                 modout(mlpareto)$BIC,
                 modout(mllogG)$BIC,modout(mfrechet)$BIC)
ks.com <- rbind(ksGlogM$p.value, 
                ksGB2$p.value, ksGLMGA$p.value, ksGIL$p.value, 
                kslomax$p.value, kslogM$p.value,
                ksfre$p.value)
ad.com <-  rbind(adGlogM$p.value, 
                 adGB2$p.value, adGLMGA$p.value, adGIL$p.value, 
                 adlomax$p.value, adlogM$p.value,
                 adfre$p.value)
cvm.com <-  rbind(cvmGlogM$p.value, 
                 cvmGB2$p.value, cvmGLMGA$p.value, cvmGIL$p.value, 
                 cvmlomax$p.value, cvmlogM$p.value,
                 cvmfre$p.value)
var95.com <- rbind(Var95.GlogM, Var95.GB2,
                   Var95.GLMGA,
                   Var95.GIL, Var95.pareto, Var95.logG,  
                   Var95.frechet)
var99.com <- rbind(Var99.GlogM, Var99.GB2,
                   Var99.GLMGA,
                   Var99.GIL, Var99.pareto, Var99.logG,  
                   Var99.frechet)
cor.com <- rbind(cor(q.GlogM$data, q.GlogM$qnorm),
                 cor(q2$data, q2$qnorm),
                 cor(q.GLMGA$data, q.GLMGA$qnorm),
                 NA,
                 cor(q.Lomax$data, q.Lomax$qnorm),
                 cor(q.logGA$data, q.logGA$qnorm),
                 cor(q.Frechet$data, q.Frechet$qnorm))
var95.dif <- (var95.com - quantile(y, 0.95))/quantile(y, 0.95)
var99.dif <- (var99.com - quantile(y, 0.99))/quantile(y, 0.99)
row.names(AIC.com) <- c('logmoyal', 'GB2', 
                        'logmoyalGA', 
                        'Inverse-Lind', 'Lomax', 'log-gamma',
                        'frechet')
colnames(npars.com) <- 'Npars'
colnames(LL.com) <- 'LL'
colnames(AIC.com) <- 'AIC'
colnames(BIC.com) <- 'BIC'
colnames(ks.com) <- 'KS'
colnames(ad.com) <- 'AD'
colnames(cvm.com) <- 'CVM'
colnames(var95.com) <- 'Var95'
colnames(var99.com) <- 'Var99'
colnames(var95.dif) <- 'var95.dif'
colnames(var99.dif) <- 'var99.dif'
colnames(cor.com) <- 'cor.residual'
results.com <- cbind(npars.com, 
                     LL.com, round(cbind(AIC.com, BIC.com), 1), 
                     ks.com, ad.com, cvm.com,
                     var95.com,
                     Evar95 = quantile(y, 0.95),
                     Evar99 = quantile(y, 0.99),
                     var99.com,
                     var95.dif,
                     var99.dif)
results.com
# results.com <- data.table(results.com)
# results.com[, var95.dif := (Var95 - Evar95)/Evar95]
# results.com[, var99.dif := (Var99 - Evar99)/Evar99]
# row.names(results.com) <- c('logmoyal', 'GB2', 
#                         'logmoyalGA', 
#                         'Inverse-Lind', 'Lomax', 'log-gamma',
#                         'frechet')
# results.com[order(abs(var99.dif))]
results.com

#write.csv(results.com, 'Tables/var - data (1).csv')
# write.csv(pars.com, 'Tables/parsEST - data (1).csv')





# ------------------------------------------------------------------
# VaR risk measures
# ------------------------------------------------------------------

