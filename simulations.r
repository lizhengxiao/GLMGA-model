
source('0. definition of logmoyal distribution.r')
source('0. definition of logmoyal-gamma distribution.r')

# 1. define the loglikelihood function of GLMGA -------------------------------------------------
LLlogmoyalGA <- function(y, pars, Xsigma, Xb) {
  sigma <- exp(Xsigma %*% pars[1:dim(Xsigma)[2]]) # introduce the covariates in parameter sigma using log link function 
  b <- exp(Xb %*%pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])]) # introduce the covariates in parameter using log link function 
  a <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
  ll <- -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  loglike <- -sum(ll)  
  return(loglike)
}

# 2: n = 200 and a = 0.5 ---------------------------------------
Nsim <- 2000 # simulation times
beta.est <- alpha.est <- matrix(0, nrow = Nsim, ncol = 2)
a.est <- sigma.est <- b.est <- matrix(0, nrow = Nsim, ncol = 1)
convergence.vector <- rep(0, length.out = Nsim)
asy.var <- matrix(0, nrow = Nsim, ncol = 5)
set.seed(112)
for(i in 1:Nsim){
  repeat{
    n <- 200      # sample size 
    beta <- c(-1, 0.5)
    alpha <- c(1, 0.5)
    x1 <- rnorm(n, 0, 1)
    x2 <- rnorm(n, 0, 1)
    Xsigma <- model.matrix(~x1)  # density matrix 
    Xb <- model.matrix(~x2)      # density matrix 
    sigma <- exp(Xsigma%*%beta)
    b <- exp(Xb%*%alpha)
    a <- 0.5  # true value of parameter a
    y <- rLMGA(n = n, sigma = sigma, a = a, b = b)

    mlogmoyalGA <- optim(fn = LLlogmoyalGA, Xsigma = Xsigma,
                          Xb = Xb,
                          y = y, hessian = T,
                          control = list(maxit = 50000),
                          method = 'Nelder-Mead',
                          par = c(beta, alpha, log(a)))
    if (mlogmoyalGA$convergence == 0) break
  }
  convergence.vector[i] <- mlogmoyalGA$convergence
  pars <- mlogmoyalGA$par
  asy.var[i,] <- diag(solve(mlogmoyalGA$hessian))*c(1,1,1,1,exp(pars[5])^2)
  beta.est[i,] <- pars[1:dim(Xsigma)[2]]
  alpha.est[i,] <- pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])]
  a.est[i,] <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
}

truevalue <- c(beta, alpha, (a))
# sample mean
apply(beta.est, 2, FUN = mean)
apply(alpha.est, 2, FUN = mean)
mean(a.est)
smean <- c(apply(beta.est, 2, FUN = mean),
  apply(alpha.est, 2, FUN = mean), mean(a.est))
# sample variance
apply(beta.est, 2, FUN = var)
apply(alpha.est, 2, FUN = var)
var(a.est)
svar <- c(apply(beta.est, 2, FUN = var), apply(alpha.est, 2, FUN = var), var(a.est))
# asymptotic variance
asvar <- apply(asy.var, 2, FUN = mean)
# ratio, bias, relative bias, mse
ratio <- svar/asvar
bias <- smean - truevalue
relbias <- smean/truevalue
beta.mat <- matrix(beta, nrow = nrow(beta.est), ncol = ncol(beta.est), byrow = T)
alpha.mat <- matrix(alpha, nrow = nrow(alpha.est), ncol = ncol(alpha.est), byrow = T)
a.mat <- matrix(a, nrow = nrow(a.est), ncol = ncol(a.est), byrow = T)
mse <- (c(apply((beta.est - beta.mat)^2, MARGIN = 2, FUN = mean),
apply((alpha.est - alpha.mat)^2/Nsim, MARGIN = 2, FUN = sum),
sum((a.est-(a.mat))^2)/Nsim))
mse


# 3. n = 100 to 1000 (a = 1)-------------------------------------------------
Nsim <- 2000 
beta.est <- alpha.est <- matrix(0, nrow = Nsim, ncol = 2)
a.est <- sigma.est <- b.est <- matrix(0, nrow = Nsim, ncol = 1)
convergence.vector <- c()
asy.var <- matrix(0, nrow = Nsim, ncol = 5)
#n.vector <- seq(from = 100, to = 1000, by = 20)
n.vector <- seq(from = 100, to = 1000, by = 500)
a.true.vector <- 1
smean <- svar <- asvar <- ratio <- bias <- relbias <- mse <- matrix(0, nrow = length(n.vector), ncol = 5)

set.seed(112)
for(j in 1:length(n.vector)){
  n <- n.vector[j]
  for(i in 1:Nsim){
    repeat{
      a <- a.true.vector
      beta <- c(-1, 0.5)
      alpha <- c(1, 0.5)
      x1 <- rnorm(n, 0, 1)
      x2 <- rnorm(n, 0, 1)
      Xsigma <- model.matrix(~ x1)
      Xb <- model.matrix(~ x2) 
      sigma <- exp(Xsigma%*%beta)
      b <- exp(Xb%*%alpha)
      y <- rLMGA(n = n, sigma = sigma, a = a, b = b)
      y[y == 0] <- 0.00001
      mlogmoyalGA <- optim(fn = LLlogmoyalGA, Xsigma = Xsigma,
                            Xb = Xb,
                            y = y, hessian = T,
                            control = list(maxit = 50000),
                            method = 'Nelder-Mead',
                            par = c(beta, alpha, log(a)))
      convergence.vector[i] <- mlogmoyalGA$convergence
      if (convergence.vector[i] == 0) break
    }
    pars <- mlogmoyalGA$par
    asy.var[i,] <- diag(solve(mlogmoyalGA$hessian))*c(1,1,1,1,exp(pars[5])^2)
    beta.est[i,] <- pars[1:dim(Xsigma)[2]]
    alpha.est[i,] <- pars[(dim(Xsigma)[2]+1):(dim(Xsigma)[2] + dim(Xb)[2])]
    a.est[i,] <- exp(pars[dim(Xsigma)[2] + dim(Xb)[2]+1])
  }
  truevalue <- c(beta, alpha, (a))
  # sample mean
  smean[j,] <- c(apply(beta.est, 2, FUN = mean),
                 apply(alpha.est, 2, FUN = mean), mean(a.est))
  # sample variance
  svar[j,] <- c(apply(beta.est, 2, FUN = var), apply(alpha.est, 2, FUN = var), var(a.est))
  # asymptotic variance
  asvar[j,] <- apply(asy.var, 2, FUN = mean)
  # ratio, bias, relative bias
  ratio[j,] <- svar[j,]/asvar[j,]
  bias[j,] <- smean[j,] - truevalue
  relbias[j,] <- smean[j,]/truevalue
  # mse
  beta.mat <- matrix(beta, nrow = nrow(beta.est), ncol = ncol(beta.est), byrow = T)
  alpha.mat <- matrix(alpha, nrow = nrow(alpha.est), ncol = ncol(alpha.est), byrow = T)
  a.mat <- matrix(a, nrow = nrow(a.est), ncol = ncol(a.est), byrow = T)
  mse[j,] <- (c(apply((beta.est - beta.mat)^2, MARGIN = 2, FUN = mean),
                apply((alpha.est - alpha.mat)^2/Nsim, MARGIN = 2, FUN = sum),
                sum((a.est-(a.mat))^2)/Nsim))
}

rbind(smean, svar, asvar, ratio, relbias, mse, bias)

