# GLMGA-model
Many catastrophic loss data are known to be heavy-tailed with extreme value.
In this article, a new family of probability distributions with domain in $\mathbb{R}^+$ is
introduced.
This class can be considered as a natural extension of the generalized logMoyal distribution.
This new family is obtained through the mixture of generalized logMoyal distribution with gamma distribution, suitable for modelling heavy tailed data. We also show some important features such as expressions of probability density function, moments, etc. Statistical inference of the model parameters is discussed using the method of maximum likelihood estimation. Parametric regression modelling is also discussed assuming that the response variable follows the generalized logMoyal-gamma distribution.
A Chinese earthquake loss data-set is used to show the applicability of the new class of distributions.

```r

# 1. pdf
dGlogM <- function(y, sigma, tau, log = FALSE) {
  if( log == FALSE){
    exp(0.5*log(tau) - 0.5*log(2*pi) - log(sigma) - (1/(2*sigma) + 1)*log(y) - tau/2*(1/y)^(1/sigma))
  } else if (log == TRUE){
    0.5*log(tau) - 0.5*log(2*pi) - log(sigma) - (1/(2*sigma) + 1)*log(y) - tau/2*(1/y)^(1/sigma)
  }
}
dGlogM(c(0.5,0.1), sigma = 2, tau = 1, log = FALSE)

# 2. cdf
pGlogM <- function(y, sigma, tau) {
  mu <- tau^sigma
  z <- 1/sqrt(2)*(mu/y)^(1/(2*sigma))
  p <- pracma::erfc(z)
  p
}
pGlogM(c(0.5,0.1), sigma = 2, tau = 1)

# 3. qf - quantile function
qGlogM <- function(u, sigma, tau) {
  mu <- tau^sigma
  mu*(sqrt(2)*pracma::erfcinv(u))^(-2*sigma)
}
qGlogM(c(0.5,0.1), sigma = 2, tau = 1)


# 4. generating random 
rGlogM <- function(n, sigma, tau) {
  u <- runif(n, min = 0, max = 1)
  r <- qGlogM(u, sigma = sigma, tau = tau)
  r
}
rGlogM(1000, sigma = 2, tau = 1)







```
