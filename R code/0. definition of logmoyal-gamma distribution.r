
# 1. pdf -------------------------------------------------------
dLMGA <- function(y, sigma, a, b, log = FALSE) {
  if( log == FALSE){
    exp(-0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b))
  } else if (log == TRUE){
    -0.5*log(2*pi) - log(sigma) + a*log(b) + lgamma(a + 0.5) - lgamma(a) - (1/(2*sigma)+1)*log(y) - (a + 0.5)*log(0.5*(1/y)^(1/sigma) + b)
  }
}
dLMGA(c(0.5,0.1), sigma = 2, a = 2, b = 3, log = FALSE)

# 2. cdf -------------------------------------------------------
pLMGA <- function(y, sigma, a, b) {
  # f <- function(y) dLMGA(y, sigma = sigma, a = a, b = b, log = FALSE)
  # FLMFA <- integrate(f, lower = 0, upper = y)$value
  # FLMFA
  z <- y^(-1/sigma)/(y^(-1/sigma) + 2*b)
  p <- 1- pbeta(z, shape1 = 0.5, shape2 = a)
  p
}
pLMGA(c(10,20), sigma = 2, a = 2, b = 3)

# 3. qf - quantile function -------------------------------------------------------
qLMGA <- function(u, sigma, a, b) {
  c <- (2*b)^(-sigma)
  #I <- pbeta(u, shape1 = 0.5, shape2 = a)
  Iinv <- qbeta(1- u, shape1 = 0.5, shape2 = a)
  c*(Iinv/(1 - Iinv))^(-sigma)
}
qLMGA(c(0.5,0.1), sigma = 2, a = 2, b = 3)

# 4. generating random-------------------------------------------------------
rLMGA <- function(n, sigma, a, b) {
  u <- runif(n, min = 0, max = 1)
  qLMGA <- Vectorize(qLMGA)
  r <- qLMGA(u, sigma = sigma, a = a, b = b)
  r
}
rLMGA(1000, sigma = 2, a = 2, b = 3)

# 5. mode, mean and variance -------------------------------------------------------
# mode function
modeLMGA <- function(sigma, a, b) {
  ((b + 2*b*sigma)/(a - sigma))^(-sigma)
}
# moment function
momLMGA <- function(h, sigma, a, b) {
  (2*b)^(-sigma*h)*beta(0.5 - h*sigma, a + h*sigma)/(beta(0.5, a))
}
# quantile function
varLMGA <- function(sigma, a, b) {
  momLMGA(h = 2, sigma = sigma, a = a, b = b) - (momLMGA(h = 1, sigma = sigma, a = a, b = b))^2
}
sigma <- 0.1; a <- 1; b <- 2

modeLMGA(sigma = sigma, a = a, b = b)
momLMGA(h = 1, sigma = sigma, a = a, b = b)
qLMGA(u = 0.5, sigma = sigma, a = a, b = b)

# 6. TVaR -------------------------------------------------------
TVaR.GLMGA <- function(p, sigma, a, b) {
  p1 <- (2*b)^(-sigma)*beta(0.5 - sigma, a + sigma)/beta(0.5, a)
  inv <- qbeta(1 - p, shape1 = 0.5, shape2 = a)
  p2 <- pbeta(inv, shape1 = 0.5 - sigma, shape2 = a + sigma)/(1-p)
  p1*p2
}
TVaR.GLMGA(p = 0.95, sigma = 0.1, a = 0.5, b = 0.1)

# 7. Excess-of-Loss reinsurance contract -------------------------------------------------------
EoL.GLMGA <- function(R, sigma, a, b){
  mean <- momLMGA(h = 1, sigma = sigma, a = a, b = b)
  z <- R^(-1/sigma)/(R^(-1/sigma) + 2*b)
  out <- (mean - R)*pbeta(z, shape1 = 0.5, shape2 = a)
  out
}
EoL.GLMGA(R = c(1, 0.4), sigma = 0.1, a = 0.5, b = 0.1)

# 8. Excess-of-Loss reinsurance contract -------------------------------------------------------
MeF.GLMGA <- function(u, sigma, a, b){
  mean <- momLMGA(h = 1, sigma = sigma, a = a, b = b)
  z <- u^(-1/sigma)/(u^(-1/sigma) + 2*b)
  I1 <- pbeta(z, shape1 = 0.5 - sigma, shape2 = a + sigma)
  I2 <- pbeta(z, shape1 = 0.5, shape2 = a)
  out <- (mean - u)*I1/I2
  out
}
MeF.GLMGA(u = c(0.4), sigma = 0.1, a = 0.5, b = 0.1)
