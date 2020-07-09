# define the log-gamma distribution
dlgamma <- function(y, alpha, beta, log = FALSE){
  if (log == FALSE){
    out <- log(y + 1)^(alpha - 1)*(y + 1)^(-(1 + beta)/beta)/(beta^alpha)/gamma(alpha)
  } else if (log == TRUE) {
    logL <- -alpha*log(beta) - (1 + beta)/beta*log(1 + y) + (-1 + alpha)*log(log(1 + y)) - lgamma(alpha)
    out <- logL
  }
  return(out)
}

plgamma <- function(y, alpha, beta){
  z <- log(y + 1)
  out <- pgamma(z, shape = alpha, scale = beta)
  return(out)
}

rlgamma <- function(n, alpha, beta){
  xsim <- rgamma(n, shape = alpha, scale = beta)
  zsim <- exp(xsim) - 1
  return(zsim)
}

qlgamma <- function(p, alpha, beta){
  x <- qgamma(p = p, shape = alpha, scale = beta)
  z <- exp(x) - 1
  return(z)
}
