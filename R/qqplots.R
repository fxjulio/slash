#' @title qqplot multivariate norm
#' @export
qqplotmvnorm <- function(dd){
  dd <- as.matrix(dd)
  p <- ncol(dd)
  
  mu <- apply(dd, 2, mean)
  Sigma <- cov(dd)
  
  delta <- apply( dd, 1, function(x)  t(x-mu)%*%solve(Sigma)%*%(x-mu) )
  
  u <- pchisq(delta, p)
  z <- qnorm(u)
  
  qqnorm(z, main="Normal")
  qqline(z)

}

#' @title qqplot multivariate slash
#' @export
qqplotmvslash <- function(dd, fit) {
  p <- ncol(dd)
  delta <- with(fit, apply( dd, 1, function(x)  t(x-mu)%*%solve(Sigma)%*%(x-mu) ))
  
  u <- pdelta.slash(delta, fit$eta, p)
  z <- qnorm(u)
  
  qqnorm(z, main="Slash")
  qqline(z)
}