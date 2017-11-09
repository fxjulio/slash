#' Envelopes for C(alpha)
#' 
#' QQplot with envelopes
#' 
#' @export
envelopeForCalpha <- function( y, distr="mvnorm", method="moment", reps=20, conf=0.95, ... ){
stopifnot( is.matrix(y), distr %in% c("mvnorm", "slash", "mvt"), method %in% c("moment", "mve", "mcd") )
p <- ncol(y)
n <- nrow(y)

conf <- 1-conf

par <- estimator( y, distr, method)
d2 <- mahalanobis(y, par$mu, par$Sigma)
if( distr == "mvnorm"){

  u <- pchisq( d2, p )
  repl <- replicate(reps, rmvnorm(n, par$mu, par$Sigma)  )
  zrepl <- apply(repl, 3, function(x){
    par <- estimator(x, "mvnorm", method)
    d2 <- mahalanobis(x, par$mu, par$Sigma)
    u <- pchisq(d2, p)
    qnorm(sort(u))
  })
}
if( distr == "slash" ){

  u <- pdelta.slash(d2, par$eta, p)
  repl <- replicate(reps, sim.slash(n, par$mu, par$Sigma, par$eta)  )
  zrepl <- apply(repl, 3, function(x){
    par <- estimator(x, "slash", method)
    d2 <- mahalanobis(x, par$mu, par$Sigma)
    u <- pdelta.slash(d2, par$eta, p)
    qnorm(sort(u))
  })
}
if( distr == "mvt"){
  Ft <- 1/(1-2*par$eta)*d2/p

  u <- pf(Ft, p, 1/par$eta )
  repl <- replicate(reps,  rmt( n, par$mu, par$Sigma, par$eta) )
  zrepl <- apply(repl, 3, function(x){
    par <- estimator(x, "mvt", method)
    d2 <- mahalanobis(x, par$mu, par$Sigma)
    Ft <- 1/(1-2*par$eta)*d2/p
    u <- pf(Ft, p, 1/par$eta )
    qnorm(sort(u))
  })
}

bands <- apply(zrepl, 1, quantile, probs=c(conf/2, 1-conf/2), names=FALSE )

z <- qnorm(u)
qxy <- qqnorm( z , ...)
lines(sort(qxy$x), bands[1,])
lines(sort(qxy$x), bands[2,])


}

fn <- function(eta, y, mu, Sigma){
  p <- ncol(y)
  d2 <- mahalanobis(y,mu,Sigma)
  ceta <- eta/(1-2*eta)
  Kpeta <- (ceta/pi)^(p/2)*gamma((1+eta*p)/(2*eta))/gamma(1/(2*eta))
  detSigma <- det(Sigma)
  
  sum(log(Kpeta) -1/2*log(detSigma) -(1+eta*p)/(2*eta)*log(1+ceta*d2))
  
}

estimator <- function( y, distr="mvnorm", method="moment"){
  if( method == "moment"){
    est <- list( center=apply(y, 2, mean), cov=cov(y)  )
  } else {
    est <- cov.rob(y, method = method) # mve or mcd 
  }
  
  if( distr == "mvt" ){
    out <- optimize(fn, c(0.001, 0.4999), y=y, mu=est$center, Sigma=est$cov, maximum = TRUE)
    eta <- out$maximum
  }
  if( distr == "slash"){
    out <- startValues.slash( y, mu0=est$center, Sigma0 = est$cov )
    eta <- out$eta0ml
  }
  if( distr == "mvnorm" ){
    eta <- 0
  }
  
  list(mu=est$center, Sigma=est$cov, eta=eta )
  
}
