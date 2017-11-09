#' Asymptotic covariance Omega
#'
#' Evaluate asymptotic covariance of 
#'
#' @param Sigma Covariance matrix estimator.
#' @param eta Shape parameter estimator.
#' @param kappa Multivariate excess of kurtosis estimator.
#'
#' @return
#' A square matrix of p + p*(p+1)/2
#'
#' @export
OmegaMat <- function( Sigma, eta, kappa=NULL ){
  p <- ncol(Sigma)
  past <- p*(p+1)/2
  
  out <- matrix(0, p+past, p+past)
  
  out[1:p, 1:p] <- Sigma
  
  kappa <- ifelse( is.null(kappa), eta^2/(1-2*eta), kappa )
   
  Dp <- D.matrix( p )
  H <- solve(t(Dp)%*%Dp, t(Dp) )
  Kp <- K.matrix( p )
  
  out[p+(1:past), p+(1:past)] <- H%*%( (kappa + 1)*(diag(p^2) + Kp)%*%(Sigma%x%Sigma) + kappa*vec(Sigma)%*%t(vec(Sigma))  )%*%t(H)
  
  out
}
