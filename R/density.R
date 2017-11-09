#
dens.slash <- function( Y, mu=rep(0,2), Sigma=diag(2), eta=0.25){
  #verificar p
  Y <- as.matrix(Y)
  p <- ncol(Y)
  ceta <- 1/(1-eta)
  K <- (ceta)^(p/2)/eta/(2*pi)^(p/2)
  if( p > 1){
    Sigma.inv <- solve(Sigma)
    detSigma <- det(Sigma)
  } else {
    Sigma.inv <- 1/Sigma
    detSigma <- Sigma
  }
  a <- p/2 + 1/eta
  
  G <- function(a,b) integrate( function(v) v^(a-1)*exp(-v*b) , 0, 1)$value
 
  apply( Y, 1, function(y){
    delta <- t(y-mu)%*%Sigma.inv%*%(y-mu)
    b <- ceta*delta/2
    Gab <- G(a, b)
    K*detSigma^(-1/2)*Gab  
  })

}