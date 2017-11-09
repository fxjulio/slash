#' @title Evaluate C(alpha) for group comparation.
#' 
#' @description
#' Evaluate C(alpha) for two group comparation.
#' 
#' @details
#' C(alpha) is based in sqrt(n)-consistent estimators.
#'
#' @param g1 matrix (n1 rows and p columns).
#' @param g2 matrix (n2 rows and p columns).
#'
#' @return
#' A \code{list} with \code{C.alpha}, \code{p.value}, \code{k1} and estimated values \code{sv}.
#'
#' @examples
#' p <- 5
#' g1 <- sim.slash(100, rep(0,p), diag(p), 0.45)
#' g2 <- sim.slash(100, rep(1,p), diag(p), 0.45)
#' twoGroups.slash(g1, g2)
#'
#' g1 <- sim.slash(100, rep(0,p), diag(p), 0.45)
#' g2 <- sim.slash(100, rep(0,p), diag(p), 0.45)
#' twoGroups.slash(g1, g2)
#'
#' @export
twoGroups.slash <- function(g1, g2){
  Y <- rbind(g1, g2)
  p <- ncol(Y)
  n <- nrow(Y)

  sv <- startValues.slash(Y)

  ###########
  # Las Informaciones
  ########


  I1 <- with(sv, infoFisher.slash( g1, mu0, Sigma0, eta0ml, block=TRUE ))$I.mu.mu
  I2 <- with(sv, infoFisher.slash( g2, mu0, Sigma0, eta0ml, block=TRUE ))$I.mu.mu


  U1 <- with(sv, score.slash(g1, mu0, Sigma0, eta0ml, block=TRUE) )$U.mu
  U2 <- with(sv, score.slash(g2, mu0, Sigma0, eta0ml, block=TRUE) )$U.mu

  

  I <- matrix(0, 2*p, 2*p)
  I[1:p,1:p]<-I1
  I[p+(1:p),p+(1:p)]<-I2

  U <- c(U1,U2)

  k1 <- p 
  C <- matrix(0, k1, 2*p)
  for( i in 1:k1){
    C[i,i+c(0,p)]<-c(1,-1)
  }

  I.inv <- solve(I)
  temp <- t(U)%*%I.inv%*%t(C)

  C.alpha <- c(temp%*%solve(C%*%I.inv%*%t(C))%*%t(temp)/n)

  p.value <- pchisq(C.alpha, k1,  lower=FALSE) #p-value
  
  list( C.alpha=C.alpha, p.value=p.value, k1=k1, sv=sv)

}

