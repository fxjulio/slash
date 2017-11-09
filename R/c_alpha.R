#' @title Evaluate C(alpha).
#' 
#' @description
#' Evaluate C(alpha) for three cases according knows parameters.
#' 
#' @details
#' C(alpha) is based in sqrt(n)-consistent estimators. Only three cases are implemented.
#' Case 1: mu is know
#' Case 2: Sigma is know
#' Case 3: 
#' Case 4: sigma_11=...=sigma_pp and sigma_ij=sigma_lk for i neq j and l neq k
#' Case 5: mu_1=...=mu_p
#' Case 6: sigma_11=...=sigma_pp
#' Case 7: mu_1=...=mu_p and sigma_11=...=sigma_pp
#' 
#'
#' @param Y Matrix of slash samples (n rows and p columns).
#' @param test \code{c("case1", "case2", "equicorr", "mequ")}
#' @param mu0 H0 proposed mu parameter. Default \code{NULL}.
#' @param Sigma0 H0 proposed Sigma parameter. Default \code{NULL}.
#' @param eta0 H0 proposed eta parameter. Default \code{NULL}.
#'
#' @return
#' A \code{list} with \code{C.alpha}, \code{p.value}, \code{k1}, estimated values \code{sv}  values
#' and \code{test} applied.
#'
#' @examples
#' n <- 1000
#' p <- 2
#' eta <- 0.25
#' Y <- sim.slash( n, rep(0, p), matrix(c(2,0,0,3),2,2), eta )
#' 
#' Calpha.slash( Y, test="case1", mu0=c(0,0) )
#' 
#' Calpha.slash( Y, test="case2", Sigma0=diag(p) )
#' 
#'
#'
#' @export
Calpha.slash <- function( Y, test=c("case1", "case2", "equicorr", "mequ", "sequ", "msequ"),
    mu0=NULL, Sigma0=NULL, eta0=NULL ){
	p <- ncol(Y)
	n <- nrow(Y)
  
  k <- p + p*(p+1)/2 + 1
  
  if( test=="case1" & !is.null(mu0) ){
    k1 <- p
    C <- matrix(0, k1, k)
    C[,1:k1] <- diag(k1)
  }
  
  if( test=="case2" & !is.null(Sigma0) ){
    k1 <- p*(p+1)/2
    C <- matrix(0, k1, k)
    C[,p+(1:k1)] <- diag(k1)
  }
  
  if( test=="equicorr" ){
    A <- diag(p)
    vechA <- as.logical(vech(A))

    p.ast <- p*(p+1)/2
    k1 <- p.ast - 2

    dind <- 1
    tind <- p+1
    ind <- rep(0, p.ast)
    for( i in 1:p.ast){
      if( vechA[i] ){
        ind[i] <- dind
        dind <- dind + 1
      } else {
        ind[i] <- tind
        tind <- tind + 1
      }	
    }

    D <- matrix( 0, k1, p.ast)
    for( i in 1:(p-1) ){
      D[i,   i+(0:1)   ] <- c(1,-1)
    }

    if( p > 2 ){
      for( i in p:k1 ){
        D[i,   i+(1:2)   ] <- c(1,-1)
      }
    }

    C <- matrix(0, k1, k)
    C[,p+(1:p.ast)] <- D[,ind]

    S <- var(Y)
		sigma2hat <- mean(diag(S))
		rhohat <- mean( S[upper.tri(S)] ) / sigma2hat
		Sigma0 <- sigma2hat*((1-rhohat)*diag(p) + rhohat*matrix(1, p, 1)%*%matrix(1, 1, p))
   
  }
  
  if( test=="mequ" ){
    k1 <- p - 1 
    C <- matrix(0, k1, k)
    for( i in 1:k1){
      C[i, i+c(0,1)] <- c(1,-1)
    }
    
    mu0 <- mean( apply(Y, 2, mean) )*rep(1,p)
  }
  
  if( test=="sequ" ){
    k1 <- p - 1 
    A <- diag(p)
    vechA <- as.logical(vech(A))
    p.ast <- p * (p + 1)/2
    k1 <- p - 1
    ind <- which(vechA)

    C <- matrix(0, k1, k)
    for (i in 1:k1) {
       C[i, p + ind[i+ (0:1)] ] <- c(1, -1)
    }
    
    Sigma0 <- cov(Y)
    diag(Sigma0) <- mean(diag(Sigma0))

  }
  
  if( test=="msequ" ){
    k1 <- 2*p - 2 

    C <- matrix(0, k1, k)
    for( i in 1:(p-1)){
      C[i, i+c(0,1)] <- c(1,-1)
    }

    A <- diag(p)
    vechA <- as.logical(vech(A))
    p.ast <- p * (p + 1)/2
    ind <- which(vechA)

    for (i in 1:(p-1)) {
       C[i + p - 1, p + ind[i+ (0:1)] ] <- c(1, -1)
    }
    
    mu0 <- mean( apply(Y, 2, mean) )*rep(1,p)
    e <- t(apply( Y, 1, function(y) y-mu0))

    Sigma0 <- cov(e)
    diag(Sigma0) <- mean(diag(Sigma0))
    

  }
  
  
  
  sv <- startValues.slash(Y, mu0, Sigma0, eta0)
  
  I <- with(sv, infoFisher.slash(Y, mu0, Sigma0, eta0ml, block=FALSE))
  U <- with(sv, score.slash(Y, mu0, Sigma0, eta0ml, block=FALSE, hessian=FALSE))
  
  if( sv$eta0ml < 1e-3 | !is.null(eta0) ) {
    U <- U[-k,,drop=FALSE]
    I <- I[,-k]
    I <- I[-k,]
    C <- C[,-k,drop=FALSE]
  } 
  
  F <- solve(I)
  
  temp <- t(U)%*%F%*%t(C)
  
  C.alpha <- c( 1/n * temp%*%solve(C%*%F%*%t(C))%*%t(temp) )
  
  
  
  
  
  
  
  
  

	# if( equicorr & is.null(Sigma0) ){
    # S <- var(Y)
		# sigma2hat <- mean(diag(S))
		# rhohat <- mean( S[upper.tri(S)] ) / sigma2hat
		# Sigma0 <- sigma2hat*((1-rhohat)*diag(p) + rhohat*matrix(1, p, 1)%*%matrix(1, 1, p))
  # } else {
    # stop('equicorr=TRUE and Sigma0=NULL is required')
  # }
  
  # if( is.null(eta0) & eta.fixed ) stop("eta0 is NULL and eta.fixed is TRUE")
  
  # sv <- startValues.slash(Y, mu0, Sigma0, eta0)
  
	# if( eta.fixed ){
		# I <- with(sv, infoFisher.slash(Y, mu0, Sigma0, eta0, block=TRUE))
		# U <- with(sv, score.slash(Y, mu0, Sigma0, eta0, block=TRUE, hessian=FALSE))
	# } else {
		# I <- with(sv, infoFisher.slash(Y, mu0, Sigma0, eta0ml, block=TRUE))
		# U <- with(sv, score.slash(Y, mu0, Sigma0, eta0ml, block=TRUE, hessian=FALSE))
	# }


	# Calpha1 <- function(){
		# F.mu.mu <- solve(I$I.mu.mu)
		# U.mu <- as.matrix(U$U.mu)
		# list( C.alpha = c(1/n * t(U.mu)%*%F.mu.mu%*%U.mu), k1 = p)
	# }

	# Calpha2 <- function(){
		# c <- with(I,  c(I.eta.eta - t(I.phi.eta)%*%solve(I.phi.phi)%*%I.phi.eta) )
		# F.phi.eta <- with(I, -solve(I.phi.phi)%*%I.phi.eta / c )
		# F.phi.phi <- with(I, solve(I.phi.phi)%*%( diag(p*(p+1)/2) - I.phi.eta%*%t(F.phi.eta) ) )
		
		# K.phi.eta <- F.phi.phi %*% as.matrix(U$U.phi) + F.phi.eta * U$U.eta

		# list( C.alpha = c( 1/n * t(K.phi.eta) %*% solve(F.phi.phi) %*% K.phi.eta ), k1 = p*(p+1)/2 )
	# }
	
	# Calpha2etaFijo <- function(){
		# F.phi.phi <- solve(I$I.phi.phi)
		# U.phi <- as.matrix(U$U.phi)

		# list( C.alpha = c( 1/n * t(U.phi) %*% F.phi.phi %*% U.phi ), k1 = p*(p+1)/2 )
	# }
	


	# #case H1
	# if( !is.null(mu0) & is.null(Sigma0) ){
		# out <- Calpha1()
		# C.alpha <- out$C.alpha
		# k1 <- out$k1
	# }
	# #case H2 eta0 NO fijo
	# if( is.null(mu0) & !is.null(Sigma0) & !eta.fixed ){
    # out <- if( sv$eta0ml < 1e-3 ) Calpha2etaFijo() else Calpha2()
		# C.alpha <- out$C.alpha
		# k1 <- out$k1

	# }
	# #case H2 eta0 fijo
	# if( is.null(mu0) & !is.null(Sigma0) & eta.fixed ){
		# out <- Calpha2etaFijo()
		# C.alpha <- out$C.alpha
		# k1 <- out$k1

	# }
	
	# #case H3 eta NO fijo
	# if( !is.null(mu0) & !is.null(Sigma0) & !eta.fixed ){
		# out1 <- Calpha1()
		# out <- if( sv$eta0ml < 1e-3 ) Calpha2etaFijo() else Calpha2()
		# C.alpha <- out1$C.alpha + out2$C.alpha
		# k1 <- out1$k1 + out2$k1
	# }
	
	# #case H3 eta fijo
	# if( !is.null(mu0) & !is.null(Sigma0) & eta.fixed ){
		# out1 <- Calpha1()
		# out2 <- Calpha2etaFijo()
		# C.alpha <- out1$C.alpha + out2$C.alpha
		# k1 <- out1$k1 + out2$k1
	# }
  
  # if( equicorr ) k1 <- k1 - 2

	if( is.na(C.alpha) )
	return(list( C.alpha=NA, p.value=NA, k1=k1, sv=sv, test=test))
	
	if( C.alpha < 0 ) stop( "C.alpha negative" )
	p.value <- pchisq(C.alpha, k1, lower.tail=FALSE)

	list( C.alpha=C.alpha, p.value=p.value, k1=k1, sv=sv, test=test)

	

}
