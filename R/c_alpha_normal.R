#' @title C(alpha) under multivariate normal distribution
#' @export

Calpha.normal <- function(Y, test=c("case1","case2","equicorr"),
    mu0=NULL, Sigma0=NULL){
	n <- nrow(Y)
	p <- ncol(Y)
	k <- p + p*(p+1)/2
  
  if( test == "equicorr" ){
    S <- var(Y)
		sigma2hat <- mean(diag(S))
		rhohat <- mean( S[upper.tri(S)] ) / sigma2hat
		Sigma0 <- sigma2hat*((1-rhohat)*diag(p) + rhohat*matrix(1, p, 1)%*%matrix(1, 1, p))
    
    mu0 <- apply(Y, 2, mean)
    
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

    if( p > 2){
      for( i in p:k1 ){
        D[i,   i+(1:2)   ] <- c(1,-1)
      }
    }

    C <- matrix(0, k1, k)
    C[,p+(1:p.ast)] <- D[,ind]
    
  }
  

	if( test == "case1" & !is.null(mu0) ){

		k1 <- p
		C <- matrix(0, k1, k)
		C[,1:p]<-diag(k1)

    e <- t(apply(Y,1,function(x) x-mu0 ))
		Sigma0<-cov(e)

	}

	if( test == "case2" & !is.null(Sigma0) ){

		k1 <- p*(p+1)/2
		C <- matrix(0, k1, k)
		C[,p+(1:k1)]<-diag(k1)
		#c <- matrixcalc::vech(Sigma0)

		mu0 <- apply(Y,2, mean)
		

	}


	I <- matrix(0, k, k)
	Sigma0.inv <- solve(Sigma0)
	I[1:p, 1:p] <- Sigma0.inv
	Dp <- Dmatrix(p)
	I[(p+1):k, (p+1):k] <- 1/2*t(Dp)%*%( Sigma0.inv %x% Sigma0.inv )%*%Dp

	I.inv <- solve(I)

	e <- Y - matrix( mu0, n, p, byrow=TRUE )
	
	U.mu <-  Sigma0.inv %*% colSums(e)

	Sigma.e.e.Sigma <- sapply( 1:n, function(i) (Sigma0.inv%*%e[i,])%*%( t(e[i,])%*%Sigma0.inv ))

	U.phi <- 1/2 * t(Dp)  %*% rowSums(Sigma.e.e.Sigma - c(Sigma0.inv) )

	U <- as.matrix( c(U.mu, U.phi) )	

	temp <- t(U)%*%I.inv%*%t(C)

	C.alpha <- c(1/n * temp %*% solve(C%*%I.inv%*%t(C)) %*% t(temp))
	p.value <- pchisq(C.alpha, k1, lower.tail=FALSE)
	
	list( sv=list(mu=mu0, Sigma=Sigma0), 
    C.alpha=C.alpha, p.value=p.value, test=test)


}