#' @title Score Slash
#'
#' @description
#' Obtain a score values for parameters.
#' 
#' @details
#' Evaluate analitic loglik derivatives.
#'
#' @param Y Sample slash.
#' @param mu Mean vector with length \code{p}.
#' @param Sigma Covariance matrix with size \code{p}x\code{p}.
#' @param eta  Parameter between 0 and 1.
#' @param block Returns a blocked score.
#' @param hessian Logical. Is the hessian matrix evaluated?
#'
#' @return
#' A \code{list} with \code{U.mu}, \code{U.Sigma} and \code{U.eta} components if \code{block=FALSE} and \code{hessian=FALSE}
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' sv <- startValues.slash(Y)
#' fit <- em.slash( Y, sv$mu0 , sv$Sigma0, sv$eta0, verbose=FALSE )
#' qqplot.slash(Y, fit$mu, fit$Sigma, fit$eta)
#' loglik.slash(Y, fit$mu , fit$Sigma, fit$eta )
#' score.slash(Y, fit$mu , fit$Sigma, fit$eta )
#' 
#' @export
score.slash <- function(Y, mu, Sigma, eta, block=TRUE, hessian=FALSE){
	Y <- t(Y)
	n <- ncol(Y)
	p <- nrow(Y)
	mu <- as.matrix(mu)
	Sigma.inv <- solve(Sigma)

	Dp <- Dmatrix(p)
	
	if( eta < 1e-3 ){

		e <- Y - matrix(mu,p,n)
		U.mu <- Sigma.inv %*% rowSums(e)
		Sigma.e.e.Sigma <- sapply( 1:n, function(i) (Sigma.inv%*%e[,i])%*%( t(e[,i])%*%Sigma.inv ))
		U.phi <- 1/2 * t(Dp)  %*% rowSums(Sigma.e.e.Sigma - c(Sigma.inv) )
		U.eta <- NA
	
	
	} else { # eta > 1e-3
	
	e <- lapply(1:n, function(i) Y[,i,drop=FALSE]-mu )
	delta <- sapply(1:n, function(i) t(e[[i]])%*%Sigma.inv%*%(e[[i]]) )

	c.eta <- (1-eta)^(-1)
	a <- p/2 + 1/eta
	
	wG.delta <- wGdelta(delta, eta, p )


	#dummy.eta <- function(eta, delta) pgamma(1, p/2 + 1/eta, 1/2*(1-eta)^(-1)*delta, log=TRUE)+ lgamma(p/2 + 1/eta) - (p/2 + 1/eta)*log(1/2*(1-eta)^(-1)*delta)
	#wG.eta <- sapply( 1:n, function(i) numDeriv::grad(dummy.eta, eta, delta=delta[i]) )
	wG.eta <- wGeta(delta, eta, p)

	U.mu <- rowSums( sapply(1:n, function(i) -2*wG.delta[i] * Sigma.inv%*%e[[i]]) )
	
	Sigma.e.e.Sigma <-  lapply(1:n, function(i) 
		 (Sigma.inv%*% e[[i]]) %*% (t(e[[i]]) %*% Sigma.inv) )  # (Ax)(x'A) es mas rapido que Axx'A
		 
	temp <- sapply( 1:n , function(i) wG.delta[i]*Sigma.e.e.Sigma[[i]] + 1/2*Sigma.inv)
	U.phi <- -t(Dp)%*%rowSums(temp)

	#temp <- sapply( 1:n, function(i) -t(Dp)%*%matrixcalc::vec(wG.delta[i]*Sigma.e.e.Sigma[[i]] + 1/2*Sigma.inv)   )
	#U.phi <- rowSums(temp)


	U.eta <- sum(-1/eta + p/2 * c.eta + wG.eta)
	
	}

	
	if( hessian ){
		wG.delta.delta <- wGdeltadelta(delta, eta, p)

		temp <- lapply(1:n, function(i)  
			2*wG.delta[i] * Sigma.inv + 4*wG.delta.delta[i]* Sigma.e.e.Sigma[[i]]   )

		L.mu.mu <- sumList(temp)/n

		temp <- lapply(1:n, function(i)  2*(    wG.delta[i]*(t(e[[i]])%x% diag(p))%*%( (Sigma.inv%x%Sigma.inv)%*%Dp ) +
				wG.delta.delta[i] * (Sigma.inv%*%e[[i]]) %*% ( t(vec(Sigma.e.e.Sigma[[i]]))%*%Dp  )      ) )

		L.mu.phi <- sumList(temp)/n

		aux <- (Sigma.inv%x%Sigma.inv)%*%Dp

		temp <- lapply(1:n, function(i) 
			1/2*t(Dp)%*% ( Sigma.inv%x%Sigma.inv )%*%Dp +
			2*wG.delta[i]*t(Dp)%*%( Sigma.inv%x%Sigma.e.e.Sigma[[i]] )%*%Dp+
			t(aux)%*%((e[[i]]%*%t(e[[i]]))%x%(e[[i]]%*%t(e[[i]])))%*%aux )
		
		L.phi.phi <- sumList(temp)/n
		
		wG.delta.eta <- wgEtaDelta(delta, eta, p)
		
		temp <- lapply(1:n, function(i)	-2*wG.delta.eta[i]*Sigma.inv%*%e[[i]] )
		L.mu.eta <- sumList(temp)/n
		
		temp <- lapply(1:n, function(i)	-wG.delta.eta[i]*t(Dp)%*%vec( Sigma.e.e.Sigma[[i]] ) )
		L.phi.eta <- sumList(temp)/n
		
		temp <- 1/eta^2 + p/2*c.eta^2 + wGetaeta( delta, eta, p )
		L.eta.eta <- sum(temp)/n
		
	}

	if( block )
		if( hessian )
			list(U=list(U.mu=U.mu, U.phi=U.phi, U.eta=U.eta), 
				L=list(L.mu.mu=L.mu.mu, L.mu.phi=L.mu.phi, L.mu.eta=L.mu.eta,
						L.phi.phi=L.phi.phi, L.phi.eta=L.phi.eta, L.eta.eta = L.eta.eta))
		else
			list(U.mu=U.mu, U.phi=U.phi, U.eta=U.eta)
	else{
		as.matrix( c(U.mu, U.phi, U.eta) )
	}
	
}
