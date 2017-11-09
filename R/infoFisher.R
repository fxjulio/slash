#' @title Fisher Information Matrix Slash
#'
#' @description
#' Obtain a Information matrix values for parameters.
#' 
#' @details
#' Evaluate analitic expected values.
#'
#' @param Y Sample slash.
#' @param mu Mean vector with length \code{p}.
#' @param Sigma Covariance matrix with size \code{p}x\code{p}.
#' @param eta  Parameter between 0 and 1.
#' @param block Returns a blocked matrix.
#'
#' @return
#' A \code{list} with \code{I.mu.mu}, \code{I.phi.phi} and \code{U.eta} components.
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' sv <- startValues.slash(Y)
#' fit <- em.slash( Y, sv$mu0 , sv$Sigma0, sv$eta0, verbose=FALSE )
#' infoFisher.slash(Y, fit$mu, fit$Sigma, fit$eta)
#' 
#' @export
infoFisher.slash <- function(Y, mu, Sigma, eta, block=TRUE){
	Y <- t(Y)
	n <- ncol(Y)
	p <- nrow(Y)
	mu <- as.matrix(mu)
	Sigma.inv <- solve(Sigma)
	
	if( eta < 1e-3 ){ # info Normal
		I.mu.mu <- Sigma.inv
		Dp <- Dmatrix(p)
		I.phi.phi <- 1/2*t(Dp)%*%( Sigma.inv %x% Sigma.inv )%*%Dp
		I.phi.eta <- matrix(0,p*(p+1)/2,1)
		I.eta.eta <- 0
	} else {

	dg <- dg.slash( eta, p )
	fg <- fg.slash( eta, p )

	I.mu.mu <- 4*dg/p*Sigma.inv  

	Dp <- Dmatrix(p)
	Kp <- Kmatrix(p)
	Np <- 1/2*(diag(p^2) + Kp)

	I.phi.phi <- t(Dp)%*%( 2*fg/p/(p+2)*(Sigma.inv %x% Sigma.inv)%*%Np + (fg/p/(p+2) - 1/4)*matrixcalc::vec(Sigma.inv)%*%t(matrixcalc::vec(Sigma.inv)) ) %*% Dp 

	hg <- hg.slash(eta, p)

	I.phi.eta <- hg/p * t(Dp) %*% vec(Sigma.inv)

	tg1 <- tg1.slash(eta, p)
	tg2 <- tg2.slash(eta, p)
	#tg <- tg.slash(eta, p)

	c.eta <- 1/(1-eta)
	I.eta.eta <- 1/eta^2 - p/eta*c.eta - (2/eta - p * c.eta )*tg1 + p^2/4 * c.eta^2 + tg2 
	#I.eta.eta <- -1/eta^2 - p/2*c.eta^2 - tg
	
	}# eta 

	if( block )
		list(I.mu.mu=I.mu.mu, I.phi.phi=I.phi.phi, I.phi.eta=I.phi.eta, I.eta.eta = I.eta.eta)
	else{
		npar <- p + p*(p+1)/2 + 1
		I <- matrix(0, npar, npar)
		I[1:p,1:p] <- I.mu.mu
		I[(p+1):(npar-1), (p+1):(npar-1)] <- I.phi.phi
		I[(p+1):(npar-1), npar] <- I.phi.eta
		I[ npar, (p+1):(npar-1)] <- t(I.phi.eta)
		I[ npar, npar ] <- I.eta.eta
		
		I
	}
		
}