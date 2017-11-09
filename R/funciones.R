#' Slash reparametrized
#'
#' The package \pkg{slash} implements routines to EM estimate.
#'
#' Several utilities for ML estimators.
#'
#' @docType package
#' @name slash
#' @import matrixcalc numDeriv gsl
#' @useDynLib slash
NULL





sumList <- function( Y ){
	if(length(Y) == 1 ) return( Y[[1]] )
	sum <- Y[[1]]
	for( i in 2:length(Y) )
		sum <- sum + Y[[i]]
	sum
}

P <- function(a, b, log=FALSE) pgamma(1, a, b, log=log)


dlP.eta <- function( p, eta, delta ){
	
	sapply( 1:length(delta), function(i){
		dummy <- function(x) pgamma(1, p/2 + 1/x, (1-x)^(-1) * delta[i] / 2, log=TRUE  )
		numDeriv::grad( dummy, eta )
	})
}

#' @title Log-likelihood of multivariate slash reparemetrized
#'
#' @description
#' Obtain a log-likelihood of a sample Slash.
#' 
#' @details
#' Compute the log-likelihood without conditioning.
#'
#' @param Y Sample slash.
#' @param mu Mean vector with length \code{p}.
#' @param Sigma Covarianze matrix with size \code{p}x\code{p}.
#' @param eta  Parameter between 0 and 1.
#'
#' @return
#' A \code{numeric}.
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' loglik.slash(Y, rep(0, p), diag(p), 1/2 )
#' 
#' @export
loglik.slash <- function( Y, mu, Sigma, eta ){
	Y <- t(Y)
	n <- ncol(Y)
	p <- nrow(Y)
	mu <- as.matrix(mu)

	Sigma.inv <- solve(Sigma)
	e <- lapply(1:n, function(i) Y[,i,drop=FALSE]-mu )
	delta <- sapply(1:n, function(i) t(e[[i]])%*%Sigma.inv%*%(e[[i]]) )

	c.eta <- (1-eta)^(-1)
	a <- p/2 + 1/eta
	b <- 1/2 * c.eta * delta
	lKp.eta <- -log(eta)- p/2*log(1-eta) - p/2*log(2*pi)

	lG.ab <- pgamma(1, a, b, log=TRUE) + lgamma(a) - a*log(b)
	detSigma <- determinant(Sigma, logarithm = TRUE)$modulus
	attr(detSigma, "logarithm") <- NULL

	n*lKp.eta -n/2*detSigma  + sum(lG.ab)
}

#' @title Start values for EM algotithm
#'
#' @description
#' Obtain start values for EM algoritm plug-in estimators of mean and covariance matrix. \emph{eta} is 
#' calculate from coeficient of multivariate kurtosis
#' 
#' @details
#' Compute Start values for EM algotithm.
#'
#' @param Y Sample slash.
#' @param mu0 Specify mu0. Default is \code{NULL}.
#'
#' @return
#' A \code{list} with \code{mu0}, \code{Sigma0} and \code{eta0} components.
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' startValues.slash(Y)
#' 
#' @export
startValues.slash <- function(Y, mu0=NULL, Sigma0=NULL, eta0=NULL){
	n <- nrow(Y)
	p <- ncol(Y)
	mu0 <- if(is.null(mu0)) apply(Y, 2, mean) else mu0
	e <- Y-matrix(mu0, n, p, byrow=TRUE)
	#Sigma0 <- if(is.null(Sigma0)) 1/(n-1)*sumList(lapply(1:n, function(i)  t(e[i,,drop=FALSE])%*%e[i,,drop=FALSE] )) else Sigma0
	Sigma0 <- if(is.null(Sigma0)) cov(e) else Sigma0
	Sigma0.inv <- solve(Sigma0)
	
	profileloglik <- function(x) loglik.slash(Y, mu0, Sigma0, x)
	eta0ml <- optimize( profileloglik, c(0,1),  maximum=TRUE,
				tol = .Machine$double.eps)$maximum
	
	# scoreeta <- function(x) score.slash(Y, mu0, Sigma0, x)$U.eta
	# eta0score <- tryCatch(
		# uniroot( scoreeta, eta0ml*c(0.8,1.2), tol = .Machine$double.eps)$root,
		# error = function(e) NA)
	
	
	#beta.2p <- sum(sapply(1:n, function(i) (e[i,,drop=FALSE]%*%Sigma0.inv%*%t(e[i,,drop=FALSE]))^2 ))/n
	beta.2p <- mean( mahalanobis( Y, mu0, Sigma0 )^2 )
	
	#a <- beta.2p/p/(p+2)
	a <- beta.2p - p*(p+2)

	#if( is.null(eta0) & a > 1){
	#	eta0 <- 1-a + sqrt(a*(a-1))
	#}
	if( is.null(eta0) & a > 0){
		eta0 <- (-a + sqrt(beta.2p*a))/(p*(p+2))
	}
	
	
	if( is.null(eta0) ){
		cat("Using pseudo ML eta estimator\n")
		eta0 <- NA
	}
	
	list( mu0=mu0, Sigma0=Sigma0,
		eta0=eta0, eta0ml=eta0ml, 
		#eta0score=eta0score, 
		beta.2p=beta.2p )
}

#' @title EM algotithm
#'
#' @description
#' Obtain MLE values from EM algoritm.
#' 
#' @details
#' Conditioning \emph{v} in step E. Use \code{startValues.slash} function to obtain start values
#'
#' @param Y Sample slash.
#' @param mu0 Mean vector with length \code{p}.
#' @param Sigma0 Covariance matrix with size \code{p}x\code{p}.
#' @param eta0  Parameter between 0 and 1.
#' @param relError Relative error.
#' @param verbose Logical. \code{TRUE} indicate trace the steps.
#'
#' @return
#' A \code{list} with \code{mu0}, \code{Sigma0} and \code{eta0} components.
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' sv <- startValues.slash(Y)
#' em.slash( Y, sv$mu0 , sv$Sigma0, sv$eta0, verbose=TRUE )
#' 
#' @export
em.slash <- function( Y, mu0, Sigma0, eta0, relError=1e-5, verbose=FALSE ){
	Y <- t(Y)
	p <- nrow(Y)
	n <- ncol(Y)
	mu0 <- as.matrix(mu0)
	
	flag <- TRUE
	iter <- 1

	max.iter <- 1000

	while(flag){
	c.eta <- 1/(1-eta0)
	Sigma.inv <- solve(Sigma0)

	e <- lapply(1:n, function(i) Y[,i,drop=FALSE]-mu0 )
	delta <- sapply(1:n, function(i) t(e[[i]])%*%Sigma.inv%*%(e[[i]]) )

	lnumP <- P( p/2+1/eta0 + 1, c.eta*delta/2, log=TRUE )
	ldenP <- P( p/2+1/eta0, c.eta*delta/2, log=TRUE)
	w <- (p*eta0 + 2)/(eta0*delta*c.eta)*exp(lnumP-ldenP)

	ci <- sapply(1:n, function(i){
		numI <- function(v) log(v)*v^(p/2 + 1/eta0 -1)*exp(-1/2*c.eta*delta[i] * v)
		denI <- function(v) v^(p/2 + 1/eta0 -1)*exp(-1/2*c.eta*delta[i] * v)
		integrate(numI, 0, 1)$value / integrate(denI, 0, 1)$value
	})

	myfun <- function( eta, delta ) P( p/2+1/eta, 1/(1-eta)*delta/2, log=TRUE)
	dlogP.deta <- sapply(delta, function(delta) grad( myfun, eta0, delta=delta ) )

	ci.asymptotic <- 1/eta0^2*( log(c.eta * delta / 2 ) - digamma(p/2 + 1/eta0) ) - (p/2 + 1/eta0)*c.eta + dlogP.deta

	if( any( is.nan(ci) )  )
		cat("eta:", eta0, "some NaN in c_i\n")

	mu1 <- rowSums(rep(w, each=p)*Y)/sum(w)

	temp <- lapply(1:n, function(i) c.eta*w[i]*(Y[,i,drop=FALSE]-mu0)%*%t(Y[,i,drop=FALSE]-mu0) )
	Sigma1 <- sumList(temp)/n
	
	sci <- sum(ci, na.rm = TRUE)

	feta <- function(x) -n*log(x)-n*p/2*log(1-x)+(p/2+1/x-1)*sci - 1/2*sum(w*delta)/(1-x)
	eta1 <- optimize( feta , c(0, 1),  maximum = TRUE )$maximum

	crit <- c()
	crit[1] <- max(abs(mu0-mu1))/max(abs(mu1))
	crit[2] <- max(abs(Sigma0 - Sigma1))/max(abs(Sigma1))
	crit[3] <- abs(eta0 - eta1)/eta1

	flag <- any(crit > relError)

	mu0 <- mu1
	Sigma0 <- Sigma1
	eta0 <- eta1
	iter <- iter+1

	if(verbose){
		curve( feta, main=paste("eta.max:", eta1))
		cat("iter:", iter,  "crit1: ", crit[1], "crit2: ", crit[2], "crit3: ", crit[3], "\n")
		flush.console()
	}
	if(iter > max.iter){
		cat("max.iter reached\n")
		break
		}
	}


	list( mu=mu1, Sigma=Sigma1, eta=eta1 )
}

#' @title qq-plot Slash
#'
#' @description
#' Obtain a qq-plot for \emph{delta}.
#' 
#' @details
#' Use distribution of \emph{delta} values
#'
#' @param Y Sample slash.
#' @param mu Mean vector with length \code{p}.
#' @param Sigma Covariance matrix with size \code{p}x\code{p}.
#' @param eta  Parameter between 0 and 1.
#'
#' @return
#' A qq-plot.
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' sv <- startValues.slash(Y)
#' fit <- em.slash( Y, sv$mu0 , sv$Sigma0, sv$eta0, verbose=FALSE )
#' qqplot.slash(Y, fit$mu, fit$Sigma, fit$eta)
#' 
#' @export
qqplot.slash <- function(Y, mu, Sigma, eta){
	Y <- t(Y)
	n <- ncol(Y)
	p <- nrow(Y)
	Sigma.inv <- solve(Sigma)
	delta <- sapply( 1:n, function(i)  t(Y[,i,drop=FALSE]-mu)%*%Sigma.inv%*%(Y[,i,drop=FALSE]-mu)   )

	qtl <- ppoints(n)
	thq <- qdelta.slash(qtl, eta, p)

	thq <- thq[order(order(delta))]

	plot(delta, thq)
	y <- quantile(delta, c(0.25, 0.75))
	x <- qdelta.slash(c(0.25, 0.75), eta, p)

	slope <- diff(x)/diff(y)
	int <- x[1L]-slope*y[1L]

	points(y, x, col=2)

	abline(int, slope)
}




#' @title Distribution, quantile and density of delta
#'
#' @description
#' Obtain distribution, quantile or density.
#' 
#' @details
#' Evaluate analitic expected values.
#'
#' @param x Value of \emph{delta}.
#' @param qtl Quantile value. Between 0 and 1.
#' @param eta  Parameter between 0 and 1.
#' @param p Dimension of \code{Y}.
#'
#' @return
#' A distribution, density or quantile value.
#'
#' @examples
#' n <- 1000
#' eta <- 1/3
#' p <- 4
#' Z <- sim.slash( n, rep(0,p), diag(p), eta )
#' delta <- sapply(1:n, function(i) Z[i,,drop=FALSE]%*%t(Z[i,,drop=FALSE]) )
#'
#' hist(delta, prob=TRUE)
#' curve(ddelta.slash(x, eta, p), add=TRUE, col=2)
#' 
#' @name delta
#' @export pdelta.slash qdelta.slash ddelta.slash
NULL

#' @rdname delta
pdelta.slash <- function(x, eta, p) {
  if(eta > 1e-4 )
    return(pchisq(x/(1-eta), p)  - exp( 1/eta*log(2) + lgamma(p/2+1/eta) - 1/eta*log(x/(1-eta)) - lgamma(p/2)  )*pchisq(x/(1-eta), p + 2/eta))
  else 
    return( pchisq(x, p)   )
  
}

#' @rdname delta
qdelta.slash <- function(qtl, eta, p) sapply( 1:length(qtl), function(i){
	xinf <- exp( log(2) + eta*lgamma(p/2 + 1/eta) -eta*lgamma(p/2) - 709*eta)
	xsup <- exp( log(2) + eta*lgamma(p/2 + 1/eta) -eta*lgamma(p/2) + 709*eta)
	uniroot( function(x) pdelta.slash(x, eta, p)-qtl[i], c(xinf, xsup)   )$root
	})

#' @rdname delta
ddelta.slash <- function(x, eta, p){
	ceta <- 1/(1-eta)
	a <- p/2+1/eta
	if( eta > 1e-4 ){
		aux <- log(ceta/2/eta) + lgamma(a) - lgamma(p/2)
		out <- sapply( x, function(x){
			b <- ceta/2*x
			P <- gamma_inc_P(a, b)
			aux2 <- exp(aux - (1/eta+1)*log(b))
			if( P > 0 & is.finite(aux2) )	aux2*P else dchisq(x, p)
		})
	}
	else
		out <- dchisq(x, p)
	
	out
}

#' @title Evaluate integrals \emph{dg} and \emph{fg}.
#'
#' @description
#' Obtain \emph{dg} and \emph{fg} values from numerical integrals.
#' 
#' @details
#' Evaluate integrals \emph{dg} and \emph{fg}.
#'
#' @param x Value of \emph{delta}.
#' @param qtl Quantile value. Between 0 and 1.
#' @param eta  Parameter between 0 and 1.
#' @param p Dimension of \code{Y}.
#'
#' @return
#' A \code{numeric}.
#'
#' @examples
#' p <- 4
#' etas <- seq(1/50, 1-1/50, length.out=200)
#' dgs <- sapply( etas, function(eta) dg.slash( eta, p ) )
#' fgs <- sapply( etas, function(eta) fg.slash( eta, p ) )
#'
#' plot( etas, dgs, xlim=c(0,1), type="l")
#' abline( h=p/4, col=2 )
#'
#' plot( etas, fgs, xlim=c(0,1), type="l")
#' abline( h=p*(p+2)/4, col=2 )
#' 
#' @name dgfg
#' @export dg.slash fg.slash
NULL

#' @rdname dgfg
dg.slash <- function(eta, p){
	if( eta < 1e-3 )  return(p/4)

	myfunc <- function(delta){
		wGdelta(delta, eta, p)^2 * delta * ddelta.slash( delta, eta, p)
	}

	tryCatch( integrate(myfunc, 0, Inf)$value, error=function(e){
		cat("dg eta:",eta,"p:",p,"a:",p/2+1/eta,"\n")
		NA
	})

}

#' @rdname dgfg
fg.slash <- function(eta, p){
	if( eta < 1e-3 ) return( p*(p+2)/4 )

	myfunc <- function(delta){
		wGdelta(delta, eta, p)^2 * delta^2 * ddelta.slash( delta, eta, p)
	}

	tryCatch( integrate(myfunc, 0, Inf)$value, error=function(e){
		cat("fg eta:",eta,"p:",p,"a:",p/2+1/eta,"\n")
		NA
	})

}

#' @title Expected values \emph{hg} and \emph{tg}
#' @name hgtg
#' @export hg.slash tg1.slash tg2.slash
NULL

#' @rdname hgtg
hg.slash <- function(eta, p){

	myfunc <- function(delta){
		wgEtaDelta(delta, eta, p) * delta * ddelta.slash( delta, eta, p)
	}
	
	tryCatch( integrate(myfunc, 0, Inf)$value, error=function(e){
		cat("hg eta:",eta,"p:",p,"a:",p/2+1/eta,"\n")
		NA
	})

}


#' @rdname hgtg
tg1.slash <- function(eta, p){

	myfunc <- function(delta){
		wGeta(delta, eta, p) * ddelta.slash( delta, eta, p)
	}

	tryCatch( integrate(myfunc, 0, Inf)$value, error=function(e){
		cat("tg1 eta:",eta,"p:",p,"a:",p/2+1/eta,"\n")
		NA
	})

}

#' @rdname hgtg
tg2.slash <- function(eta, p){

	myfunc <- function(delta){
		wGeta(delta, eta, p)^2 * ddelta.slash( delta, eta, p)
	}

	tryCatch( integrate(myfunc, 0, Inf)$value, error=function(e){
		cat("tg2 eta:",eta,"p:",p,"a:",p/2+1/eta,"\n")
		NA
	})

}

tg.slash <- function(eta, p){

	myfunc <- function(delta){
		wGetaeta(delta, eta, p) * ddelta.slash( delta, eta, p)
	}
	
	tryCatch( integrate(myfunc, 0, Inf)$value, error=function(e){
		cat("tg eta:",eta,"p:",p,"a:",p/2+1/eta,"\n")
		NA
	})
}

wGeta <- function(delta, eta, p){

#	myfunc <- function(eta, delta) pgamma(1, p/2 + 1/eta, 1/2*(1-eta)^(-1)*delta, log=TRUE)+ lgamma(p/2 + 1/eta) - (p/2 + 1/eta)*log(1/2*(1-eta)^(-1)*delta)
#	sapply( delta, function(delta) grad(myfunc, eta, delta=delta) )
	a <- p/2 + 1/eta
	ceta <- 1/(1-eta)
  
	
	out <- sapply( ceta*delta/2, function(b){
		if( b <= a  )
      return(1/Gserie(a,b)*(-1/eta^2*H(a,b)-b*ceta*Gserie(a+1,b)))
    else {
      I <- digami(b,a)
      return( -1/eta^2*( digamma(a) + I[3]/I[6] - log(b) ) + 
        b*ceta*(I[1]/I[6] - a/b) )
    }
	})
	
	# f0 <- 1/eta^2/a 
	# f1 <- ceta/2*( 1/eta^2/(a+1) - a/eta^2/(a+1)^2 - a*ceta/(a+1) )
	# f2 <- (ceta/2)^2*(2*(a^3*ceta*eta^2 + 3*a^2*ceta*eta^2 + a^2 + 2*a*ceta*eta^2 + a - 1)/(eta^2*(a^5 + 7*a^4 + 19*a^3 + 25*a^2 + 16*a + 4)))
	# f3 <- (ceta/2)^3*( 6*(a^5*ceta*eta^2 + 5*a^4*ceta*eta^2 + a^4 + 5*a^3*ceta*eta^2 + 2*a^3 - 5*a^2*ceta*eta^2 - 5*a^2 - 6*a*ceta*eta^2 - 8*a + 2)/(eta^2*(a^8 + 14*a^7 + 83*a^6 + 272*a^5 + 539*a^4 + 662*a^3 + 493*a^2 + 204*a + 36)))
	# f4 <- (ceta/2)^4*( 24*(a^8*ceta*eta^2 + 9*a^7*ceta*eta^2 + a^7 + 19*a^6*ceta*eta^2 + 5*a^6 - 43*a^5*ceta*eta^2 - 11*a^5 - 216*a^4*ceta*eta^2 - 96*a^4 - 254*a^3*ceta*eta^2 - 138*a^3 - 44*a^2*ceta*eta^2 + 33*a^2 + 48*a*ceta*eta^2 + 114*a - 12)/(eta^2*(a^12 + 25*a^11 + 279*a^10 + 1837*a^9 + 7945*a^8 + 23775*a^7 + 50477*a^6 + 76631*a^5 + 82602*a^4 + 61700*a^3 + 30344*a^2 + 8832*a + 1152)) )
	# f5 <- (ceta/2)^5*120*(a^10*ceta*eta^2 + 11*a^9*ceta*eta^2 + a^9 + 14*a^8*ceta*eta^2 + 6*a^8 - 264*a^7*ceta*eta^2 - 36*a^7 - 1323*a^6*ceta*eta^2 - 356*a^6 - 2121*a^5*ceta*eta^2 - 771*a^5 - 64*a^4*ceta*eta^2 + 334*a^4 + 2614*a^3*ceta*eta^2 + 2510*a^3 + 1372*a^2*ceta*eta^2 + 1264*a^2 - 240*a*ceta*eta^2 - 984*a + 48)/(eta^2*(a^15 + 36*a^14 + 589*a^13 + 5806*a^12 + 38542*a^11 + 182440*a^10 + 636002*a^9 + 1662628*a^8 + 3286613*a^7 + 4914332*a^6 + 5515889*a^5 + 4567166*a^4 + 2702844*a^3 + 1080392*a^2 + 261120*a + 28800))
	
	# flag <- !is.finite(out) | delta^5*abs(f5)/120 < 1e-5
	
	# out[flag] <- f0 + delta[flag]*f1 + delta[flag]^2*f2/2 + delta[flag]^3*f3/6 +  delta[flag]^4*f4/24 
	out
}

lGserie <- function (a, b) 
{
    if (b <= a  ) {
        out <- 1/a
        flag <- TRUE
        i <- 1
        while (flag) {
		aux <- sapply( i, function(i) i*log(b) - sum(log(a+(0:i))) )
            flag <- exp(aux)/out > .Machine$double.eps/2
            i <- i + 1
            out <- out + exp(aux)
        }
        return( -b + log(out) )
    }
    else {
        return( -a*log(b) + lgamma(a) )
    }
}


Gserie <- function(a,b){
	
	if( b <= a  ){

	out <- exp(-b)/a
  if( out == 0 ){
    warning("Gserie a:", a," b:", b)
    return(NA)
  }
	flag <- TRUE
	i <- 1
	while( flag ){
		aux <- sapply( i, function(i) exp(-b)*b^i/prod(a + (0:i)) )
		flag <- aux/out  > .Machine$double.eps/2
		i <- i + 1
		out <- out + aux 
	}
		return(out)
	}
	
	else {
		return(  exp(-a*log(b)+lgamma(a)) )
	}
}

G <- function(a,b) integrate( function(v) v^(a-1)*exp(-b*v), 0, 1)$value
H <- function(a,b) integrate( function(v) log(v)*v^(a-1)*exp(-b*v), 0, 1)$value
F <- function(a,b) integrate( function(v) log(v)^2*v^(a-1)*exp(-b*v), 0, 1)$value

wGetaeta <- function(delta, eta, p){

#	myfunc <- function(eta, delta) pgamma(1, p/2 + 1/eta, 1/2*(1-eta)^(-1)*delta, log=TRUE)+ lgamma(p/2 + 1/eta) - (p/2 + 1/eta)*log(1/2*(1-eta)^(-1)*delta)
#	sapply( delta, function(delta) hessian(myfunc, eta, delta=delta) )
	a <- p/2 + 1/eta
	ceta <- 1/(1-eta)
	
	f0 <- -(2*a*eta - 1)/(a**2*eta**4)
	f1 <- (ceta/2)*(-2*(a**3*ceta**2*eta**4 + 3*a**2*ceta**2*eta**4 + 3*a*ceta**2*eta**4 - a*ceta*eta**2 + a*eta + ceta**2*eta**4 - ceta*eta**2 + eta - 1)/(eta**4*(a + 1)**3))
	f2 <- (ceta/2)^2*(2*(a**5*ceta**2*eta**4 + 6*a**4*ceta**2*eta**4 + 4*a**4*ceta*eta**2 - 2*a**4*eta + 13*a**3*ceta**2*eta**4 + 16*a**3*ceta*eta**2 - 8*a**3*eta + 3*a**3 + 12*a**2*ceta**2*eta**4 + 16*a**2*ceta*eta**2 - 8*a**2*eta + 6*a**2 + 4*a*ceta**2*eta**4 - 4*a*ceta*eta**2 + 2*a*eta - 4*a - 8*ceta*eta**2 + 4*eta - 10)/(eta**4*(a + 1)**4*(a + 2)**3))
	f3 <- (ceta/2)^3*(12*(a**8*ceta**2*eta**4 + 11*a**7*ceta**2*eta**4 + 3*a**7*ceta*eta**2 - a**7*eta + 46*a**6*ceta**2*eta**4 + 24*a**6*ceta*eta**2 - 8*a**6*eta + 2*a**6 + 86*a**5*ceta**2*eta**4 + 54*a**5*ceta*eta**2 - 18*a**5*eta + 10*a**5 + 49*a**4*ceta**2*eta**4 - 30*a**4*ceta*eta**2 + 10*a**4*eta - 4*a**4 - 61*a**3*ceta**2*eta**4 - 267*a**3*ceta*eta**2 + 89*a**3*eta - 94*a**3 - 96*a**2*ceta**2*eta**4 - 318*a**2*ceta*eta**2 + 106*a**2*eta - 152*a**2 - 36*a*ceta**2*eta**4 - 78*a*ceta*eta**2 + 26*a*eta - 28*a + 36*ceta*eta**2 - 12*eta + 58)/(eta**4*(a + 1)**5*(a + 2)**3*(a + 3)**3))
	f4 <- (ceta/2)^4*(24*(3*a**12*ceta**2*eta**4 + 57*a**11*ceta**2*eta**4 + 8*a**11*ceta*eta**2 - 2*a**11*eta + 432*a**10*ceta**2*eta**4 + 120*a**10*ceta*eta**2 - 30*a**10*eta + 5*a**10 + 1536*a**9*ceta**2*eta**4 + 592*a**9*ceta*eta**2 - 148*a**9*eta + 55*a**9 + 1479*a**8*ceta**2*eta**4 + 152*a**8*ceta*eta**2 - 38*a**8*eta + 90*a**8 - 8259*a**7*ceta**2*eta**4 - 9672*a**7*ceta*eta**2 + 2418*a**7*eta - 1294*a**7 - 35514*a**6*ceta**2*eta**4 - 41096*a**6*ceta*eta**2 + 10274*a**6*eta - 7952*a**6 - 63342*a**5*ceta**2*eta**4 - 75600*a**5*ceta*eta**2 + 18900*a**5*eta - 17946*a**5 - 56832*a**4*ceta**2*eta**4 - 55368*a**4*ceta*eta**2 + 13842*a**4*eta - 11631*a**4 - 19848*a**3*ceta**2*eta**4 + 17664*a**3*ceta*eta**2 - 4416*a**3*eta + 19581*a**3 + 4032*a**2*ceta**2*eta**4 + 48576*a**2*ceta*eta**2 - 12144*a**2*eta + 34596*a**2 + 3456*a*ceta**2*eta**4 + 17088*a*ceta*eta**2 - 4272*a*eta + 10848*a - 2304*ceta*eta**2 + 576*eta - 4944)/(eta**4*(a + 1)**6*(a + 2)**4*(a + 3)**3*(a + 4)**3))
	f5 <- (ceta/2)^5*(240*(2*a**15*ceta**2*eta**4 + 52*a**14*ceta**2*eta**4 + 5*a**14*ceta*eta**2 - a**14*eta + 528*a**13*ceta**2*eta**4 + 105*a**13*ceta*eta**2 - 21*a**13*eta + 3*a**13 + 2212*a**12*ceta**2*eta**4 + 695*a**12*ceta*eta**2 - 139*a**12*eta + 48*a**12 - 2688*a**11*ceta**2*eta**4 - 805*a**11*ceta*eta**2 + 161*a**11*eta + 107*a**11 - 76244*a**10*ceta**2*eta**4 - 37735*a**10*ceta*eta**2 + 7547*a**10*eta - 2796*a**10 - 397156*a**9*ceta**2*eta**4 - 239135*a**9*ceta*eta**2 + 47827*a**9*eta - 27341*a**9 - 1093924*a**8*ceta**2*eta**4 - 736295*a**8*ceta*eta**2 + 147259*a**8*eta - 107880*a**8 - 1672530*a**7*ceta**2*eta**4 - 1040175*a**7*ceta*eta**2 + 208035*a**7*eta - 164925*a**7 - 1023568*a**6*ceta**2*eta**4 + 262510*a**6*ceta*eta**2 - 52502*a**6*eta + 196716*a**6 + 858228*a**5*ceta**2*eta**4 + 3282370*a**5*ceta*eta**2 - 656474*a**5*eta + 1201702*a**5 + 1993712*a**4*ceta**2*eta**4 + 4646500*a**4*ceta*eta**2 - 929300*a**4*eta + 1766112*a**4 + 1271216*a**3*ceta**2*eta**4 + 2151080*a**3*ceta*eta**2 - 430216*a**3*eta + 621686*a**3 + 197760*a**2*ceta**2*eta**4 - 535680*a**2*ceta*eta**2 + 107136*a**2*eta - 777864*a**2 - 57600*a*ceta**2*eta**4 - 524640*a*ceta*eta**2 + 104928*a*eta - 513216*a + 28800*ceta*eta**2 - 5760*eta + 85152)/(eta**4*(a + 1)**7*(a + 2)**4*(a + 3)**3*(a + 4)**3*(a + 5)**3))
	f6 <- (ceta/2)^6*(720*(5*a**20*ceta**2*eta**4 + 195*a**19*ceta**2*eta**4 + 12*a**19*ceta*eta**2 - 2*a**19*eta + 3080*a**18*ceta**2*eta**4 + 396*a**18*ceta*eta**2 - 66*a**18*eta + 7*a**18 + 21660*a**17*ceta**2*eta**4 + 4596*a**17*ceta*eta**2 - 766*a**17*eta + 189*a**17 - 13240*a**16*ceta**2*eta**4 + 5844*a**16*ceta*eta**2 - 974*a**16*eta + 1267*a**16 - 1664610*a**15*ceta**2*eta**4 - 432372*a**15*ceta*eta**2 + 72062*a**15*eta - 15491*a**15 - 16457300*a**14*ceta**2*eta**4 - 5592552*a**14*ceta*eta**2 + 932092*a**14*eta - 358181*a**14 - 90848220*a**13*ceta**2*eta**4 - 36418248*a**13*ceta*eta**2 + 6069708*a**13*eta - 3079451*a**13 - 319263455*a**12*ceta**2*eta**4 - 139660248*a**12*ceta*eta**2 + 23276708*a**12*eta - 14153125*a**12 - 687385605*a**11*ceta**2*eta**4 - 274277004*a**11*ceta*eta**2 + 45712834*a**11*eta - 28619535*a**11 - 597615740*a**10*ceta**2*eta**4 + 117081660*a**10*ceta*eta**2 - 19513610*a**10*eta + 52120958*a**10 + 1315449960*a**9*ceta**2*eta**4 + 2406582804*a**9*ceta*eta**2 - 401097134*a**9*eta + 532248482*a**9 + 5633681610*a**8*ceta**2*eta**4 + 7218458052*a**8*ceta*eta**2 - 1203076342*a**8*eta + 1639737278*a**8 + 9654420060*a**7*ceta**2*eta**4 + 11072757012*a**7*ceta*eta**2 - 1845459502*a**7*eta + 2519432054*a**7 + 8913285360*a**6*ceta**2*eta**4 + 7737348576*a**6*ceta*eta**2 - 1289558096*a**6*eta + 1130470836*a**6 + 3475125600*a**5*ceta**2*eta**4 - 2110748832*a**5*ceta*eta**2 + 351791472*a**5*eta - 2631944040*a**5 - 1060045920*a**4*ceta**2*eta**4 - 8069797248*a**4*ceta*eta**2 + 1344966208*a**4*eta - 4695349152*a**4 - 1503455040*a**3*ceta**2*eta**4 - 4967097408*a**3*ceta*eta**2 + 827849568*a**3*eta - 2377625280*a**3 - 360806400*a**2*ceta**2*eta**4 - 296144640*a**2*ceta*eta**2 + 49357440*a**2*eta + 509235840*a**2 + 31104000*a*ceta**2*eta**4 + 458887680*a*ceta*eta**2 - 76481280*a*eta + 583073280*a - 12441600*ceta*eta**2 + 2073600*eta - 52427520)/(eta**4*(a + 1)**8*(a + 2)**5*(a + 3)**4*(a + 4)**3*(a + 5)**3*(a + 6)**3))

	out <- sapply( delta, function(delta){
		b <- ceta*delta/2
		-1/G(a,b)^2*(H(a,b)/eta^2+b*ceta*G(a+1,b))^2 + 
		1/G(a,b)*(2/eta^3*H(a,b)+1/eta^2*(1/eta^2*F(a,b)+b*ceta*H(a+1,b))-
		2*b*ceta^2*G(a,b)+b*ceta*(1/eta^2*H(a+1,b)+b*ceta*G(a+2,b)) )
	
		# I <- digami(b,a)
		
		# 2/eta^3*( digamma(a) + I[3]/I[6] - log(b)) +
		# 1/eta^4*( trigamma(a) - (I[3]/I[6])^2 + I[4]/I[6]  ) +
		# -2*b*c.eta/eta^2*( -I[1]*I[3]/I[6]^2 + I[5]/I[6] - 1/b ) +
		# 2*b*c.eta^2*( I[1]/I[6] - a/b ) + 
		# (b*c.eta)^2*( -(I[1]/I[6])^2 + I[2]/I[6] + a/b^2 )		
	})
	
	flag <- !is.finite(out) | delta^6*abs(f6)/720 < 1e-5
	out[flag] <- f0 + delta[flag]*f1 + delta[flag]^2*f2/2 +
		delta[flag]^3*f3/6 + delta[flag]^4*f4/24 + delta[flag]^5*f5/120
	out

}


wgEtaDelta <- function(delta, eta, p){
	#myfunc <- function(eta, delta, p) wGdelta(delta, eta, p)
	#sapply(delta, function(delta) grad( myfunc, eta, delta=delta, p=p ) )
	a <- p/2 + 1/eta
	ceta <- 1/(1-eta)
	out <- 	sapply(ceta* delta/2, function(b){
	if( b <= a )
		return( ceta/2*((-b*ceta*Gserie(a + 1, b) - H(a, b)/eta^2)*Gserie(a + 1, b)/Gserie(a, b)^2 + (b*ceta*Gserie(a + 2, b) - ceta*Gserie(a + 1, b) + H(a + 1, b)/eta^2)/Gserie(a, b)) )
	else {
		I <- digami(b, a)
        	return( ceta^2/2 * (I[1]/I[6] - a/b) + -1/eta^2 * ceta/2 * (-1/I[6]^2 * 
            I[1] * I[3] + I[5]/I[6] - 1/b) + b * ceta^2/2 * (-(I[1]/I[6])^2 + 
            I[2]/I[6] + a/b^2) )
	}
	})	
	
	# f0 <- (ceta/2)*(-(a^2*ceta*eta^2 + a*ceta*eta^2 - 1)/(eta^2*(a + 1)^2))
	# f1 <- (ceta/2)^2*(2*(a^3*ceta*eta^2 + 3*a^2*ceta*eta^2 + a^2 + 2*a*ceta*eta^2 + a - 1)/(eta^2*(a + 1)^3*(a + 2)^2))
	# f2 <- (ceta/2)^3*(6*(a^5*ceta*eta^2 + 5*a^4*ceta*eta^2 + a^4 + 5*a^3*ceta*eta^2 + 2*a^3 - 5*a^2*ceta*eta^2 - 5*a^2 - 6*a*ceta*eta^2 - 8*a + 2)/(eta^2*(a + 1)^4*(a + 2)^2*(a + 3)^2))
	# f3 <- (ceta/2)^4*(24*(a^8*ceta*eta^2 + 9*a^7*ceta*eta^2 + a^7 + 19*a^6*ceta*eta^2 + 5*a^6 - 43*a^5*ceta*eta^2 - 11*a^5 - 216*a^4*ceta*eta^2 - 96*a^4 - 254*a^3*ceta*eta^2 - 138*a^3 - 44*a^2*ceta*eta^2 + 33*a^2 + 48*a*ceta*eta^2 + 114*a - 12)/(eta^2*(a + 1)^5*(a + 2)^3*(a + 3)^2*(a + 4)^2))
	# f4 <- (ceta/2)^5*(120*(a^10*ceta*eta^2 + 11*a^9*ceta*eta^2 + a^9 + 14*a^8*ceta*eta^2 + 6*a^8 - 264*a^7*ceta*eta^2 - 36*a^7 - 1323*a^6*ceta*eta^2 - 356*a^6 - 2121*a^5*ceta*eta^2 - 771*a^5 - 64*a^4*ceta*eta^2 + 334*a^4 + 2614*a^3*ceta*eta^2 + 2510*a^3 + 1372*a^2*ceta*eta^2 + 1264*a^2 - 240*a*ceta*eta^2 - 984*a + 48)/(eta^2*(a + 1)^6*(a + 2)^3*(a + 3)^2*(a + 4)^2*(a + 5)^2))
	# f5 <- (ceta/2)^6*(720*(a^14*ceta*eta^2 + 18*a^13*ceta*eta^2 + a^13 + 63*a^12*ceta*eta^2 + 12*a^12 - 876*a^11*ceta*eta^2 - 44*a^11 - 10131*a^10*ceta*eta^2 - 1424*a^10 - 44172*a^9*ceta*eta^2 - 8871*a^9 - 81847*a^8*ceta*eta^2 - 19467*a^8 + 24060*a^7*ceta*eta^2 + 22586*a^7 + 384342*a^6*ceta*eta^2 + 195802*a^6 + 635874*a^5*ceta*eta^2 + 337556*a^5 + 312804*a^4*ceta*eta^2 + 90201*a^4 - 139704*a^3*ceta*eta^2 - 285180*a^3 - 121392*a^2*ceta*eta^2 - 169796*a^2 + 8640*a*ceta*eta^2 + 56640*a - 1440)/(eta^2*(a + 1)^7*(a + 2)^4*(a + 3)^3*(a + 4)^2*(a + 5)^2*(a + 6)^2))
	# f6 <- (ceta/2)^7*(5040*(a^16*ceta*eta^2 + 20*a^15*ceta*eta^2 + a^15 + 24*a^14*ceta*eta^2 + 13*a^14 - 2644*a^13*ceta*eta^2 - 123*a^13 - 30704*a^12*ceta*eta^2 - 3379*a^12 - 149126*a^11*ceta*eta^2 - 24415*a^11 - 241630*a^10*ceta*eta^2 - 54255*a^10 + 851386*a^9*ceta*eta^2 + 221681*a^9 + 5106891*a^8*ceta*eta^2 + 1660121*a^8 + 9966514*a^7*ceta*eta^2 + 3693954*a^7 + 5468822*a^6*ceta*eta^2 + 1162222*a^6 - 8453902*a^5*ceta*eta^2 - 8005562*a^5 - 12172668*a^4*ceta*eta^2 - 10382010*a^4 - 2151768*a^3*ceta*eta^2 + 459744*a^3 + 1869264*a^2*ceta*eta^2 + 4348488*a^2 - 60480*a*ceta*eta^2 - 642240*a + 8640)/(eta^2*(a + 1)^8*(a + 2)^4*(a + 3)^3*(a + 4)^2*(a + 5)^2*(a + 6)^2*(a + 7)^2))
	
	# flag <- !is.finite(out) | delta^6*abs(f6)/720 < 1e-5
	# out[flag] <- f0 + delta[flag]*f1 + delta[flag]^2*f2/2 +
		# delta[flag]^3*f3/6 + delta[flag]^4*f4/24 + delta[flag]^5*f5/120
	out
	
}


wGdelta <- function(delta, eta, p){
#	c.eta <- (1-eta)^(-1)
#	a <- p/2 + 1/eta
#   	-a/delta * exp( pgamma(1, a+1, c.eta*delta/2, log=TRUE) - pgamma(1, a, c.eta*delta/2, log=TRUE) )
	a <- p/2 + 1/eta
	c.eta <- 1/(1-eta)
	b <- c.eta*delta/2
	
	f0 <- (c.eta/2)*-a/(a+1)
	f1 <- (c.eta/2)^2*(a/((a + 1)^2*(a + 2)))
	f2 <- (c.eta/2)^3*(2*a*(a - 1)/((a + 1)^3*(a + 2)*(a + 3)))
	f3 <- (c.eta/2)^4*(6*a*(a^3 - a^2 - 6*a + 2)/((a + 1)^4*(a + 2)^2*(a + 3)*(a + 4)))
	f4 <- (c.eta/2)^5*(24*a*(a - 1)*(a^3 - 3*a^2 - 14*a + 2)/((a + 1)^5*(a + 2)^2*(a + 3)*(a + 4)*(a + 5)))
	f5 <- (c.eta/2)^6*(120*a*(a^7 - 3*a^6 - 49*a^5 - 57*a^4 + 222*a^3 + 264*a^2 - 198*a + 12)/((a + 1)^6*(a + 2)^3*(a + 3)^2*(a + 4)*(a + 5)*(a + 6)))
	f6 <- (c.eta/2)^7*(720*a*(a - 1)*(a^7 - 7*a^6 - 81*a^5 - 37*a^4 + 766*a^3 + 1048*a^2 - 390*a + 12)/((a + 1)^7*(a + 2)^3*(a + 3)^2*(a + 4)*(a + 5)*(a + 6)*(a + 7)))
	
	out <- -exp(log(gamma_inc_P(a+1,b)*a)-log(gamma_inc_P(a,b)*delta))
	flag <- !is.finite(out) | delta^6*abs(f6)/720 < 1e-5
	
	out[flag] <- f0 + delta[flag]*f1 + delta[flag]^2*f2/2 +
		delta[flag]^3*f3/6 + delta[flag]^4*f4/24 + delta[flag]^5*f5/120
	out
	
	# sapply( delta, function(delta){
		# b <- c.eta*delta/2
		# I <- tryCatch( digami(b,a), error=function(e) NULL)
		# out <- NA
		# if( !is.null(I) ){
			# out <- c.eta/2*( I[1]/I[6] - a/b )
		# }
		# if( !is.finite(out) | is.null(I) )
			# out <- -c.eta/2*a/(a+1) + delta * (c.eta/2)^2*( a/(a+2) - (a/(a+1))^2 )
		# out	
	# })
		
	# else {
		# I <- ainc.gamma( b[flag], a )
		# out[flag] <- c.eta/2*( I[1,]/I[2,] - a/b[flag] )
	# }
}

# wGdelta2 <- function(delta, eta, p){
# #	c.eta <- (1-eta)^(-1)
# #	a <- p/2 + 1/eta
# #   	-a/delta * exp( pgamma(1, a+1, c.eta*delta/2, log=TRUE) - pgamma(1, a, c.eta*delta/2, log=TRUE) )
	# a <- p/2 + 1/eta
	# c.eta <- 1/(1-eta)
	# #sapply( delta, function(delta){
		
		# b <- c.eta*delta/2  ## b es la x
		
		# I6 <- ainc.gamma( b, a)
		# I1 <- grad( ainc.gamma, b, a=a, side=c(+1,rep(NA,length(delta)-1)) )
		
		# c.eta/2*( I1/I6 - a/b )
	# #})
# }


wGdeltadelta <- function(delta, eta, p){
	a <- p/2 + 1/eta
	c.eta <- 1/(1-eta)
	out <- sapply( delta, function(delta){
		
		b <- c.eta*delta/2
		I <- digami(b,a)
		
		(c.eta/2)^2*( -(I[1]/I[6])^2 + I[2]/I[6] + a/b^2 )
	})
	
	f0 <- (c.eta/2)^2*(a/((a + 1)^2*(a + 2)))
	f1 <- (c.eta/2)^3*(2*a*(a - 1)/((a + 1)^3*(a + 2)*(a + 3)))
	f2 <- (c.eta/2)^4*(6*a*(a^3 - a^2 - 6*a + 2)/((a + 1)^4*(a + 2)^2*(a + 3)*(a + 4)))
	f3 <- (c.eta/2)^5*(24*a*(a - 1)*(a^3 - 3*a^2 - 14*a + 2)/((a + 1)^5*(a + 2)^2*(a + 3)*(a + 4)*(a + 5)))
	f4 <- (c.eta/2)^6*(120*a*(a^7 - 3*a^6 - 49*a^5 - 57*a^4 + 222*a^3 + 264*a^2 - 198*a + 12)/((a + 1)^6*(a + 2)^3*(a + 3)^2*(a + 4)*(a + 5)*(a + 6)))
	f5 <- (c.eta/2)^7*(720*a*(a - 1)*(a^7 - 7*a^6 - 81*a^5 - 37*a^4 + 766*a^3 + 1048*a^2 - 390*a + 12)/((a + 1)^7*(a + 2)^3*(a + 3)^2*(a + 4)*(a + 5)*(a + 6)*(a + 7)))
	f6 <- (c.eta/2)^8*(5040*a*(a^11 - 8*a^10 - 172*a^9 - 354*a^8 + 3265*a^7 + 13498*a^6 + 1164*a^5 - 46836*a^4 - 23650*a^3 + 38356*a^2 - 6096*a + 96)/((a + 1)^8*(a + 2)^4*(a + 3)^2*(a + 4)^2*(a + 5)*(a + 6)*(a + 7)*(a + 8)))
	
	out <- -exp(log(gamma_inc_P(a+1,b)*a)-log(gamma_inc_P(a,b)*delta))
	flag <- !is.finite(out) | delta^6*abs(f6)/720 < 1e-5
	
	out[flag] <- f0 + delta[flag]*f1 + delta[flag]^2*f2/2 +
		delta[flag]^3*f3/6 + delta[flag]^4*f4/24 + delta[flag]^5*f5/120
	out


}



dg.lange <- function(eta, p){
	c.eta <- (1-eta)^(-1)
	a <- p/2 + 1/eta

	list( dg=a*eta^(-1)*( (a+1)*(a+2) + eta^(-1) + 1 )/(eta^(-1) + 1 )/(a+1)^2/(a+2), error = a^(-2))

}

fg.lange <- function(eta, p){
	c.eta <- (1-eta)^(-1)
	a <- p/2 + 1/eta

	list( fg=p*(p+2)*a*( (a+1)*(a+2) + eta^(-1)  )/(a+1)^2/(a+2), error = a^(-2))

}




