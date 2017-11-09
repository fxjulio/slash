#' @title Simulate multivariate slash reparemetrized
#'
#' @description
#' Return a multivariate slash simulation.
#' 
#' @details
#' Let v distruted beta with parameter 1/eta and 1. And Z distributed multivariate normal
#' with mean zero and  covariance matrix \emph{Sigma}. We obtain a sample slash \emph{Y} with
#' Y = mu + sqrt(1-eta)/sqrt(v) Z
#'
#'
#' @param n Number of simulates
#' @param mu Mean vector with length \code{p}.
#' @param Sigma Covarianze matrix with size \code{p}x\code{p}.
#' @param eta  Parameter between 0 and 1.
#'
#' @return
#' A matrix with \code{n} rows and \code{p} .
#'
#' @examples
#' n <- 1000
#' p <- 4
#' Y <- sim.slash( n, rep(0, p), diag(p), 1/2 )
#' plot(Y, asp=1)
#' 
#' @export
sim.slash <- function(n, mu, Sigma, eta){
	if( length(mu) != nrow(Sigma) ) stop("p dimension is not the same in mu and Sigma")
	if( eta < 0 | eta > 1 ) stop("eta value is not in (0,1)")

	p <- length(mu)
	v <- rbeta(n, 1/eta, 1)
	A <- chol(Sigma, pivot = TRUE)  # A'A = Sigma
	A <- A[,order(attr(A, "pivot"))]
	d <- sqrt(1-eta)
	
	Z <- matrix(rnorm(p*n),n)%*%A
	Y <- matrix(mu, n, p, byrow=TRUE) + d*Z/sqrt(v)

	Y
}