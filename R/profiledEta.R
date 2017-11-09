#' Profiled eta
#' @param y Matrix n by p
#' @export
profiledEta <- function( y ){
	n <- nrow(y)
	p <- ncol(y)
	past <- p*(p+1)/2
	sv <- startValues.slash( y )
	stopifnot( sv$eta0ml < 0.5  )
	Omega <- OmegaMat(sv$Sigma0, sv$eta0ml)
	Info <- infoFisher.slash( y, sv$mu0, sv$Sigma0, sv$eta0ml, block = FALSE )
	Fetaeta <- Info[p + past + 1, p + past + 1]
	Fetagamma <- Info[p + past + 1, 1:(p+past)]

	sigma2eta <- c(1/Fetaeta + 1/Fetaeta^2 * t(Fetagamma) %*% Omega %*% Fetagamma)
	se <- sqrt(sigma2eta / n )

	list( eta=sv$eta0ml, se=se)
}
