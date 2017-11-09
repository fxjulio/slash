#' Compute de I(x,p) derivates
#' @export
digami <- function(x,p){
	gplog <- lgamma(p)
	psip <- digamma(p)
	psidp <- trigamma(p)
	out <- .Fortran("DIGAMI", 
		d=double(6), 
		x=as.double(x),
		p=as.double(p),
		gplog=as.double(gplog),
		gp1log=as.double(log(p)+gplog),
		psip=as.double(psip),
		psip1=as.double(1/p + psip),
		psidp=as.double(psidp),
		psidp1=as.double(psidp - 1/p^2),
		ifault=as.integer(0) )
	if( out$ifault  ) stop("Error en digami")
	
	as.double(out$d)
}