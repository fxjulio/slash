Kmatrix <- function (r, c = r) {
    if (missing(r)) 
        stop("argument r is missing")
    if (!is.numeric(r)) 
        stop("argument r is not numeric")
    if (r != trunc(r)) 
        stop("argument r is not an integer")
    if (r < 2) 
        stop("argument r is less than 2")
    if (!is.numeric(c)) 
        stop("argument c is not numeric")
    if (c != trunc(c)) 
        stop("argument c is not an integer")
    if (c < 2) 
        stop("argument c is less than 2")
	
	H <- H.matrices(r, c)
    p <- r * c
	
	#void Kmatrix( double *K, double *H, int *r, int *c )

	out <- .C("Kmatrix",
		K = double( p*p  ),
		H = as.double(unlist(H)),
		r = as.integer(r),
		c = as.integer(c), PACKAGE="slash"  )

	return( matrix(out$K, p, p) )
}
