Dmatrix <- function(n){
	if (missing(n)) 
        stop("argument n is missing")
    if (!is.numeric(n)) 
        stop("argument n is not numeric")
    if (n != trunc(n)) 
        stop("argument n is not an integer")
    if (n < 2) 
        stop("argument n is less than 2")
	p <- n*(n+1)/2
	nsq <- n*n
	Dt <- matrix(0, nrow = p, ncol = nsq)
	T <- unlist(matrixcalc::T.matrices(n))
	u <- matrixcalc::u.vectors(n)

	out <- .C("Dmatrix",
		Dt = double( n*(n+1)/2 * n*n  ),
		I = as.double(u$I),
		k = as.integer(u$k),
		T = as.double(T),
		n = as.integer(n), PACKAGE="slash"  )
	matrix(out$Dt, nsq, p, byrow=TRUE) # t(Dp)
}
