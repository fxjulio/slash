#' Bootstrap samples of profiled eta
#'
#' @param y Observation matrix. n by p. 
#' @param B Number of samples. Default B=1000.
#' @param method String. Parametric (P) or no parametric (NP)
#' @export
bootEta <- function(y, B=1000, method="P"){
  stopifnot(method %in% c("P", "NP"))
  n <- nrow(y)
  
  if( method == "NP")
	  return(replicate(B, {
  		yprime <- y[sample(n, replace = TRUE),]
  		svprime <- startValues.slash(yprime)
  		svprime$eta0
	  })  )
  if( method == "P"){
    sv <- startValues.slash(y)
    return(
      replicate(B,{
        yprime <- sim.slash(n, sv$mu0, sv$Sigma0, sv$eta0ml)
        svprime <- startValues.slash(yprime)
        svprime$eta0ml
      })
    )
  }
}
