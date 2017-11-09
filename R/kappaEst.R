#' Multivariate excess of kurtosis
#' 
#' Evaluate multivariate excess of kurtosis according emphirical kappa
#' or expected kappa.
#'
#' @param y Matrix n rows by p columns.
#' @param eta Numeric. Default NULL
#'
#' @return A numeric value
#'
#' @export
kappaEst <- function( y, eta=NULL ){
  if( is.null(eta) ){
    delta2 <- mahalanobis(y, colMeans(y), cov(y))
    kappa <- mean(delta2^2)/(p*(p+2)) - 1
  } else {
    kappa <- eta^2/(1-2*eta)
  }
  kappa
}
