#' Test Lagrange Multipliers
#' 
#' @param y Numeric vector
#' @param p ARCH order
#' @param q GARCH order
#' 
#' @return A list
#' \describe{
#'   \item{LM}{Statistic}
#'   \item{p.value}{p value}
#'   \item{df}{Degree of freedom}
#' }
#' 
#' @references 
#' Lee, J. 1991. A Lagrange multiplier test for GARCH models. Economics Letters 37: 265-271.
#' 
#' @examples 
#' set.seed(1984)
#' y <- rnorm(500)
#' LMgarch.test(y, 1, 1)
#' # H0: alpha_1=beta_1=0 
#' # is accepted
#' @export
LMgarch.test <- function( y, p, q){
#H0: alpha_1=...=alpha_q=beta_1=beta_p=0
#H1: There exist at least one alpha_i an beta_j >0, for i=1,...,q,   j=1,...,p
fit <- lm( y ~ 1  )
e <- resid(fit)
n <- length(e)

What <- matrix(0, n-q, q+1)
for( t in (q+1):n) What[t-q,] <- c(1, e[t-(1:q)]^2)

m <- max(p,q)
sigma2hat <- sum(e[(m+1):n]^2)/(n-m)

f0 <- e[(q+1):n]^2/sigma2hat - 1


LM <- 1/2*(t(f0)%*%What)%*%solve(t(What)%*%What, t(What)%*%f0)
list( LM=c(LM), p.value=pchisq(c(LM),q,lower=FALSE), df=q)
}




