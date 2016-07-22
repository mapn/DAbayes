#' @title Compute the logarithm of the pdf of normal random variable with mean u and cariance matrix sigma
#' 
#' @param x normal random variable (vector) 
#'
#' @param u mean vector of x
#' 
#' @param sigma covariance matrix of x
#'
#' @return a real number representing the logarithm of pdf for x ~ N(u, sigma)
#' 
#' @description This function computes the logarithm of the density function for 
#' normal distribution with mean and covariance matrix 
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @export
#' 
#' @keywords models
#' 
#' @examples 
#' set.seed(1234)
#' n <- 10
#' x <- matrix(rnorm(n), nrow=n, ncol=1)
#' u <- matrix(rnorm(n), nrow=n, ncol=1)
#' sigma <- exp(-as.matrix(dist(rnorm(n))))
#' out <- lnnormpdf(x, u, sigma)



lnnormpdf <- function(x,u,sigma){
  
  if(!is.matrix(x)){
  	x <- as.matrix(x, length(x), 1)
  }

  if(!is.matrix(u)){
  	u <- as.matrix(u, length(u), 1)
  }

  if(!is.matrix(sigma)){
  	n <- as.integer(sqrt(length(sigma)))
  	sigma <- as.matrix(sigma, n, n)
  }

    Rsig <- chol(sigma)
    temp <- backsolve(Rsig,x-u,transpose=TRUE)
    out <- -sum(log(diag(Rsig))) - t(temp)%*%temp/2-length(x)/2*log(2*pi)
  
  return(as.vector(out))
}
