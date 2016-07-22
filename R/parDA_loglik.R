#' @title Compute the loglikelihood function for parallel computing
#' 
#' @param outW a 7 by 1 list containing precomputed quantities from 
#' function parcomputeW(...)
#'
#' @param theta a list of 3 elements containing parameters in the MCMC algorithm
#'
#' @param N an integer representing the number of ensemble members
#'
#' @param Lj an m by 1 vector containing the number of runs for each forcing scenario
#'
#' @return a real number
#'
#' @description This function evaluates the likelihood function for the Bayesian
#' statistical model, and it is used in parallel computing of the MCMC algorithm.
#'
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @keywords models
#' 
#' @export
#' 
#' @seealso update_loglik



parDA_loglik <- function(outW,theta,N,Lj){
### Input Arguments:
# outW: a list of global quantities 
# theta: a list of parameters:
# beta: an mx1 vector containing regression coefficients for m forcing scenarios
# logsigma: log of sigma
# lambda: a rx1 vector containing eigenvalues 
# N: number of temperature ensemble members
# Lj: an mx1 vector containg the number of runs for each forcing scenario


## get parameters
beta <- theta$beta
logsigma <- theta$logsigma
lambda <- theta$lambda

beta <- as.matrix(beta, length(beta), 1)  
logsigma <- as.vector(logsigma)
lambda <- as.matrix(lambda, length(lambda), 1)

r <- length(lambda)

## get other quantities
w <- outW$w
ybar <- outW$ybar
xbar <- outW$xbar
trYpWiY <- outW$trYpWiY
XpWiX <- outW$XpWiX
XpWiybar <- outW$XpWiybar
B <- outW$B

## initial computations
K <- diag(exp(as.vector(lambda)),r,r)
K_chol <- diag(sqrt(exp(as.vector(lambda))),r,r)
sigma <- exp(logsigma)
g_beta <- 1 + sum(as.vector(beta)^2/Lj)
wt <- w + N*g_beta*sigma^2
BpWtiB <- t(B)%*%diag(wt^(-1))%*%B

# log-determinant
temp <- diag(rep(1,r),r,r) + N*g_beta*K_chol%*%BpWtiB%*%t(K_chol)
temp <- (temp + t(temp))/2
logdet_SigmaY <- (N-1)*sum(log(w)) + sum(log(wt)) + 2*sum(log(diag(chol(temp))))

# the matrix A
Kinv <- diag(exp(as.vector(lambda))^(-1),r,r)
temp <- diag(wt^(-1))%*%B%*%solve(Kinv/(g_beta*N)+BpWtiB, diag(rep(1,r),r,r))
A <- K - K%*%t(B)%*%temp

# summaries involving Y
Ytone <- N*(ybar-xbar%*%beta)
v <- Ytone / w
vt <- Ytone / wt
traceYtpWinvYt <- as.vector(trYpWiY + N*t(beta)%*%XpWiX%*%beta - 
                  N*2*t(beta)%*%XpWiybar)

# the quadratic form in the exponent
yJGy <- as.vector(g_beta*(t(v)%*%B%*%A%*%t(B)%*%vt + 
        sigma^2*t(v)%*%vt - sigma^2*t(v)%*%temp%*%t(B)%*%vt))
quadrformexp <- traceYtpWinvYt - yJGy

# put the pieces together
loglikelihood <- -0.5*logdet_SigmaY - 0.5*quadrformexp

out <- list(loglikelihood,quadrformexp)
names(out) <- c("loglik","chisq")
return(out)


}








