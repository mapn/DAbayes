#' @title Compute loglikelihood function in the model
#' 
#' @param outW a 7 by 1 list containing precomputed quantities associated with W from 
#'       the output of function computeW(...)
#'
#' @param theta a list of 5 elements containing parameters in MCMC algorithm
#'
#' @param n an integer representing the number of grid cells over the globe
#' 
#' @param N an integer representing the number of ensemble members
#' 
#' @param Lj an m by 1 vector containing the number of runs for each 
#'        forcing scenario
#' 
#' @return a real rumber
#' 
#' @description This function computes the logarithm of the likelihood function 
#' for the Bayesian statistical model.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @keywords models
#' 
#' @export

update_loglik <- function(outW,theta,n,N,Lj)
{

##### Input Arguments:
# gind: an integer containing the index for a vector of discrete gamma priors
# outW: a 7x1 list containing precomputed quantities associated with W from 
#       the output of function computeW(...)
# beta: an mx1 vector containing regression coefficients for m forcing scenarios
# logsigma: log of sigma
# lambda: a rx1 vector containing eigenvalues 
# n: number of grid cells over the globe
# N: number of temperature ensemble members
# Lj: an mx1 vector containg the number of runs for each forcing scenario

gind <- theta$gind
beta <- theta$beta
logsigma <- theta$logsigma
lambda <- theta$lambda



beta <- as.matrix(beta, length(beta), 1)  
logsigma <- as.vector(logsigma)
lambda <- as.matrix(lambda, length(lambda), 1)

# get precomputed quantities associated with W
logdet_Wti <- outW[[1]]
trYpWtiY <- outW[[2]]
XpWtiX <- outW[[3]]
BpWtiB <- outW[[4]]
XpWtiYone <- outW[[5]]
BpWtiYone <- outW[[6]]
BpWtiX <- outW[[7]]

#m <- dim(XpWtiX[[1]])[1]
r <- dim(BpWtiB[[1]])[1]

# other initial computations
K <- diag(exp(as.vector(lambda))) # rxr matrix
sigma <- exp(logsigma)
g_beta <- sqrt(1+sum(as.vector(beta)^2/Lj))  

# log determinnat in loglikelihood fun
temp <- diag(rep(1,r)) + N*g_beta*K%*%BpWtiB[[gind]]/sigma^2
temp <- (temp + t(temp))/2
logdet_SigmaY <- N*(2*n*log(sigma) - logdet_Wti[gind]) + 2*sum(log(diag(chol(temp))))

# compute matrix A
A <- solve(K,diag(rep(1,r)))/(g_beta*N) + BpWtiB[[gind]]/sigma^2
RA <- chol(A)
tmp <- backsolve(RA,t(g_beta*K%*%(BpWtiB[[gind]]/sigma^2)))
Atilde <- g_beta*K - t(backsolve(RA,tmp)) 

# the quadratic form in the exponent
traceYtpWinvYt <- (trYpWtiY[gind] + N*sum((beta%*%t(beta))*XpWtiX[[gind]])
	              - 2*t(beta)%*%XpWtiYone[ ,gind]) / sigma^2
#traceYtpWinvYt <- (trYpWtiY[gind] + N*sum(diag(beta%*%t(beta)%*%XpWtiX[[gind]]))
#	              - 2*t(beta)%*%XpWtiYone[ ,gind]) / sigma^2
BpWinvYtone <- (BpWtiYone[ ,gind] - N*BpWtiX[[gind]]%*%beta) / sigma^2
quadform <- traceYtpWinvYt - t(BpWinvYtone)%*%Atilde%*%BpWinvYtone

# put pieces together
loglikelihood <- -0.5*logdet_SigmaY - 0.5*quadform

return(as.vector(loglikelihood))


}


















