#' @title Gibbs sampling with Metropolis-Hastings algorithm for the Bayesian model
#' of detection and attribution problems
#'
#' @param outW a 7 by 1 list containing precomputed quantities associated with W from 
#'       the output of function computeW(...)
#' 
#' @param theta a list of 3 elements contained in parameters in the MCMC algorithm
#' 
#' @param prior_theta a list of 2 elements, each of whcih is the quantity 
#'       in the prior distribution of lambda
#'  
#' @param ADtheta a list of 6 elements, each of which is the parameter in
#'      the corresponding proposal distribution
#'      
#' @param N an integer representing the number of ensemble members
#' 
#' @param Lj an m by 1 vector containing the number of runs for each 
#'        forcing scenario
#' 
#' @param niter an integer specifying the total number of MCMC iterations
#' 
#' @param start_adapting an integer specifiying when to adapt proposal in 
#'        the MCMC algorithm
#'
#' @return a list of 6 elements containing posterior quantities of 
#'         parameters, log-likelihood, chisq statistics, and prior:  
#'         
#' beta: a matrix holds the posterior samples for the parameter beta with each 
#'              row corresponding to each beta
#'              
#' logsigma: a vector holds the posterior samples for the parameter log of sigma
#' 
#' lambda: a vector holds the posterior samples for the parameter lambda      
#'    
#' loglik: a vector holds the log-likelihood evaluated with updated parameters
#' 
#' chisq: a vector holds chisquare statistics for residual consistency test
#' 
#' prior: a vector holds the prior density evaluated with updated parameters
#'
#' @description This function implements Gibbs sampling with Metropolis-Hasting
#' algorithm to sample from posterior distributions for the proposed Bayesian 
#' statistical model. This function is used in parallel computing of MCMC algorithm
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @keywords models
#'
#' @export
#' 
#' @seealso DA_GIBBS
#' 



parDA_GIBBS <- function(outW,theta,prior_theta,ADtheta,N,Lj,niter=20000,start_adapting=500)
{

m <- length(theta$beta)
r <- length(theta$lambda)

Res <- vector(mode="list",6)
names(Res) <- c("beta","logsigma","lambda","loglik","chisq","prior")
Res$beta <- matrix(0,m,niter)
Res$logsigma <-rep(0,niter)
Res$lambda <- matrix(0,r,niter)
Res$loglik <- rep(NA,niter)
Res$chisq <- rep(NA,niter)
Res$prior <- rep(NA,niter)

out <- parDA_loglik(outW,theta,N,Lj)
loglik_old <- out[[1]]
chisq_old <- out[[2]]


for(i in 1:niter)
{
	out <- parDA_theta_MH(outW,theta,prior_theta,ADtheta,
		   loglik_old,chisq_old,N,Lj)
	theta <- out$theta
	loglik_old <- out$loglik
	chisq_old <- out$chisq
	prior <- out$prior
	#if(i > start_adapting){
	#	ADtheta <- DA_ADtheta_update(theta,ADtheta,i-1)
	#}

	# update proposal covariance matrix for beta
	if (i>start_adapting){
		out <- update_var(theta$beta,ADtheta$cov_beta,ADtheta$mean_beta,i-1)
		ADtheta$cov_beta <- out[[1]]
		ADtheta$mean_beta <- out[[2]]
	}

	# update proposal covariance matrix for logsigma
	out <- update_var(theta$logsigma,ADtheta$cov_logsigma,ADtheta$mean_logsigma,i-1)
	ADtheta$cov_logsigma <- out[[1]]
	ADtheta$mean_logsigma <- out[[2]]


	# update proposal covariance matrix for lambda
	out <- update_var(theta$lambda,ADtheta$cov_lambda,ADtheta$mean_lambda,i-1)
	ADtheta$cov_lambda <- out[[1]]
	ADtheta$mean_lambda <- out[[2]]


	  
	Res$beta[ ,i] <- theta$beta
	Res$logsigma[i] <- theta$logsigma
	Res$lambda[ ,i] <- theta$lambda
	Res$loglik[i] <- loglik_old
	Res$chisq[i] <- chisq_old
	Res$prior[i] <- prior
	#Res$ADtheta[[i]] <- ADtheta
	
}


return(Res)

}
