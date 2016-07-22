#' @title Update parameters in the MCMC algorithm
#' 
#' @param outW a 7 by 1 list containing precomputed quantities  
#'       associated with W from the output of function computeW(...)
#'
#' @param theta a list of 4 elements containing parameters in the MCMC algorithm
#'
#'
#' @param prior_theta parameters in prior distributions
#'
#' @param ADtheta a 6 by 1 list containing quantities to adjust mean and 
#'        covariance in proposal distributions
#'        
#' @param tol a very small value to avoid numerical instability
#'
#' @param n number of grid cells over the globe
#' 
#' @param N number of ensemble members
#' 
#' @param Lj an m by 1 vector containing the number of runs for each 
#'        forcing scenario
#'
#' @return a list of 4 elements containing updated parameters in the MCMC
#' iteration
#' 
#' @description This function updates parameters (gamma, beta, sigma, lambda)  
#' for the Bayesian statistical model. gamma is updated with a 
#' discrete uniform prior; beta is updated with the noninformative prior 1
#'  and random-walk proposal; log of sigma (variance in W)
#' is updated with the normal prior and random-walk proposal; 
#' lambda is updated with the normal prior and random-walk proposal.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @keywords models
#' 
#' @seealso DA_GIBBS function

DA_theta_MH <- function(outW,theta,prior_theta,ADtheta,n,N,Lj,tol=1e-10)
{

	#m <- length(theta$beta)
	m <- length(theta[['beta']])
	#r <- length(theta$lambda)
	r <- length(theta[["lambda"]])
	ng <- length(outW[[1]])


	# prior information for parameters
	lambda_mean <- prior_theta$lambda_mean
	lambda_var <- prior_theta$lambda_var

	# quantities for adaptive MCMC
	mean_beta <- ADtheta$mean_beta
	cov_beta <- ADtheta$cov_beta
	mean_logsigma <- ADtheta$mean_logsigma
	cov_logsigma <- ADtheta$cov_logsigma
	mean_lambda <- ADtheta$mean_lambda
	cov_lambda <- ADtheta$cov_lambda

	# compute variance for proposal dist.
	Delta.beta <- cov_beta*2.38^2/m + tol*diag(rep(1,m))
	Delta.logsigma <- cov_logsigma*2.38^2/length(cov_logsigma) + 
	    			  tol*diag(rep(1,length(cov_logsigma)))
	Delta.lambda <- cov_lambda*2.38^2/r + tol*diag(rep(1,r))


	################ update gamma ################

	# calculate log-likelihood for each possible gamma value
	loglik_gamma <- rep(0,ng)
	temp_theta <- theta
	
	for(j in 1:ng)
	{
	  temp_theta$gind <- j
		loglik_gamma[j] <- update_loglik(outW,temp_theta,n,N,Lj)
	}
	#sample from the (discrete) posterior dist assuming a uniform prior
	gind <- sample(1:ng,1,replace=TRUE,prob=exp(loglik_gamma))
	loglik_old <- loglik_gamma[gind]

	theta$gind <- gind

	##############################################

	############## update beta ###################

	beta <- theta$beta 
	prop_theta <- theta

	# propose new beta
	beta_prop <- beta + t(chol(Delta.beta))%*%rnorm(m)
	prop_theta$beta <-beta_prop

	# update likelihood
	loglik_prop <- update_loglik(outW,prop_theta,n,N,Lj)
	# (log) prior ratio
	logpriorratio <- 0   # assuming improper uniform prior on R^m
	# calculate acceptance prob
	logpropratio <- 0    # symmetric random-walk proposal

	MH_ratio <- loglik_prop-loglik_old + logpriorratio + logpropratio

	if(log(runif(1)) < MH_ratio)
	{
		theta$beta <- prop_theta$beta
		loglik_old <- loglik_prop
	}


	######################################################


	############ update log of sigma (variance in W) ###########
	logsigma <- theta$logsigma
	prop_theta <- theta
	# propose new sigma
	logsigma_prop <- logsigma + sqrt(Delta.logsigma)*rnorm(1)
	prop_theta$logsigma <- logsigma_prop
	# calculate log-likelihood with proposed value
	loglik_prop <- update_loglik(outW,prop_theta,n,N,Lj)
	# (log) prior ratio
	logpriorratio <-2*logsigma - 2*logsigma_prop
	# calculate acceptance prob
	logpropratio <- 0 # (symmetric random-walk proposal)

	MH_ratio <- loglik_prop-loglik_old + logpriorratio + logpropratio

	if (log(runif(1)) < MH_ratio)
	{
		theta$logsigma <- prop_theta$logsigma
		loglik_old <- loglik_prop
	}



	############ update lambda (log diagonal elements of K) #########

	lambda <- theta$lambda
	prop_theta <- theta
	# propose new lambda
	lambda_prop <- lambda + t(chol(Delta.lambda))%*%rnorm(r)
	prop_theta$lambda <- lambda_prop
	# calculate log-likelihood with proposed value
	loglik_prop <- update_loglik(outW,prop_theta,n,N,Lj)
	# (log) prior ratio
	logpriorratio <- (sum((lambda-lambda_mean)^2) - 
					sum((lambda_prop-lambda_mean)^2)) / (2*lambda_var)
	# calculate acceptance prob
	logpropratio <- 0   # symmetric random-walk proposal
	# MH ratio
	MH_ratio <- loglik_prop-loglik_old + logpriorratio + logpropratio

	if(log(runif(1)) < MH_ratio)
	{
		theta$lambda <- prop_theta$lambda
		loglik_old <- loglik_prop
	}


	return(theta)



}
