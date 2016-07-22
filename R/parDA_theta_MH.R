#' @title Update parameters from the model in the MCMC algorithm
#' 
#' @param outW a 7 by 1 list containing precomputed quantities associated with W 
#'
#' @param theta a list of 3 elements contained in parameters in MCMC algorithm
#'
#' @param prior_theta quantities (parameters) in prior dist. 
#'
#' @param ADtheta a 6 by 1 list containing quantities to adjust mean and 
#'        covariance in proposal dist
#'
#' @param loglik_old a real number representing likelihood for given parameters
#' 
#' @param chisq_old a real number representing the quadratic form
#'
#' @param tol a very small value to avoid numerical instability
#'
#' @param N number of temperature ensemble members
#' 
#' @param Lj an m by 1 vector containing the number of runs for each 
#'        forcing scenario
#'
#' @return a list of updated parameters, loglikelihood, chisq statistics (quadratic form),
#' and prior
#' 
#' 
#' @description This function updates parameters for the Bayesian statistical 
#' model. gamma is updated with a discrete uniform prior; beta is updated with
#' noninformative prior 1 and random-walk proposal; log of sigma (variance in W)
#' is updated with normal prior and random-walk proposal; lambda is updated 
#' with normal prior and random-walk proposal. This function is used in 
#' parallel computing of the MCMC algorithm.
#' 
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @keywords models
#' 
#' @export
#' @seealso DA_theta_MH




parDA_theta_MH <- function(outW,theta,prior_theta,ADtheta,loglik_old,chisq_old,N,Lj,tol=1e-10)
{

	m <- length(theta$beta)
	r <- length(theta$lambda)


	# prior information for parameters
	beta_mean <- prior_theta$beta_mean
	beta_var <- prior_theta$beta_var
	lambda_mean <- prior_theta$lambda_mean
	lambda_var <- prior_theta$lambda_var
	logsigma_mean <- prior_theta$logsigma_mean
	logsigma_var <- prior_theta$logsigma_var

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
	Delta.lambda <- cov_lambda*2.38^2/r + tol*diag(rep(1,r),r,r)


	############## update beta ###################

	beta <- theta$beta 
	prop_theta <- theta

	# propose new beta
	beta_prop <- beta + t(chol(Delta.beta))%*%rnorm(m)
	prop_theta$beta <-beta_prop

	# update likelihood
	out <- parDA_loglik(outW,prop_theta,N,Lj)
	loglik_prop <- out$loglik
	chisq_prop <- out$chisq

	# (log) prior ratio
	logpriorratio <- lnnormpdf(beta_prop,beta_mean,beta_var) - lnnormpdf(beta,beta_mean,beta_var)
	#logpriorratio <- 0
	# calculate acceptance prob
	logpropratio <- 0    # symmetric random-walk proposal

	MH_ratio <- loglik_prop-loglik_old + logpriorratio + logpropratio

	if(log(runif(1)) < MH_ratio){
		theta$beta <- prop_theta$beta
		loglik_old <- loglik_prop
		chisq_old <- chisq_prop
	}


	######################################################


	############ update log of sigma (variance in C) ###########
	logsigma <- theta$logsigma
	prop_theta <- theta
	# propose new sigma
	logsigma_prop <- logsigma + sqrt(Delta.logsigma)*rnorm(1)
	prop_theta$logsigma <- logsigma_prop
	# calculate log-likelihood with proposed value
	out <- parDA_loglik(outW,prop_theta,N,Lj)
	loglik_prop <- out$loglik
	chisq_prop <- out$chisq
	# (log) prior ratio
	logpriorratio <- lnnormpdf(logsigma_prop,logsigma_mean,logsigma_var) - lnnormpdf(logsigma,logsigma_mean,logsigma_var)
	# calculate acceptance prob
	logpropratio <- 0 # (symmetric random-walk proposal)

	MH_ratio <- loglik_prop-loglik_old + logpriorratio + logpropratio

	if (log(runif(1)) < MH_ratio){
		theta$logsigma <- prop_theta$logsigma
		loglik_old <- loglik_prop
		chisq_old <- chisq_prop
	}



	############ update lambda (log diagonal elements of K) #########

	lambda <- theta$lambda
	prop_theta <- theta
	# propose new lambda
	lambda_prop <- lambda + t(chol(Delta.lambda))%*%rnorm(r)
	prop_theta$lambda <- lambda_prop
	# calculate log-likelihood with proposed value
	out <- parDA_loglik(outW,prop_theta,N,Lj)
	loglik_prop <- out$loglik
	chisq_prop <- out$chisq
	# (log) prior ratio
	logpriorratio <- lnnormpdf(lambda_prop,lambda_mean,lambda_var) - lnnormpdf(lambda,lambda_mean,lambda_var)
	# calculate acceptance prob
	logpropratio <- 0   # symmetric random-walk proposal
	# MH ratio
	MH_ratio <- loglik_prop-loglik_old + logpriorratio + logpropratio

	if(log(runif(1)) < MH_ratio){
		theta$lambda <- prop_theta$lambda
		loglik_old <- loglik_prop
		chisq_old <- chisq_prop
	}

	prior <- -lnnormpdf(beta,beta_mean,beta_var)*
	         exp(-exp(lnnormpdf(theta$logsigma,logsigma_mean,logsigma_var)))*
	         exp(lnnormpdf(theta$lambda,lambda_mean,lambda_var))

	out <- list(theta,loglik_old,chisq_old,prior)
	names(out) <- c("theta","loglik","chisq","prior")

	return(out)



}
