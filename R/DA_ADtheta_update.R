# Author: Pulong Ma University of Cincinnati
# Date: May 18, 2016
# Last Modified by: Pulong Ma
# Last Modified Date: May 18, 2016
# Last Modified time: 21:04:11

#' @title Adaptively adjust mean and covariance matrix in the
#'  proposal distributions
#' 
#' @param theta a list of 3 continuous parameters updated in the MCMC algorithm
#' 
#' @param ADtheta a list of 6 parameters in the corresponding proposal
#'      distributions for parameters updated in the MCMC algorithm
#'      
#' @param oldits an integer of the index of the MCMC iteration
#' 
#' @return a list of updated parameters in proposal distributions
#' 
#' @description The function modifies the parameters in proposal distributions for 
#' corresponding parameters at each iteration to obtain better mixing propertities
#' in the MCMC algorithm.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @export
#' 
#' @keywords models
#'





DA_ADtheta_update <- function(theta,ADtheta,oldits)
{

	# get parameters
	beta <- theta$beta
	logsigma <- theta$logsigma
	lambda <- theta$lambda

	# get quantities for adaptive MCMC
	mean_beta <- ADtheta$mean_beta
	cov_beta <- ADtheta$cov_beta
	mean_logsigma <- ADtheta$mean_logsigma
	cov_logsigma <- ADtheta$cov_logsigma
	mean_lambda <- ADtheta$mean_lambda
	cov_lambda <- ADtheta$cov_lambda

	# update proposal covariance matrix for beta
	out <- update_var(beta,cov_beta,mean_beta,oldits)
	cov_beta <- out[[1]]
	mean_beta <- out[[2]]

	# update proposal covariance matrix for logsigma
	out <- update_var(logsigma,cov_logsigma,mean_logsigma,oldits)
	cov_logsigma <- out[[1]]
	mean_logsigma <- out[[2]]


	# update proposal covariance matrix for lambda
	out <- update_var(lambda,cov_lambda,mean_lambda,oldits)
	cov_lambda <- out[[1]]
	mean_lambda <- out[[2]]


	ADtheta$mean_beta <- mean_beta 
	ADtheta$cov_beta <- cov_beta 
	ADtheta$mean_logsigma <- mean_logsigma 
	ADtheta$cov_logsigma <- cov_logsigma 
	ADtheta$mean_lambda <- mean_lambda 
	ADtheta$cov_lambda <- cov_lambda 

	return(ADtheta)





}
