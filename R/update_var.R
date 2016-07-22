#' @title Adjust quantities in proposal distributions
#' 
#' @param newvec updated parameters in prior distribution
#' 
#' @param oldcov covariance matrix in proposal distribution
#' 
#' @param oldmean mean in proposal distribution
#' 
#' @param oldits index of the previous MCMC iteration 
#' 
#' @return updated quantities in proposal distributions
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @export
#' 
#' @keywords models
#' 
#' @description This function specifies how we adjust mean and covariance matrix
#' in the proposal distribution
#' 

update_var <- function(newvec,oldcov,oldmean,oldits)
{
	newmean <- (oldits*oldmean + newvec) / (oldits+1)
	newcov <- (oldits*oldcov + oldits*(oldmean%*%t(oldmean)) 
			  - (oldits+1)*(newmean%*%t(newmean)) + newvec%*%t(newvec)) / (oldits+1)
	out <- list(newcov,newmean)	
	#names(out) <- c("newcov", "newmean")
	return(out)	  
}

