#' @title Compute n by n covariance matrix W over n grid cells
#'
#' @param neighborhood_matrix an n by n neighborhood matrix for n grid cells
#'
#' @param gamma a real valued parameter in covariance matrix W
#'
#' @return an n by n covariance  matrix
#'
#' @description This function computes the covariance matrix of the ensemble 
#' members for the measured variable, which is a fixed diagonal 
#' matrix with diagonal elements equal to the empirical variance of the 
#' corresponding ensemble memebers. Further work will allow W to be a more 
#' general matrix (e.g., based on a GMRF) with unknown parameters
#'
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @export
#' 
#' @keywords models
#' 
#' @examples 
#' n <- 30  # number of grid cells over the globe
#' neighborhood_matrix <- Matrix::sparseMatrix(1:n, 1:n, x=rep(1,n))
#' gamma <- 0.1
#' W <- calc_prec_matrix(neighborhood_matrix, gamma)
#' Matrix::image(W, xlab="", ylab="", cex=0.5, cex.lab=0.6)




calc_prec_matrix <- function(neighborhood_matrix, gamma)
{
	prec_matrix <- gamma*neighborhood_matrix

	return(prec_matrix)
}

