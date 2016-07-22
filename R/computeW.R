#' @title Obtain quantities associated with W
#' 
#' @param y an n by N matrix of ensembles of measured climate variable
#' (e.g., temperature) increases
#'
#' @param x a list of m elements, each of which is model-output climate
#' variable (e.g., temperature) increases under a specific forcing 
#' scenario
#'
#' @param B an n by r basis function matrix
#'
#' @param neighborhood_matrix an n by n matrix
#'
#' @param gamma_prior a vector of possible discrete values for the
#'        prior distribution of gamma
#' 
#' @return a list of 7 elements
#'
#' @description This function computes quantities associated with the 
#' covariance matrix W, and these quantities are used to evaluate  
#' the likelihood function for the proposed statistical model.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @export
#' 
#' @keywords models
#' 
#' @examples 
#' 
#' 
#' #################### simulate data ########################
#' set.seed(1234)
#' n <- 30 # number of spatial grid cells on the globe
#' N <- 10 # number of ensemble members
#' m <- 3 # number of forcing scenarios
#' Lj <- c(5, 3, 7) # number of runs for each scenario
#' L0 <- 8 # number of control runs without any external forcing scenario
#' trend <- 30
#' DAdata <- simDAdata(n, N, m, Lj, L0, trend)
#' # ensembles of the measured variable
#' y <- DAdata[[1]]
#' # model outputs for the measured variable under different forcing scenarios
#' x <- DAdata[[2]]
#' # model outputs for the measured variable without any external forcing scenario
#' x0 <- DAdata[[3]]
#' #################### end of simulation ####################
#' 
#' 
#' 
#' ###########################################################
#' ### preprocessing 
#' # center the data
#' y <- y - mean(y)
#' for(j in 1:m){
#'   x[[j]] <- x[[j]]-mean(x[[j]])
#' }
#' # construct basis function matrix B with principal components
#' r <- 3
#' empiricalcovmat <- cov(t(x0))
#' 
#' # compute first r eigenvalues and eigenvectors
#' temp <- RSpectra::eigs_sym(empiricalcovmat, r, which="LM")
#' B <- temp$vectors
#' K_hat <- temp$values
#' 
#' 
#' lambda_mean <- log(K_hat)
#' lambda_var <- var(lambda_mean)/3^2
#' 
#' 
#' ### pre-computation for W
#' 
#' # possible values for the discrete prior of gamma
#' gamma_prior <- c(0.5,0.99,1,1.01,2,5)
#' ng <- length(gamma_prior)
#' 
#' # neighborhood matrix (replace by real one)
#' neighborhood_matrix <- Matrix::sparseMatrix(1:n, 1:n, x=rep(1,n))
#' \dontrun{
#' outW <- computeW(y,x,B,neighborhood_matrix,gamma_prior)
#' }



computeW <- function(y,x,B,neighborhood_matrix,gamma_prior)
{

m <- length(x)
r <- dim(B)[2]
ng <- length(gamma_prior)

avg <- function(x) {return(apply(x, 1, mean))}
Xbar <- sapply(x, avg)


logdet_Wti <- rep(0, ng)
trYpWtiY <- rep(0, ng)

XpWtiX <- vector(mode="list", ng)
BpWtiB <- vector(mode="list", ng)
XpWtiYone <- as(matrix(0, nrow=m, ncol=ng), "dgCMatrix")
BpWtiYone <- as(matrix(0, nrow=r, ncol=ng), "dgCMatrix")
#XpWtiYone <- matrix(0, nrow=m, ncol=ng)
#BpWtiYone <- matrix(0, nrow=r, ncol=ng)
BpWtiX <- vector(mode="list", ng)


for(j in 1:ng)
{
	# W_tilde_inverse: a sparse matrix
	Wti <- calc_prec_matrix(neighborhood_matrix, gamma_prior[j])
	RWti <- chol(Wti)
	# log determinant of W_tilde_inverse
	logdet_Wti[j] <- 2*sum(log(diag(RWti)))  ##
	RWtiY <- RWti %*% y
	# trace of y'*W_tilde_inverse*y
	trYpWtiY[j] <- sum(RWtiY*RWtiY)  ## 
	RWtiXbar <- RWti %*% Xbar
	# Xbar'*W_tilde_inverse*Xbar
	XpWtiX[[j]] <- t(RWtiXbar) %*% RWtiXbar  ##
	RWtiB <- RWti %*% B
	BpWtiB[[j]] <- t(RWtiB) %*% RWtiB  ##
	ysum <- apply(y,1,sum)
	# Xbar'*W_tilde_inverse*y*(vector of ones)
	XpWtiYone[ ,j] <- t(Xbar) %*% Wti %*% ysum  ##
	# B'*W_tilde_inverse*y*(vector of ones)
	BpWtiYone[ ,j] <- t(B) %*% Wti %*% ysum  ##
	BpWtiX[[j]] <- as.matrix(t(B) %*% Wti %*% Xbar)  ##

}

out <- list(logdet_Wti, trYpWtiY, XpWtiX, BpWtiB, 
	   as.matrix(XpWtiYone), as.matrix(BpWtiYone), BpWtiX)

return(out)

}
