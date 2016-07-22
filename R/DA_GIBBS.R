
#' @title Gibbs sampling with Metropolis-Hastings algorithm for the Bayesian model
#' of detection and attribution problems
#'
#' @param outW a 7 by 1 list containing precomputed quantities associated with W from 
#'       the output of function computeW(...)
#' 
#' @param theta a list of 4 elements containing parameters in MCMC algorithm
#' 
#' @param prior_theta a list of 2 elements, each of whcih is the quantity 
#'       in the prior distribution of lambda
#'  
#' @param ADtheta a list of 6 elements, each of which contains parameters
#'       in the corresponding proposal distribution
#'      
#' @param n number of grid cells over the globe
#' 
#' @param N number of ensemble members
#' 
#' @param Lj an m by 1 vector containing the number of runs for each 
#'        forcing scenario
#' 
#' @param niter an integer specifying the total number of MCMC iterations
#' 
#' @return a list of 5 elements containing posterior quantities of 
#'         parameters and log-likelihood:
#'         
#'         gind: a vector holds the posterior samples for the parameter gamma
#'         
#'         beta: a matrix holds the posterior samples for the parameter beta with each 
#'              row corresponding to each beta
#'              
#'         logsigma: a vector holds the posterior samples for the parameter log of sigma
#'         
#'         lambda: a vector holds the posterior samples for the parameter lambda
#'         
#'         loglik: a vector holds the log-likelihood evaluated with updated parameters
#' 
#' @description This function implements Gibbs sampling with Metropolis-Hasting
#' algorithm to sample from posterior distributions for the proposed Bayesian 
#' statistical model.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#'
#' @export
#' 
#' @keywords models
#' 
#' @examples 
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
#' outW <- computeW(y,x,B,neighborhood_matrix,gamma_prior)
#' 
#' # MCMC output quantities and initial values of parameters
#' beta <- rep(0,m)
#' lambda <- rep(0,r)
#' logsigma <- 0
#' gamma_ind <- 0
#' 
#' # parameters and initial values for MCMC proposals
#' mean_beta <- beta
#' cov_beta <- 0.1*diag(rep(1, m))
#' mean_logsigma <- logsigma
#' cov_logsigma <- 0.1
#' mean_lambda <- lambda
#' cov_lambda <- lambda_var*diag(rep(1, r))
#' 
#' 
#' 
#' theta <- list(gamma_ind,beta,logsigma,lambda)
#' names(theta) <- c("gind","beta","logsigma","lambda")
#' prior_theta <- list(lambda_mean,lambda_var)
#' names(prior_theta) <- c("lambda_mean","lambda_var")
#' ADtheta <- list(mean_beta,cov_beta,mean_logsigma,cov_logsigma, 
#'                mean_lambda,cov_lambda)
#' names(ADtheta) <- c("mean_beta","cov_beta","mean_logsigma",
#'                     "cov_logsigma","mean_lambda","cov_lambda")
#' 
#' niter <- 10000
#' 
#' \dontrun{
#' ## run GIBBS sampling with Metroplis Hastings algorithm
#' Res <- DA_GIBBS(outW,theta,prior_theta,ADtheta,n,N,Lj,niter)
#' ## check convergence for posterior samples
#' plot(Res$gind,xlab="Iteration")
#' plot(Res$beta[1,],xlab="Iteration", type="l")
#' plot(Res$beta[2,],xlab="Iteration", type="l")
#' plot(Res$beta[3,],xlab="Iteration", type="l")
#' plot(Res$logsigma,xlab="Iteration", type="l")
#' plot(Res$lambda[1,],xlab="Iteration", type="l")
#' plot(Res$lambda[2,],xlab="Iteration", type="l")
#' plot(Res$lambda[3,],xlab="Iteration", type="l")
#' }
#' 



DA_GIBBS <- function(outW,theta,prior_theta,ADtheta,n,N,Lj,niter=20000)
{

#m <- length(theta$beta)
#r <- length(theta$lambda)
m <- length(theta[["beta"]])
r <- length(theta[["lambda"]])

Res <- vector(mode="list",4)
names(Res) <- c("gind","beta","logsigma","lambda")
Res$gind <- rep(0,niter)
Res$beta <- matrix(0,m,niter)
Res$logsigma <-rep(0,niter)
Res$lambda <- matrix(0,r,niter)
Res[[5]] <- rep(NA,niter)
names(Res) <- c("gind","beta","logsigma","lambda", "loglik")


for(i in 1:niter)
{
	theta <- DA_theta_MH(outW,theta,prior_theta,ADtheta,n,N,Lj)
	ADtheta <- DA_ADtheta_update(theta,ADtheta,i)
	
	#print(paste("iteration:", i), quote=FALSE)
    loglik <- update_loglik(outW,theta,n,N,Lj)
  

	Res$gind[i] <- theta$gind
	Res$beta[ ,i] <- theta$beta
	Res$logsigma[i] <- theta$logsigma
	Res$lambda[ ,i] <- theta$lambda
	Res$loglik[i] <- loglik
	
}


return(Res)

}
